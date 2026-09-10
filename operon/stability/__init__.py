"""
Operon Calculator -- Module D: mRNA stability (spec section 9).
**Status: features implemented; the final fitted rate is deliberately
gated.** Spec section 9.2. RBSForge explicitly does not include this --
it is entirely this module's responsibility. Consumes coupled TIRs from
``operon.coupling`` as ribosome-protection inputs, but is otherwise
independent.

Implements the 2021 biophysical model's two physically-grounded
features -- ribosome-protection distance (``Δx``) and RNase-accessibility
(``rnase_score``) -- in full. The final linear/exponential combination
into ``log(mRNA)`` / ``k_decay`` / ``t_half`` needs ``{a, b}`` coefficients
**fit on real qPCR data** (Cetnar & Salis 2021's own calibration set);
this module does not invent them. ``mrna_stability`` always returns the
features; it only returns a rate if the caller supplies real
coefficients explicitly via ``coefficients=``.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

from rbsforge.thermo.fold import BuiltinFolder

from operon.coupling import to_physical
from operon.core import Assembled, OperonHost

# Directionally-motivated (not fit) unpaired-position weights: A/U are the
# 9.4-fold destabilizing motif, G/C is protective (spec section 9.1/9.2).
# Same honesty tier as operon.elongation's placeholder codon weights.
_RNASE_MOTIF_WEIGHT = {"A": 1.0, "U": 1.0, "T": 1.0, "G": 0.3, "C": 0.3}

_FIVE_PRIME_HAIRPIN_SEARCH_NT = 20  # how far into the 5' UTR to look for a protective hairpin
_FIVE_PRIME_HAIRPIN_MIN_STEM = 4


@dataclass
class StabilityFeatures:
    dx: float  # ribosome-protection distance, nt of naked mRNA between elongating ribosomes
    rnase_score: float  # unpaired-AU-richness score over the scored 5' UTR window
    has_5p_hairpin: bool


@dataclass
class StabilityResult:
    features: StabilityFeatures
    per_cistron_dx: List[float]
    log_mrna_level: Optional[float] = None  # section 9.2's log(mRNA) level estimate, not a rate
    k_decay: Optional[float] = None
    t_half: Optional[float] = None
    coefficients_used: Optional[Dict[str, float]] = None


def ribosome_protection_dx(tir_au: float, host: OperonHost) -> float:
    """``dx = r_elong / r_init - L_footprint`` (section 9.2): average nt
    of naked (ribosome-free) mRNA between consecutively initiating
    ribosomes on one CDS. Smaller (even negative, meaning ribosomes queue)
    means more protection, hence slower decay.

    ``tir_au`` is converted via ``operon.coupling.to_physical`` -- the
    same documented identity seam coupling uses, for the same reason (no
    published au -> s^-1 formula exists to invent one here either).
    """
    r_init = to_physical(tir_au, host)
    if r_init <= 0:
        return float("inf")
    return host.r_elong_nominal_nt_s / r_init - host.ribosome_footprint_nt


def _unpaired_positions(window: str) -> List[bool]:
    fold = BuiltinFolder().fold(window)
    paired = set()
    for i, j in fold.pairs:
        paired.add(i)
        paired.add(j)
    return [i not in paired for i in range(len(window))]


def rnase_accessibility_score(utr_and_cds_start: str) -> float:
    """Sum of ``motif_w[nt] for each unpaired position`` over the folded
    5' UTR (+ ~30 nt of CDS1, per section 9.2's window) -- long unpaired
    AU tracts score high (destabilizing); folded/GC-rich regions score
    low (protective).
    """
    seq = utr_and_cds_start.upper().replace("T", "U")
    unpaired = _unpaired_positions(seq)
    return sum(_RNASE_MOTIF_WEIGHT.get(nt, 0.0) for nt, is_unpaired in zip(seq, unpaired) if is_unpaired)


def has_5p_hairpin(utr: str, min_stem: int = _FIVE_PRIME_HAIRPIN_MIN_STEM) -> bool:
    """Whether the 5' UTR folds into a hairpin of at least ``min_stem`` bp
    whose 5'-most paired base is within the first
    ``_FIVE_PRIME_HAIRPIN_SEARCH_NT`` nt -- section 9.4's "small 5' UTR
    hairpin" stabilizer / section 9.3's "stable RppH 5' end" proxy.
    """
    seq = utr.upper().replace("T", "U")
    window = seq[:_FIVE_PRIME_HAIRPIN_SEARCH_NT]
    fold = BuiltinFolder().fold(window)
    if not fold.pairs:
        return False

    pairs_by_i = dict(fold.pairs)
    # An outermost stem pair has no (i-1, j+1) pair immediately enclosing it.
    outermost = [(i, j) for i, j in fold.pairs if pairs_by_i.get(i - 1) != j + 1]
    for i, j in outermost:
        stem = 0
        while pairs_by_i.get(i + stem) == j - stem:
            stem += 1
        if stem >= min_stem:
            return True
    return False


def mrna_stability(
    assembled: Assembled,
    tirs: List[float],
    host: OperonHost,
    coefficients: Optional[Dict[str, float]] = None,
) -> StabilityResult:
    """Ribosome-protection + RNase-accessibility features (section 9.2)
    for an assembled operon, aggregated across cistrons by minimum
    protection (the least-translated cistron is the RNase entry point --
    section 9.2's stated preference over a length-weighted average).

    Returns only the features (``k_decay``/``t_half`` are ``None``)
    unless ``coefficients`` is supplied: a dict with keys ``a0, a1, a2,
    a3, k0, b1, b2`` (section 9.2's formula) from a real fit against
    qPCR data. This module does not ship a default -- see the module
    docstring.
    """
    per_cistron_dx = [ribosome_protection_dx(tir, host) for tir in tirs]
    dx = min(per_cistron_dx) if per_cistron_dx else float("inf")

    utr_end = assembled.starts[0] if assembled.starts else 0
    scored_window_end = min(len(assembled.mrna), utr_end + 30)
    utr = assembled.mrna[:utr_end]
    scored_window = assembled.mrna[:scored_window_end]

    rnase_score = rnase_accessibility_score(scored_window)
    hairpin = has_5p_hairpin(utr)

    features = StabilityFeatures(dx=dx, rnase_score=rnase_score, has_5p_hairpin=hairpin)
    result = StabilityResult(features=features, per_cistron_dx=per_cistron_dx)

    if coefficients is not None:
        a0, a1, a2, a3 = coefficients["a0"], coefficients["a1"], coefficients["a2"], coefficients["a3"]
        k0, b1, b2 = coefficients["k0"], coefficients["b1"], coefficients["b2"]
        log_mrna = a0 + a1 * (-rnase_score) + a2 * (-dx) + a3 * (1.0 if hairpin else 0.0)
        k_decay = k0 * math.exp(b1 * rnase_score + b2 * dx)
        result.log_mrna_level = log_mrna
        result.k_decay = k_decay
        result.t_half = math.log(2) / k_decay if k_decay > 0 else float("inf")
        result.coefficients_used = dict(coefficients)

    return result
