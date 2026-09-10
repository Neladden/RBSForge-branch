"""
Operon Calculator -- Module C: translation elongation rate (TER) and
synonymous codon recoding (spec section 8). **Status: implemented**
(with a clearly-labeled placeholder *E. coli* codon-weight table --
see ``codon_tables.py``; swap in real data before trusting absolute
numbers).

Independent of coupling; operates on one CDS at a time. Used as a score
in Predict mode and as a mutator in Design mode.
"""

from typing import Dict

from .codon_tables import ECOLI_GENOMIC, ECOLI_HIGHLY_TRANSLATED
from .genetic_code import GENETIC_CODE, SYNONYMS, translate

__all__ = [
    "recode_cds",
    "mutate_codon",
    "ter_nt_s",
    "translate",
    "ECOLI_HIGHLY_TRANSLATED",
    "ECOLI_GENOMIC",
]


def _shuffled(items, rng):
    """Fisher-Yates shuffle using only ``rng.random()`` (in [0, 1)), so any
    ``rng`` satisfying this module's minimal interface works, not just a
    full ``random.Random``."""
    out = list(items)
    for i in range(len(out) - 1, 0, -1):
        j = int(rng.random() * (i + 1))
        out[i], out[j] = out[j], out[i]
    return out


def _weighted_choice(options, weights, rng):
    """Roulette-wheel selection using only ``rng.random()`` -- the one
    method every ``rng`` this module accepts must provide."""
    total = sum(weights)
    if total <= 0:
        return options[int(rng.random() * len(options))]
    threshold = rng.random() * total
    cumulative = 0.0
    for option, weight in zip(options, weights):
        cumulative += weight
        if cumulative >= threshold:
            return option
    return options[-1]


def recode_cds(aa: str, table: Dict[str, Dict[str, float]], rng) -> str:
    """Recode a protein sequence to RNA, sampling each position's codon
    from ``table[amino_acid]`` (a ``{codon: weight}`` dict) weighted by
    ``w[c]`` (section 8.2). ``rng`` is anything exposing ``.random()``
    returning a float in [0, 1) (a ``random.Random`` instance satisfies
    this) -- pass a seeded one for reproducible recoding.

    Deterministic (argmax) recoding -- always the highest-weight codon --
    is what "Highly Translated" reduces to when every non-preferred
    synonym shares one placeholder weight (see ``codon_tables.py``); for
    an explicit deterministic pick regardless of table shape, use a
    zero-temperature ``rng`` (one whose ``.random()`` always returns 0.0)
    rather than adding a separate code path here.
    """
    codons = []
    for residue in aa.upper():
        weights_by_codon = table.get(residue)
        if not weights_by_codon:
            raise ValueError(f"no codon weights for amino acid {residue!r}")
        options = list(weights_by_codon.keys())
        weights = list(weights_by_codon.values())
        codons.append(_weighted_choice(options, weights, rng))
    return "".join(codons)


def mutate_codon(cds_nt: str, table: Dict[str, Dict[str, float]], rng) -> str:
    """Design-mode GA neighborhood operator (section 8.3, 18.4): pick one
    random codon position whose amino acid has synonyms, swap it to
    another synonym biased by ``table``, preserving the protein sequence
    exactly. A single-position mutator, not a whole-CDS rewrite -- call
    this repeatedly from the outer search, don't loop it internally.
    """
    seq = cds_nt.upper().replace("T", "U")
    if len(seq) % 3 != 0:
        raise ValueError(f"CDS length {len(seq)} is not a multiple of 3")

    n_codons = len(seq) // 3
    positions = _shuffled(list(range(n_codons)), rng)
    for pos in positions:
        codon = seq[pos * 3 : pos * 3 + 3]
        aa = GENETIC_CODE.get(codon)
        if aa is None or aa == "*":
            continue
        synonyms = SYNONYMS.get(aa, [])
        if len(synonyms) <= 1:
            continue
        weights_by_codon = table.get(aa, {c: 1.0 for c in synonyms})
        options = [c for c in synonyms if c != codon]
        weights = [weights_by_codon.get(c, 1.0) for c in options]
        new_codon = _weighted_choice(options, weights, rng)
        return seq[: pos * 3] + new_codon + seq[pos * 3 + 3 :]
    # every codon was a stop or synonym-free (e.g. Met/Trp only) -- no legal move
    return seq


def ter_nt_s(cds_nt: str, host) -> float:
    """Average translation elongation rate in nt/s, as a harmonic mean
    over per-codon elongation rates (section 8.1) -- one slow codon
    dominates dwell time, so use the harmonic mean, not the arithmetic
    one::

        r_elong_codon[c] = r_max * w[c] / w_max          # codons/s
        TER_cds = 3 / mean(1/r_elong_codon[c] for c in cds)   # nt/s

    ``r_max`` is derived from ``host.r_elong_nominal_nt_s`` (~60 nt/s for
    *E. coli*) divided by 3 nt/codon to get codons/s, matching the
    spec's "~20 aa/s ~= 60 nt/s" figure. ``w`` comes from
    ``host.codon_weights_balanced`` if set, else
    ``host.codon_weights_high``, else the placeholder
    ``ECOLI_GENOMIC`` table.
    """
    table = getattr(host, "codon_weights_balanced", None) or getattr(host, "codon_weights_high", None) or ECOLI_GENOMIC
    r_max_codons_per_s = getattr(host, "r_elong_nominal_nt_s", 60.0) / 3.0
    w_max = max(w for weights in table.values() for w in weights.values())

    seq = cds_nt.upper().replace("T", "U")
    if len(seq) % 3 != 0:
        raise ValueError(f"CDS length {len(seq)} is not a multiple of 3")

    inverse_rates = []
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        aa = GENETIC_CODE.get(codon)
        if aa is None or aa == "*":
            continue
        weight = table.get(aa, {}).get(codon)
        if weight is None:
            raise ValueError(f"no codon weight for {codon!r} ({aa}) in the supplied table")
        rate_codons_per_s = r_max_codons_per_s * weight / w_max
        inverse_rates.append(1.0 / rate_codons_per_s)

    if not inverse_rates:
        raise ValueError("no scorable (non-stop) codons in cds_nt")

    return 3.0 * len(inverse_rates) / sum(inverse_rates)
