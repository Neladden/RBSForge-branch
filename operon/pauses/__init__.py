"""
Operon Calculator -- Module G: ribosomal pause sites (spec section 12).
**Status: implemented** for three of the four documented signals
(internal SD/anti-SD, slow-codon runs, polyproline/stall motifs); stable
mRNA hairpins in the ribosomal E/P site window are explicitly skipped
(the spec itself calls this one "optional, weaker evidence").

Independent motif + structure scanner inside CDSs. The lab never
published a standalone Pause Calculator; this is the operational
definition the spec assembles from adjacent tools.
"""

from dataclasses import dataclass
from typing import List, Optional

from rbsforge.thermo.duplex import best_hybridization
from rbsforge.constants import celsius_to_kelvin

from operon.elongation import ECOLI_HIGHLY_TRANSLATED
from operon.elongation.genetic_code import GENETIC_CODE, translate

# Internal-SD scan window: long enough to hold a full anti-SD-length
# duplex plus a little slack, short enough that best_hybridization's own
# O(window^2) alignment search stays cheap over a whole CDS.
_INTERNAL_SD_WINDOW_SLACK = 6
_INTERNAL_SD_MIN_DUPLEX_LENGTH = 5  # a shorter/weaker minimum than RBS-strength SD search (4) would over-report
_INTERNAL_SD_MAX_DELTA_G = -1.0  # kcal/mol; only a net-favorable (and not marginally so) duplex is a real pause candidate

_SLOW_CODON_WEIGHT_THRESHOLD = 1.0  # strictly below "preferred" in the active codon table
_SLOW_CODON_RUN_MIN = 2  # spec section 12: "a run of >=2-3 is a pause"

_POLYPROLINE_RUN_MIN = 2  # "PP" at minimum; PPP/PPG are the named examples, not the only ones counted


@dataclass
class PauseHit:
    kind: str  # "internal_sd" | "slow_codon_run" | "polyproline"
    start: int  # 0-based nt offset into cds_nt
    end: int  # exclusive
    score: float  # higher = stronger pause signal; meaning differs by kind, see notes
    notes: str = ""


def _scan_internal_sd(cds_nt: str, host, temperature_c: Optional[float] = None) -> List[PauseHit]:
    asd = host.pack.asd
    temperature_k = celsius_to_kelvin(temperature_c if temperature_c is not None else host.pack.t_growth_c)
    window_len = len(asd) + _INTERNAL_SD_MIN_DUPLEX_LENGTH + _INTERNAL_SD_WINDOW_SLACK
    seq = cds_nt.upper().replace("T", "U")
    n = len(seq)
    hits: List[PauseHit] = []

    def no_spacing_bias(_aligned_spacing: int) -> float:
        return 0.0

    i = 0
    while i < n:
        window = seq[i : i + window_len]
        if len(window) < _INTERNAL_SD_MIN_DUPLEX_LENGTH:
            break
        duplex, _total = best_hybridization(
            window,
            asd,
            start_codon_index=len(window),
            spacing_energy_fn=no_spacing_bias,
            temperature_k=temperature_k,
            min_length=_INTERNAL_SD_MIN_DUPLEX_LENGTH,
        )
        if duplex.length >= _INTERNAL_SD_MIN_DUPLEX_LENGTH and duplex.delta_g_hybrid <= _INTERNAL_SD_MAX_DELTA_G:
            start = i + duplex.mrna_start
            end = i + duplex.mrna_end
            hits.append(
                PauseHit(
                    kind="internal_sd",
                    start=start,
                    end=end,
                    score=-duplex.delta_g_hybrid,  # more negative hybrid dG -> larger (stronger) score
                    notes=f"SD-like '{duplex.mrna_segment}' vs anti-SD '{duplex.antisd_segment}'",
                )
            )
            i = end  # don't re-report the same duplex from overlapping windows
        else:
            i += 1
    return hits


def _scan_slow_codon_runs(cds_nt: str, host) -> List[PauseHit]:
    table = getattr(host, "codon_weights_balanced", None) or getattr(host, "codon_weights_high", None) or ECOLI_HIGHLY_TRANSLATED
    seq = cds_nt.upper().replace("T", "U")
    n_codons = len(seq) // 3
    hits: List[PauseHit] = []
    run_start = None
    for pos in range(n_codons):
        codon = seq[pos * 3 : pos * 3 + 3]
        aa = GENETIC_CODE.get(codon)
        weight = table.get(aa, {}).get(codon) if aa else None
        is_slow = weight is not None and weight < _SLOW_CODON_WEIGHT_THRESHOLD
        if is_slow:
            if run_start is None:
                run_start = pos
        else:
            if run_start is not None and pos - run_start >= _SLOW_CODON_RUN_MIN:
                hits.append(
                    PauseHit(kind="slow_codon_run", start=run_start * 3, end=pos * 3, score=float(pos - run_start),
                              notes=f"{pos - run_start} consecutive non-preferred codons")
                )
            run_start = None
    if run_start is not None and n_codons - run_start >= _SLOW_CODON_RUN_MIN:
        hits.append(
            PauseHit(kind="slow_codon_run", start=run_start * 3, end=n_codons * 3, score=float(n_codons - run_start),
                      notes=f"{n_codons - run_start} consecutive non-preferred codons")
        )
    return hits


def _scan_polyproline(cds_nt: str) -> List[PauseHit]:
    protein = translate(cds_nt)
    hits: List[PauseHit] = []
    i = 0
    n = len(protein)
    while i < n:
        if protein[i] != "P":
            i += 1
            continue
        j = i
        while j < n and protein[j] == "P":
            j += 1
        if j - i >= _POLYPROLINE_RUN_MIN:
            hits.append(
                PauseHit(kind="polyproline", start=i * 3, end=j * 3, score=float(j - i), notes=f"{j - i}x proline run")
            )
        i = j
    return hits


def scan_pauses(cds_nt: str, host, calc=None, temperature_c: Optional[float] = None) -> List[PauseHit]:
    """Positions and scores of candidate ribosomal pause sites inside one
    CDS (section 12): internal SD/anti-SD hybridization, slow-codon runs,
    and polyproline/stall motifs, all in ``cds_nt``-relative nt
    coordinates.

    ``calc`` is accepted for interface symmetry with the other scanners
    (spec section 25's ``scan_pauses(cds_nt, host, calc)``) but not used
    directly -- the internal-SD signal reuses
    ``rbsforge.thermo.duplex.best_hybridization`` directly rather than
    going through an ``RBSCalculator`` (there is no start codon to anchor
    a normal RBS prediction on here).
    """
    hits: List[PauseHit] = []
    hits.extend(_scan_internal_sd(cds_nt, host, temperature_c))
    hits.extend(_scan_slow_codon_runs(cds_nt, host))
    hits.extend(_scan_polyproline(cds_nt))
    hits.sort(key=lambda h: h.start)
    return hits
