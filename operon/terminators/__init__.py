"""
Operon Calculator -- Module H: internal transcriptional terminators
(spec section 13). **Status: implemented** for the presence layer of
both intrinsic and rho-dependent detection; the optional Chen-class
strength scoring (dG_U/dG_L) is not implemented -- see ``CLAUDE.md``.

There is no Salis Terminator Calculator; this is a scan plus a parts-
database append, not a fitted rate model of its own. Intrinsic
(rho-independent) detection follows Kingsford, Ayanbule & Salzberg,
Genome Biology 8:R22 (2007), https://doi.org/10.1186/gb-2007-8-2-r22.
Rho-dependent detection follows the RhoTermPredict approach (Di Salvo
et al. 2019, https://doi.org/10.1186/s12859-019-2704-x), simplified.
"""

from dataclasses import dataclass
from typing import List, Optional

from rbsforge.thermo.fold import BuiltinFolder

_COMPLEMENT = str.maketrans("ACGTU", "TGCAA")
_WC_PAIRS = {"AT", "TA", "AU", "UA", "GC", "CG"}
_WOBBLE_PAIRS = {"GT", "TG", "GU", "UG"}

# Intrinsic terminator geometry (spec section 13.2)
_STEM_RANGE = range(15, 3, -1)  # try longer stems first: 15 down to 4
_LOOP_RANGE = range(3, 9)  # 3-8 nt
_U_TRACT_WINDOW = 8
_U_TRACT_MIN_COUNT = 5  # crude gate: >=5 T/U in the 8 nt after the stem
_HAIRPIN_MFE_MAX = -7.0  # kcal/mol; crude gate from spec section 13.2

# Rho-dependent Rut-site geometry (spec section 13.3), simplified presence layer
_RUT_WINDOW = 78
_RUT_C_SPACING_RANGE = (11, 13)
_RUT_PAUSE_SEARCH = 100  # nt downstream to look for a pause (hairpin or pyrimidine tract)
_RUT_PYRIMIDINE_TRACT_MIN = 6


def _revcomp(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


def _pairable(a: str, b: str) -> bool:
    pair = a + b
    return pair in _WC_PAIRS or pair in _WOBBLE_PAIRS


@dataclass
class TermHit:
    kind: str  # "intrinsic" | "rho"
    start: int  # 0-based, forward-strand coordinates
    end: int  # exclusive
    strand: int  # +1 or -1
    score: Optional[float] = None  # kind-specific: intrinsic -> hairpin MFE; rho -> C count in the Rut window
    notes: str = ""


def _find_hairpin_at(seq: str, i: int):
    """Longest stem-loop starting exactly at ``i``, or ``None``."""
    n = len(seq)
    for stem_len in _STEM_RANGE:
        for loop_len in _LOOP_RANGE:
            j = i + stem_len + loop_len
            end = j + stem_len
            if end > n:
                continue
            arm5 = seq[i : i + stem_len]
            arm3 = seq[j:end]
            if all(_pairable(arm5[k], arm3[stem_len - 1 - k]) for k in range(stem_len)):
                return stem_len, loop_len, end
    return None


def _scan_intrinsic_one_strand(seq: str, folder: BuiltinFolder) -> List[TermHit]:
    n = len(seq)
    hits = []
    i = 0
    while i < n:
        found = _find_hairpin_at(seq, i)
        if found is None:
            i += 1
            continue
        stem_len, loop_len, hairpin_end = found
        u_tract = seq[hairpin_end : hairpin_end + _U_TRACT_WINDOW]
        u_count = u_tract.count("U") + u_tract.count("T")
        if u_count < _U_TRACT_MIN_COUNT:
            i += 1
            continue
        hairpin_seq = seq[i:hairpin_end]
        mfe = folder.fold(hairpin_seq).delta_g
        if mfe <= _HAIRPIN_MFE_MAX:
            hits.append(
                TermHit(
                    kind="intrinsic", start=i, end=hairpin_end + len(u_tract), strand=1,
                    score=mfe, notes=f"stem {stem_len} bp, loop {loop_len} nt, {u_count}/{len(u_tract)} U-tract, hairpin MFE {mfe:.2f}",
                )
            )
            i = hairpin_end
        else:
            i += 1
    return hits


def scan_intrinsic_terminators(
    dna: str, annotated_span: Optional[tuple] = None, annotated_window: int = 50
) -> List[TermHit]:
    """Rho-independent (hairpin + U-tract) terminator hits, both strands
    (section 13.2). A hit is "internal" (not the user-supplied 3'
    terminator) if its span doesn't fall within ``annotated_window`` nt of
    ``annotated_span`` on the same strand; pass ``annotated_span=None``
    (the default) to get every hit, unfiltered -- e.g. for scoring the
    intended terminator itself.
    """
    folder = BuiltinFolder()
    seq = dna.upper().replace("T", "U")
    n = len(seq)

    fwd_hits = _scan_intrinsic_one_strand(seq, folder)

    rev_seq = _revcomp(seq)
    rev_hits_local = _scan_intrinsic_one_strand(rev_seq, folder)
    rev_hits = [
        TermHit(kind="intrinsic", start=n - h.end, end=n - h.start, strand=-1, score=h.score, notes=h.notes)
        for h in rev_hits_local
    ]

    all_hits = fwd_hits + rev_hits
    if annotated_span is None:
        return all_hits

    a_start, a_end = annotated_span
    window_start, window_end = a_start - annotated_window, a_end + annotated_window
    return [h for h in all_hits if not (h.start >= window_start and h.end <= window_end)]


def _scan_rho_one_strand(seq: str, strand: int, n_total: int) -> List[TermHit]:
    n = len(seq)
    hits = []
    lo, hi = _RUT_C_SPACING_RANGE
    min_c, max_c = _RUT_WINDOW // hi, _RUT_WINDOW // lo
    for i in range(0, max(0, n - _RUT_WINDOW + 1)):
        window = seq[i : i + _RUT_WINDOW]
        c_count = window.count("C")
        g_count = window.count("G")
        if g_count == 0 or c_count / g_count <= 1:
            continue
        if not (min_c <= c_count <= max_c):
            continue
        downstream = seq[i + _RUT_WINDOW : i + _RUT_WINDOW + _RUT_PAUSE_SEARCH]
        has_pause = _has_pyrimidine_tract(downstream) or _find_hairpin_at(downstream, 0) is not None
        if not has_pause:
            continue
        start_fwd = i if strand == 1 else n_total - (i + _RUT_WINDOW)
        end_fwd = (i + _RUT_WINDOW) if strand == 1 else n_total - i
        hits.append(
            TermHit(
                kind="rho", start=start_fwd, end=end_fwd, strand=strand, score=float(c_count),
                notes=f"Rut candidate: {c_count} C / {g_count} G in {_RUT_WINDOW} nt, pause found downstream",
            )
        )
    return hits


def _has_pyrimidine_tract(seq: str) -> bool:
    run = 0
    for ch in seq:
        if ch in "CTU":
            run += 1
            if run >= _RUT_PYRIMIDINE_TRACT_MIN:
                return True
        else:
            run = 0
    return False


def scan_rho_terminators(dna: str) -> List[TermHit]:
    """Rho-dependent terminator hits via a simplified Rut-site + pause
    scan, both strands (section 13.3): a 78-nt window with C/G > 1 and a
    C count consistent with "a C roughly every 11-13 nt", followed within
    100 nt by a pause (a pyrimidine tract or an intrinsic-like hairpin).

    A simplified presence layer, not a trained RhoTermPredict/OPLS-DA
    classifier -- see ``CLAUDE.md``.
    """
    seq = dna.upper().replace("T", "U")
    n = len(seq)
    fwd_hits = _scan_rho_one_strand(seq, 1, n)
    rev_hits = _scan_rho_one_strand(_revcomp(seq), -1, n)
    return fwd_hits + rev_hits
