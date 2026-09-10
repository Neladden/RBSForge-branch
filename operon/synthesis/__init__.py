"""
Operon Calculator -- Module J: DNA synthesis complexity and restriction
sites (spec section 15). **Status: implemented** for the numeric hard-
rule table (15.1) and the RE-site hard constraint; hairpin/G-quadruplex/
i-motif detection is explicitly out of scope here -- see ``CLAUDE.md``.

Sources: Halper, Hossain & Salis, ACS Synth. Biol. 9:1436 (2020),
https://doi.org/10.1021/acssynbio.9b00460 (Synthesis Success Calculator).
The hard rules already define a feasible set on their own (section
15.1); the trained random-forest ``P(synthesis success)`` score they
describe as a later refinement is not implemented.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from operon.repeats import find_repeats, scan_motifs

HOMOPOLYMER_GC_MAX = 8   # no run of 9+ C or G
HOMOPOLYMER_AT_MAX = 12  # no run of 13+ A or T
TRINUCLEOTIDE_TANDEM_MAX_COPIES = 5   # no 6x period-3 repeat
DINUCLEOTIDE_TANDEM_MAX_COPIES = 9    # no 10x period-2 repeat
GC_WINDOW_20 = (20, 0.15, 0.90)
GC_WINDOW_100 = (100, 0.28, 0.76)
GC_TERMINAL_30 = (30, 0.24, 0.76)  # only checked when len(dna) > 60

# Local repeat density (section 15.1); coverage is derived from
# find_repeats, whose direct/inverted hits are seed(k)-length, not
# extended -- see CLAUDE.md. This makes density estimates conservative
# (a real, longer repeat is undercounted), not exact.
DENSITY_WINDOW_70_MAX = 0.90
DENSITY_WINDOW_500_MAX = 0.60
DENSITY_WHOLE_FRAGMENT_MAX = 0.69
DENSITY_SINGLE_REPEAT_MAX = 0.40


@dataclass
class SynthesisReport:
    homopolymer_hits: List[Tuple[str, int, int]] = field(default_factory=list)  # (base, start, end)
    tandem_hits: List[Tuple[int, int, int]] = field(default_factory=list)  # (period, start, end)
    gc_window_hits: List[Tuple[int, int, int, float]] = field(default_factory=list)  # (start, end, window, frac)
    density_hits: List[str] = field(default_factory=list)  # human-readable
    forbidden_re_hits: List[dict] = field(default_factory=list)

    @property
    def is_feasible(self) -> bool:
        return not (
            self.homopolymer_hits
            or self.tandem_hits
            or self.gc_window_hits
            or self.density_hits
            or self.forbidden_re_hits
        )


def _homopolymer_hits(seq: str) -> List[Tuple[str, int, int]]:
    hits = []
    n = len(seq)
    i = 0
    while i < n:
        j = i
        while j < n and seq[j] == seq[i]:
            j += 1
        run_len = j - i
        limit = HOMOPOLYMER_GC_MAX if seq[i] in "GC" else HOMOPOLYMER_AT_MAX
        if run_len > limit:
            hits.append((seq[i], i, j))
        i = j
    return hits


def _tandem_hits(seq: str, period: int, max_copies: int) -> List[Tuple[int, int, int]]:
    # A run of a single repeated character trivially satisfies any
    # period's recurrence (seq[j] == seq[j-period] for constant regions)
    # -- that's homopolymer territory (_homopolymer_hits), not a genuine
    # dinucleotide/trinucleotide *unit* repeat, so skip degenerate units.
    n = len(seq)
    hits = []
    i = 0
    while i + period < n:
        j = i + period
        while j < n and seq[j] == seq[j - period]:
            j += 1
        span = j - i
        copies = span / period
        unit = seq[i : i + period]
        if copies > max_copies and len(set(unit)) > 1:
            hits.append((period, i, j))
            i = j
        else:
            i += 1
    return hits


def _gc_window_hits(seq: str, window: int, lo: float, hi: float) -> List[Tuple[int, int, int, float]]:
    n = len(seq)
    if n < window:
        return []
    prefix = [0] * (n + 1)
    for i, ch in enumerate(seq):
        prefix[i + 1] = prefix[i] + (1 if ch in "GC" else 0)
    hits = []
    for i in range(0, n - window + 1):
        frac = (prefix[i + window] - prefix[i]) / window
        if frac < lo or frac > hi:
            hits.append((i, i + window, window, frac))
    return hits


def _repeat_density_hits(seq: str) -> List[str]:
    n = len(seq)
    report = find_repeats(seq, k=12)
    coverage = [False] * n
    max_single_repeat = 0
    for hit in report.all_hits():
        for start in hit.positions:
            end = min(n, start + hit.length)
            for p in range(start, end):
                coverage[p] = True
            max_single_repeat = max(max_single_repeat, end - start)

    hits = []
    if n and max_single_repeat / n > DENSITY_SINGLE_REPEAT_MAX:
        hits.append(f"a single repeat covers {max_single_repeat}/{n} nt (> {DENSITY_SINGLE_REPEAT_MAX:.0%} of the fragment)")

    def window_density_ok(window: int, max_frac: float) -> bool:
        if n < window:
            return True
        covered = sum(coverage[:window])
        if covered / window > max_frac:
            return False
        for i in range(window, n):
            covered += coverage[i] - coverage[i - window]
            if covered / window > max_frac:
                return False
        return True

    if not window_density_ok(70, DENSITY_WINDOW_70_MAX):
        hits.append(f"some 70-nt window exceeds {DENSITY_WINDOW_70_MAX:.0%} repeat coverage")
    if not window_density_ok(500, DENSITY_WINDOW_500_MAX):
        hits.append(f"some 500-nt window exceeds {DENSITY_WINDOW_500_MAX:.0%} repeat coverage")
    if n and sum(coverage) / n > DENSITY_WHOLE_FRAGMENT_MAX:
        hits.append(f"the whole fragment exceeds {DENSITY_WHOLE_FRAGMENT_MAX:.0%} repeat coverage")
    return hits


def synthesis_score(dna: str, forbidden_re_sites: Optional[List[str]] = None) -> SynthesisReport:
    """Homopolymer, tandem-repeat, windowed-GC%, and repeat-density
    features against the section 15.1 hard-rule table, plus a forbidden-
    restriction-site hard constraint (``n_forbidden_re > 0`` is invalid --
    section 15.2).

    These are meant to gate a Design-mode search (reject/repair a
    candidate that fails any of them), not just score it -- see
    ``report.is_feasible``.
    """
    seq = dna.upper()
    report = SynthesisReport()

    report.homopolymer_hits = _homopolymer_hits(seq)
    report.tandem_hits = (
        _tandem_hits(seq, 2, DINUCLEOTIDE_TANDEM_MAX_COPIES)
        + _tandem_hits(seq, 3, TRINUCLEOTIDE_TANDEM_MAX_COPIES)
    )

    window20, lo20, hi20 = GC_WINDOW_20
    window100, lo100, hi100 = GC_WINDOW_100
    report.gc_window_hits = _gc_window_hits(seq, window20, lo20, hi20) + _gc_window_hits(seq, window100, lo100, hi100)
    if len(seq) > 60:
        term_window, term_lo, term_hi = GC_TERMINAL_30
        head, tail = seq[:term_window], seq[-term_window:]
        for label, window_seq, offset in (("head", head, 0), ("tail", tail, len(seq) - term_window)):
            gc_frac = sum(1 for c in window_seq if c in "GC") / term_window
            if gc_frac < term_lo or gc_frac > term_hi:
                report.gc_window_hits.append((offset, offset + term_window, term_window, gc_frac))

    report.density_hits = _repeat_density_hits(seq)

    if forbidden_re_sites:
        report.forbidden_re_hits = scan_motifs(seq, forbidden_re_sites)

    return report
