"""
Operon Calculator -- Module I: repeats and genetic stability (spec
section 14), plus transposon insertion / phage att site scanning.
**Status: implemented** (``find_repeats`` reports minimum-length,
k-mer-level hits -- see ``CLAUDE.md`` for what "not maximally extended"
means and why that's a deliberate scope boundary, not a gap).

Sources: Hossain et al., Nat. Biotechnol. 38:1466 (2020),
https://doi.org/10.1038/s41587-020-0584-2 (Nonrepetitive Parts
Calculator); the Operon Calculator's published objective is "fewer
repetitive DNA sequences above 12 bp".
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

_COMPLEMENT = str.maketrans("ACGTU", "TGCAA")

# Tandem-repeat periods to scan (spec section 14.1's "tandem, period >= 1").
# Capped at 6 (homopolymer through hexanucleotide) -- this is the
# synthesis/stability-relevant range (Module J's own hard rules only go
# up to trinucleotide/dinucleotide tandems); longer periods are covered
# by the direct-repeat scan instead once span >= k.
_TANDEM_PERIODS = range(1, 7)


def _revcomp(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


@dataclass
class RepeatHit:
    kind: str  # "direct" | "inverted" | "tandem" | "terminal"
    unit: str  # the repeated k-mer (direct/inverted) or tandem period unit
    length: int  # len(unit) for direct/inverted; total span for tandem/terminal
    positions: List[int] = field(default_factory=list)  # forward-strand start positions
    rc_positions: Optional[List[int]] = None  # inverted only: positions of the reverse complement
    period: Optional[int] = None  # tandem/terminal only
    copies: Optional[int] = None  # tandem/terminal only


@dataclass
class RepeatReport:
    direct: List[RepeatHit] = field(default_factory=list)
    inverted: List[RepeatHit] = field(default_factory=list)
    tandem: List[RepeatHit] = field(default_factory=list)
    terminal: List[RepeatHit] = field(default_factory=list)

    def longest(self) -> int:
        lengths = [h.length for group in (self.direct, self.inverted, self.tandem, self.terminal) for h in group]
        return max(lengths, default=0)

    def all_hits(self) -> List[RepeatHit]:
        return [*self.direct, *self.inverted, *self.tandem, *self.terminal]


def _tandem_runs(seq: str, min_span: int):
    """Yield (start, period, copies) for every maximal tandem run (period
    in ``_TANDEM_PERIODS``) whose total span (period * copies) is >=
    ``min_span``. Each run is reported once, at its leftmost start,
    keyed by the smallest period that explains it."""
    n = len(seq)
    reported = [False] * n  # position already covered by a shorter-period run
    for period in _TANDEM_PERIODS:
        i = 0
        while i + period < n:
            if reported[i]:
                i += 1
                continue
            j = i + period
            while j < n and seq[j] == seq[j - period]:
                j += 1
            span = j - i
            copies = span / period
            if span >= min_span and copies >= 2:
                yield i, period, span
                for p in range(i, j):
                    reported[p] = True
                i = j
            else:
                i += 1


def find_repeats(dna: str, k: int = 12) -> RepeatReport:
    """Direct, reverse-complement (inverted), tandem, and terminal
    repeats of length >= ``k``, both strands (section 14.1).

    Seeds with ``k``-mers (an O(L) dict build); direct and inverted hits
    are reported at the seed (``k``-mer) length, one hit per repeated
    unit bundling every position it occurs at -- not extended to the
    maximal exact match a pair of occurrences might share (see
    ``CLAUDE.md``). Tandem/terminal repeats are found independently by a
    short-period run scan (periods 1-6), not derived from the k-mer index.
    """
    seq = dna.upper()
    n = len(seq)
    report = RepeatReport()
    if n < k:
        return report

    kmer_positions: Dict[str, List[int]] = {}
    for i in range(n - k + 1):
        kmer_positions.setdefault(seq[i : i + k], []).append(i)

    seen_inverted_units = set()
    for kmer, positions in kmer_positions.items():
        if len(positions) >= 2:
            report.direct.append(RepeatHit(kind="direct", unit=kmer, length=k, positions=list(positions)))

        rc = _revcomp(kmer)
        if rc == kmer or kmer in seen_inverted_units:
            continue
        rc_positions = kmer_positions.get(rc)
        if rc_positions:
            seen_inverted_units.add(kmer)
            seen_inverted_units.add(rc)
            report.inverted.append(
                RepeatHit(kind="inverted", unit=kmer, length=k, positions=list(positions), rc_positions=list(rc_positions))
            )

    for start, period, span in _tandem_runs(seq, min_span=k):
        unit = seq[start : start + period]
        copies = span // period + (1 if span % period else 0)
        hit = RepeatHit(kind="tandem", unit=unit, length=span, positions=[start], period=period, copies=copies)
        report.tandem.append(hit)

    # Terminal: tandem of length >= 5 at either end of the fragment (a
    # lower bar than the general k-bp tandem threshold above).
    for start, period, span in _tandem_runs(seq, min_span=5):
        if start == 0 or start + span == n:
            unit = seq[start : start + period]
            copies = span // period + (1 if span % period else 0)
            report.terminal.append(
                RepeatHit(kind="terminal", unit=unit, length=span, positions=[start], period=period, copies=copies)
            )

    return report


def scan_motifs(dna: str, motifs: List[str], max_mismatches: int = 0) -> list:
    """Exact (or up-to-``max_mismatches``-mismatch) both-strand hits for a
    host-specific IS-element / phage att-site motif table (section 14.3).

    Returns a list of ``{"motif": str, "position": int, "strand": int,
    "mismatches": int}`` dicts. Brute-force (fine for a handful of short
    motifs against operon-scale DNA); not meant for genome-scale input.
    """
    seq = dna.upper()
    n = len(seq)
    hits = []
    for motif in motifs:
        m = motif.upper()
        L = len(m)
        if L == 0 or L > n:
            continue
        rc_m = _revcomp(m)
        for strand, needle in ((1, m), (-1, rc_m)):
            for i in range(n - L + 1):
                window = seq[i : i + L]
                mismatches = sum(1 for a, b in zip(window, needle) if a != b)
                if mismatches <= max_mismatches:
                    hits.append({"motif": motif, "position": i, "strand": strand, "mismatches": mismatches})
    return hits
