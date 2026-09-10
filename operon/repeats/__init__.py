"""
Operon Calculator -- Module I: repeats and genetic stability (spec
section 14), plus transposon insertion / phage att site scanning.

Sources: Hossain et al., Nat. Biotechnol. 38:1466 (2020),
https://doi.org/10.1038/s41587-020-0584-2 (Nonrepetitive Parts
Calculator); the Operon Calculator's published objective is "fewer
repetitive DNA sequences above 12 bp".

Status: planned, not yet implemented.
"""


def find_repeats(dna: str, k: int = 12):
    """Direct, reverse-complement, tandem, and terminal repeats of length
    >= ``k``, both strands (section 14.1). Seed with k-mers (canonical
    strand + RC as the key), extend, merge overlapping matches -- O(L)
    expected with a dict of 4^k keys.
    """
    raise NotImplementedError("Module I (repeats) is a prepared subsection, not yet implemented")


def scan_motifs(dna: str, motifs: list) -> list:
    """Exact (or 1-mismatch) both-strand hits for a host-specific IS-
    element / phage att-site motif table (section 14.3)."""
    raise NotImplementedError("Module I (repeats) is a prepared subsection, not yet implemented")
