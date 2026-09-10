"""
Operon Calculator -- Module C: translation elongation rate (TER) and
synonymous codon recoding (spec section 8).

Independent of coupling; operates on one CDS at a time. Used as a score
in Predict mode and as a mutator in Design mode.

Status: planned, not yet implemented.
"""


def recode_cds(aa: str, table: dict, rng) -> str:
    """Recode a protein sequence to DNA using a synonymous codon table
    (Highly Translated or Balanced, section 8.2). A neighborhood operator
    for the Design-mode GA (section 8.3), not a one-shot rewrite: a single
    greedy CAI pass fights every other objective (repeats, promoters,
    HTISC, RNase sites).
    """
    raise NotImplementedError("Module C (elongation/recoding) is a prepared subsection, not yet implemented")


def ter_nt_s(cds_nt: str, host) -> float:
    """Average translation elongation rate in nt/s, as a harmonic mean
    over per-codon elongation rates (section 8.1) -- one slow codon
    dominates dwell time, so use the harmonic mean, not the arithmetic
    one::

        TER_cds = 3 / mean(1 / r_elong_codon[c] for c in cds)
    """
    raise NotImplementedError("Module C (elongation/recoding) is a prepared subsection, not yet implemented")
