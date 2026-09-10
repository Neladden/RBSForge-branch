"""Numerical feature calculators for the sigma70 promoter energy model."""

from .tables import DNA_DNA_HYBRID, GROOVE_ACCESS, PERSISTENCE, RNA_DNA_HYBRID


def groove_width(seq: str) -> float:
    """Sum of the minor-groove accessibility proxy over non-overlapping
    dinucleotide steps of ``seq`` (used on each 12-nt half of the UP
    element)."""
    total = 0
    for i in range(0, len(seq) - 1, 2):
        total += GROOVE_ACCESS[seq[i:i + 2]]
    return float(total)


def rigidity(seq: str) -> float:
    """Mean sequence-dependent bending rigidity over non-overlapping
    dinucleotide steps of ``seq``, evaluated across the UP element, the
    -35 hexamer, and the first 14 nt of the spacer."""
    total = 0.0
    for i in range(0, len(seq), 2):
        total += PERSISTENCE[seq[i:i + 2]]
    return total / len(seq)


def dna_rna_hybrid_energy(itr: str) -> tuple[float, float, float]:
    """DNA:DNA and RNA:DNA hybrid free energies over the first 15 nt of the
    initial transcribed region (the length of a long abortive transcript),
    stepped as non-overlapping dinucleotides.

    Returns ``(dg_dna, dg_rna, dg_hybrid)`` where ``dg_hybrid = dg_dna -
    dg_rna`` is the R-loop free-energy differential fed into the model as
    ``dG_ITR``'s numerical feature.
    """
    dg_dna = 0.0
    dg_rna = 0.0
    window = itr[0:15]
    for i in range(0, len(window) - 1, 2):
        pair = window[i:i + 2]
        dg_dna += DNA_DNA_HYBRID[pair]
        dg_rna += RNA_DNA_HYBRID[pair]
    return dg_dna, dg_rna, dg_dna - dg_rna
