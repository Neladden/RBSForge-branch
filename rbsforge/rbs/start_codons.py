"""Start codon detection and delta_G_start lookup."""
from __future__ import annotations

from typing import Dict, Iterable, List, Tuple

from ..seqtools import validate_rna

# v1.0 scans AUG/GUG/UUG only (CUG is tabulated but excluded from the
# scan -- model summary section 1 and 4.4). Callers that want CUG
# candidates too can pass a wider `allowed_codons` explicitly.
DEFAULT_SCANNED_CODONS = ("AUG", "GUG", "UUG")


def find_start_codons(
    sequence: str,
    allowed_codons: Iterable[str] = DEFAULT_SCANNED_CODONS,
) -> List[Tuple[int, str]]:
    """Return every (0-based index, codon) occurrence of an allowed start
    codon in `sequence`, scanning all three frames (a start codon near an
    RBS need not be in one fixed reading frame relative to the transcript
    start)."""
    seq = validate_rna(sequence)
    allowed = {validate_rna(c) for c in allowed_codons}
    hits: List[Tuple[int, str]] = []
    for i in range(len(seq) - 2):
        codon = seq[i : i + 3]
        if codon in allowed:
            hits.append((i, codon))
    return hits


def start_codon_energy(codon: str, energies: Dict[str, float]) -> float:
    """delta_G_start for a given codon, using the HostPack's table with a
    conservative fallback for unrecognized/non-canonical codons."""
    codon = validate_rna(codon)
    if codon in energies:
        return energies[codon]
    return max(energies.values(), default=0.0) + 1.0
