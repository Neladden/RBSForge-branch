"""Small, dependency-free sequence utilities shared by the RBS scoring
engine.

Sequences are accepted as either DNA (T) or RNA (U) and are normalized to
RNA internally, since the thermodynamic parameter tables (nearest-neighbor
stacking energies, loop penalties) are RNA duplex/folding parameters -- the
mRNA is, after all, RNA, and the 16S rRNA anti-Shine-Dalgarno tail is RNA.
"""
from __future__ import annotations

_DNA_TO_RNA = str.maketrans("Tt", "Uu")
_COMPLEMENT_RNA = str.maketrans("ACGUacgu", "UGCAugca")

_VALID_RNA = frozenset("ACGU")


def to_rna(sequence: str) -> str:
    """Upper-case and normalize a DNA/RNA sequence to RNA (U instead of T)."""
    return sequence.strip().upper().translate(_DNA_TO_RNA)


def validate_rna(sequence: str) -> str:
    """Normalize to RNA and raise ValueError on any non-ACGU(T) character."""
    rna = to_rna(sequence)
    bad = set(rna) - _VALID_RNA
    if bad:
        raise ValueError(
            f"Sequence contains non-nucleotide characters: {sorted(bad)!r}"
        )
    return rna


def complement(sequence: str) -> str:
    """Return the RNA complement (not reversed) of an RNA/DNA sequence."""
    return to_rna(sequence).translate(_COMPLEMENT_RNA)


def reverse_complement(sequence: str) -> str:
    """Return the reverse-complement of an RNA/DNA sequence, as RNA."""
    return complement(sequence)[::-1]


def gc_content(sequence: str) -> float:
    rna = to_rna(sequence)
    if not rna:
        return 0.0
    gc = sum(1 for base in rna if base in "GC")
    return gc / len(rna)


def is_watson_crick(a: str, b: str) -> bool:
    """True if bases a and b form a canonical Watson-Crick pair (A-U, G-C)."""
    pair = a.upper() + b.upper()
    return pair in ("AU", "UA", "GC", "CG")


def is_wobble(a: str, b: str) -> bool:
    """True if bases a and b form a G-U wobble pair."""
    pair = a.upper() + b.upper()
    return pair in ("GU", "UG")


def sliding_windows(sequence: str, width: int):
    """Yield (start_index, window) for every window of `width` in sequence."""
    for i in range(0, len(sequence) - width + 1):
        yield i, sequence[i : i + width]
