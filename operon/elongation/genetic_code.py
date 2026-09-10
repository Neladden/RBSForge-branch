"""The standard genetic code (RNA codon -> one-letter amino acid, '*' for
stop). Universal across nearly all organisms (mitochondrial and a handful
of other exceptions aside, none of which this package targets) -- unlike
codon *usage bias*, this table is a fact of the code itself, not
organism-specific data to source or refit.
"""

GENETIC_CODE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = frozenset(c for c, aa in GENETIC_CODE.items() if aa == "*")

# amino acid -> every synonymous codon (built once, from the table above)
SYNONYMS = {}
for _codon, _aa in GENETIC_CODE.items():
    SYNONYMS.setdefault(_aa, []).append(_codon)
del _codon, _aa


def translate(cds_nt: str) -> str:
    """Translate an in-frame RNA/DNA CDS (from its start codon) to a
    one-letter amino acid string, stopping at (and excluding) the first
    stop codon. Raises ``ValueError`` if the length isn't a multiple of 3
    or an unrecognized codon (non-ACGU/T character) is hit."""
    seq = cds_nt.upper().replace("T", "U")
    if len(seq) % 3 != 0:
        raise ValueError(f"CDS length {len(seq)} is not a multiple of 3")
    protein = []
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        aa = GENETIC_CODE.get(codon)
        if aa is None:
            raise ValueError(f"unrecognized codon {codon!r} at nt offset {i}")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)
