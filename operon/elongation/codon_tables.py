"""
*E. coli* codon-usage weight tables (spec section 8.2).

**These are simplified, illustrative weights, not a measured codon-usage
table.** They encode, per amino acid, the single codon most commonly
cited as *E. coli*'s preferred choice in highly expressed genes (the same
consensus picks widely used across codon-optimization tools) at weight
``1.0``, with every other synonym at a flat, lower placeholder weight
(``0.3``) -- i.e. "prefer this one codon" rather than a real relative-
frequency ranking among the rest. Real codon bias is driven by tRNA
abundance (Ikemura's rule) and is not simply "prefer G/C at the third
position"; getting the *relative order of the non-preferred synonyms*
right requires an actual measured usage table, which this is not.

Before using this for a real design, replace it with real data --
e.g. the Codon Usage Database (https://www.kazusa.or.jp/codon/) or a
codon usage table computed directly from your expression host's own
highly-expressed genes (ribosomal proteins are the standard choice) and
whole genome. ``operon.core.OperonHost.codon_weights_high`` /
``codon_weights_balanced`` exist precisely so a caller can supply real
data instead of this placeholder -- see spec section 8.2's `table_high`
(highly-expressed-gene usage) vs `table_genomic` (whole-genome usage).
"""

from .genetic_code import SYNONYMS

# amino acid -> single most commonly cited *E. coli* K-12 preferred codon.
_PREFERRED_CODON = {
    "A": "GCG", "R": "CGC", "N": "AAC", "D": "GAU", "C": "UGC",
    "Q": "CAG", "E": "GAA", "G": "GGC", "H": "CAU", "I": "AUC",
    "L": "CUG", "K": "AAA", "M": "AUG", "F": "UUU", "P": "CCG",
    "S": "AGC", "T": "ACC", "W": "UGG", "Y": "UAU", "V": "GUG",
}

NON_PREFERRED_WEIGHT = 0.3


def _build_weight_table():
    table = {}
    for aa, codons in SYNONYMS.items():
        if aa == "*":
            continue
        preferred = _PREFERRED_CODON.get(aa)
        table[aa] = {c: (1.0 if c == preferred else NON_PREFERRED_WEIGHT) for c in codons}
    return table


# "Highly Translated": always the max-w codon (spec section 8.2) -- with
# only two weight tiers here, that's just _PREFERRED_CODON restated as
# per-codon weights.
ECOLI_HIGHLY_TRANSLATED = _build_weight_table()

# "Balanced" mixes toward a target TIR between this and the genomic
# table (section 8.2's `choose_codon`). Without a real whole-genome usage
# table, use the same weights as a documented, honestly-flat prior --
# the mixing math still works, it just has nothing genuinely different
# to mix toward yet.
ECOLI_GENOMIC = _build_weight_table()
