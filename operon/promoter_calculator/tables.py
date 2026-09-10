"""
Dinucleotide-level physical parameter tables used by the promoter energy
model's four numerical features (UP-element groove width/rigidity, ITR
DNA:RNA hybrid strength).

All four tables are complete over the 16 possible DNA dinucleotides (no
reverse-complement fallback is needed at lookup time).

Sources:

- ``GROOVE_ACCESS``: minor-groove width accessibility proxy per dinucleotide
  step, as used for the UP-element distal/proximal subregions in La Fleur,
  Hossain & Salis, Nat. Commun. 13:5159 (2022).
- ``PERSISTENCE``: sequence-dependent DNA bending rigidity (persistence
  length proxy), Geggier & Vologodskii, PNAS 107:15421 (2010).
- ``DNA_DNA_HYBRID``: DNA:DNA nearest-neighbor duplex free energies
  (SantaLucia-class parameters) used for the initial transcribed region
  (ITR) DNA-strand term.
- ``RNA_DNA_HYBRID``: RNA:DNA hybrid nearest-neighbor free energies
  (Sugimoto et al., Biochemistry 34:11211, 1995), used for the ITR
  R-loop / abortive-transcript term.
"""

GROOVE_ACCESS = {
    "CG": 43,
    "CA": 42, "TG": 42,
    "GG": 42, "CC": 42,
    "GC": 25,
    "GA": 22, "TC": 22,
    "TA": 14,
    "AG": 9, "CT": 9,
    "AA": 5, "TT": 5,
    "AC": 4, "GT": 4,
    "AT": 0,
}

PERSISTENCE = {
    "AA": 50.4, "TT": 50.4,
    "AC": 55.4, "GT": 55.4,
    "AG": 51.0, "CT": 51.0,
    "AT": 40.9,
    "CA": 46.7, "TG": 46.7,
    "CC": 41.7, "GG": 41.7,
    "CG": 56.0,
    "GA": 54.4, "TC": 54.4,
    "GC": 44.6,
    "TA": 44.7,
}

DNA_DNA_HYBRID = {
    "AA": -1.00, "TT": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45, "TG": -1.45,
    "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28,
    "GA": -1.30, "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.42, "CC": -1.42,
}

RNA_DNA_HYBRID = {
    "TT": -1.0, "TG": -2.1, "TC": -1.8, "TA": -0.9,
    "GT": -0.9, "GG": -2.1, "GC": -1.7, "GA": -0.9,
    "CT": -1.3, "CG": -2.7, "CC": -2.9, "CA": -1.1,
    "AT": -0.6, "AG": -1.5, "AC": -1.6, "AA": -0.2,
}
