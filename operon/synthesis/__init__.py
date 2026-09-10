"""
Operon Calculator -- Module J: DNA synthesis complexity and restriction
sites (spec section 15).

Sources: Halper, Hossain & Salis, ACS Synth. Biol. 9:1436 (2020),
https://doi.org/10.1021/acssynbio.9b00460 (Synthesis Success Calculator).
The hard rules (homopolymer runs, windowed GC%, tandem repeats) already
define a feasible set and should be used as GA constraints, not just
scores; Module I (repeats) does most of this work already.

Status: planned, not yet implemented.
"""


def synthesis_score(dna: str):
    """Homopolymer, windowed-GC%, and repeat-density features against the
    section 15.1 hard-rule table. A user-supplied restriction-site list is
    a hard constraint (``n_forbidden_re > 0`` is invalid), not part of
    this score.
    """
    raise NotImplementedError("Module J (synthesis) is a prepared subsection, not yet implemented")
