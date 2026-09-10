"""
Operon Calculator -- Design mode: inverse RBS design (spec section 5,
Appendix A) and the multi-objective search over RBS nucleotides and
synonymous codons (spec section 18).

``design_rbs`` is the missing piece that blocks a faithful Operon Design
mode: RBSForge is Predict-only. It is a v1.0-style simulated annealer
(Salis, Mirsky & Voigt, Nat. Biotechnol. 27:946 (2009),
https://doi.org/10.1038/nbt.1568) plus IUPAC constraints (Salis, Methods
Enzymol. 498:19 (2011), section 4.4), calling ``rbsforge.predict`` as the
energy oracle. It is mono-cistronic -- never trust its returned TIR for a
CDS after the first; always re-score coupled CDSs against the assembled
mRNA (section 5.6, known failure mode 2).

The outer ``design`` search is NSGA-II (or any Pareto GA) over the
objective vector in section 18.2, not a weighted sum -- scalarizing hides
the TIR-vs-repeats-vs-HTISC tradeoffs the Operon Calculator's web tool
surfaces as a Pareto set.

Status: implemented.
"""

from .design_rbs import DesignResult, design_rbs
from .nsga2 import Design, design

__all__ = ["design_rbs", "DesignResult", "design", "Design"]
