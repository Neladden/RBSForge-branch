"""
Operon Calculator -- Module A: operon assembly (spec section 3).

The only place sequence is concatenated. Independent of every biophysical
model. Layout::

    [promoter][5' UTR / RBS_1][CDS_1][intergenic_1][RBS_2][CDS_2] ... [terminator]

Status: planned, not yet implemented.
"""

from operon.core import Assembled, Operon


def assemble(operon: Operon) -> Assembled:
    """Concatenate promoter + [RBS_i + CDS_i] + terminator into one
    assembled DNA/mRNA molecule, computing each junction's intergenic
    distance ``d`` (section 3.1) along the way.

    ``d = start_{i+1} - stop_i_end``; negative means overlap.
    """
    raise NotImplementedError("Module A (assembly) is a prepared subsection, not yet implemented")
