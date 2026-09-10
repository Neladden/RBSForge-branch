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

Status: planned, not yet implemented.
"""


def design_rbs(cds_nt: str, target_tir: float, hostpack, constraint_iupac=None,
                constant_upstream=None, rbs_init=None, rbs_len_range=(20, 35),
                max_iter: int = 10000, temperature_c=None):
    """Simulated-annealing inverse RBS design against ``rbsforge.predict``
    (Appendix A). Objective::

        dG_target = 2.222 * (7.824 - ln(target_tir))
        O(seq)    = |predict(constant_upstream + rbs + cds_nt).dG_total - dG_target|

    Metropolis + the v1.0 RT-annealing schedule (RT_init=0.6,
    every 50 moves: RT /= 2 if accept_ratio > 0.20, RT *= 2 if < 0.01),
    tol=0.25 kcal/mol, max_iter default 10000. Returns best-so-far even if
    tol is missed, with ``hit_tol: bool``.
    """
    raise NotImplementedError("design_rbs is a prepared subsection, not yet implemented")


def design(operon, calc, objectives, constraints) -> list:
    """NSGA-II Pareto search over RBS_i nucleotides (within their IUPAC
    masks) and synonymous codon choices, re-assembling and re-scoring
    (coupling included) every candidate (section 18.3). Returns a
    non-dominated set of designs (5-20, de-duplicated by 12-mer identity),
    not a single winner.
    """
    raise NotImplementedError("Design mode is a prepared subsection, not yet implemented")
