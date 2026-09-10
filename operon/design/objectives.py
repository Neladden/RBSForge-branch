"""
Ready-made objective/constraint functions for ``operon.design.design``,
built from the now-implemented modules (spec section 18.2's `f1..f12`).
Not exhaustive -- a real Design run will usually mix these with problem-
specific ones -- but enough to run `design()` end to end against real
scanners rather than only synthetic test objectives.

Every objective returns a float to **minimize**; every constraint
returns ``True`` for "satisfied."
"""

import math
from typing import Dict, List, Optional

from operon.core import Assembled, Operon
from operon.elongation.genetic_code import translate
from operon.htisc import DEFAULT_TAU_HTISC, scan_htisc
from operon.repeats import find_repeats
from operon.synthesis import synthesis_score


def tir_error_objective(target_tirs: Dict[str, float], calc):
    """``f1`` (section 18.2): sum of |log(actual v1_style_rate) -
    log(target)| over every CDS with a target in ``target_tirs`` (keyed
    by ``CDS.id``), scored mono-cistronically via ``calc.predict_one`` on
    the assembled mRNA. A leaderless start (or a CDS with no scorable
    result) contributes a large fixed penalty rather than crashing.

    Mono-cistronic, not coupling-aware -- for a multi-cistron operon
    where junction coupling matters, build an equivalent objective around
    ``operon.coupling.coupled_tirs`` instead (that needs a Vienna-backed
    ``calc.folder`` and is more expensive per evaluation, which is why
    this simpler default doesn't call it automatically).
    """
    leaderless_penalty = 20.0  # kcal/mol-equivalent log-TIR-error units; large relative to a realistic target

    def objective(operon: Operon, assembled: Assembled) -> float:
        error = 0.0
        for cds, start in zip(operon.cds_list, assembled.starts):
            target = target_tirs.get(cds.id)
            if target is None:
                continue
            result = calc.predict_one(assembled.mrna, start)
            if result.leaderless or not result.v1_style_rate or result.v1_style_rate <= 0:
                error += leaderless_penalty
                continue
            error += abs(math.log(result.v1_style_rate) - math.log(target))
        return error

    return objective


def repeat_length_objective(k: int = 12):
    """``f7``: ``max(0, longest_repeat - k)`` over the assembled DNA."""

    def objective(operon: Operon, assembled: Assembled) -> float:
        report = find_repeats(assembled.dna, k=k)
        return float(max(0, report.longest() - k))

    return objective


def htisc_count_objective(calc, tau_htisc: float = DEFAULT_TAU_HTISC):
    """``f6``: count of highly translated internal start codons."""

    def objective(operon: Operon, assembled: Assembled) -> float:
        hits = scan_htisc(assembled.mrna, assembled.starts, calc, tau_htisc=tau_htisc)
        return float(len(hits))

    return objective


def synthesis_feasibility_constraint(forbidden_re_sites: Optional[List[str]] = None):
    """Hard constraint: no homopolymer/tandem/GC-window violation and no
    forbidden restriction site (section 15, section 18.2's hard-
    constraint list)."""

    def constraint(operon: Operon, assembled: Assembled) -> bool:
        return synthesis_score(assembled.dna, forbidden_re_sites=forbidden_re_sites).is_feasible

    return constraint


def protein_sequence_constraint(reference_operon: Operon):
    """Hard constraint: every CDS still translates to the same protein as
    in ``reference_operon`` (section 18.2's "protein sequence unchanged").
    Compares by position (``operon.cds_list[i]`` against
    ``reference_operon.cds_list[i]``), not by ``CDS.id``, so it still
    works if a mutation operator ever renames a CDS.
    """
    reference_proteins = [translate(cds.aa_or_nt) for cds in reference_operon.cds_list]

    def constraint(operon: Operon, assembled: Assembled) -> bool:
        if len(operon.cds_list) != len(reference_proteins):
            return False
        for cds, ref_protein in zip(operon.cds_list, reference_proteins):
            try:
                if translate(cds.aa_or_nt) != ref_protein:
                    return False
            except ValueError:
                return False
        return True

    return constraint
