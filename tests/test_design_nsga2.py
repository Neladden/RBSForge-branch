import random
import unittest

from operon.assembly import assemble
from operon.core import CDS, Operon, OperonHost
from operon.design import Design, design
from operon.design.nsga2 import _constrained_dominates, _dominates, _fast_nondominated_sort, _Individual, _mutate_operon
from operon.design.objectives import (
    protein_sequence_constraint,
    repeat_length_objective,
    synthesis_feasibility_constraint,
    tir_error_objective,
)
from operon.elongation.genetic_code import translate
from rbsforge.hostpack import get_hostpack
from rbsforge.rbs.calculator import RBSCalculator


def make_operon():
    # non-repetitive promoter/RBS/terminator flanks: the RBS in particular
    # must not itself already violate synthesis_feasibility_constraint
    # (e.g. a long A-homopolymer), or every individual in the initial
    # population starts infeasible and elitism has nothing to preserve.
    host = OperonHost(pack=get_hostpack("ecoli"))
    cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAUAA", target_tir=1000.0)
    return Operon(
        promoter_dna="CCACGTGCATGCATCGTAGC", cds_list=[cds1], rbs_list=["AGGAGGACAACTAAACGTGCATGCA"],
        terminator_dna="GCTAGCTACGATGC", host=host, intergenic_policy="free",
    )


class DominanceTest(unittest.TestCase):
    def test_strictly_better_in_all_objectives_dominates(self):
        self.assertTrue(_dominates([1.0, 1.0], [2.0, 2.0]))

    def test_equal_does_not_dominate(self):
        self.assertFalse(_dominates([1.0, 1.0], [1.0, 1.0]))

    def test_mixed_does_not_dominate(self):
        self.assertFalse(_dominates([1.0, 3.0], [2.0, 2.0]))

    def test_feasible_always_dominates_infeasible(self):
        feasible = _Individual(operon=None, objective_values=[100.0], feasible=True, n_violations=0)
        infeasible = _Individual(operon=None, objective_values=[0.0], feasible=False, n_violations=1)
        self.assertTrue(_constrained_dominates(feasible, infeasible))
        self.assertFalse(_constrained_dominates(infeasible, feasible))

    def test_fewer_violations_wins_among_infeasible(self):
        one_violation = _Individual(operon=None, objective_values=[0.0], feasible=False, n_violations=1)
        two_violations = _Individual(operon=None, objective_values=[0.0], feasible=False, n_violations=2)
        self.assertTrue(_constrained_dominates(one_violation, two_violations))


class NondominatedSortTest(unittest.TestCase):
    def test_front_zero_contains_only_nondominated(self):
        pop = [
            _Individual(operon=None, objective_values=[1.0, 4.0], feasible=True),
            _Individual(operon=None, objective_values=[4.0, 1.0], feasible=True),
            _Individual(operon=None, objective_values=[3.0, 3.0], feasible=True),  # dominated by neither of the above
            _Individual(operon=None, objective_values=[5.0, 5.0], feasible=True),  # dominated by all three above
        ]
        fronts = _fast_nondominated_sort(pop)
        self.assertEqual(set(fronts[0]), {0, 1, 2})
        self.assertEqual(fronts[-1], [3])


class MutateOperonTest(unittest.TestCase):
    def test_mutation_preserves_protein_when_it_touches_a_codon(self):
        operon = make_operon()
        original_protein = translate(operon.cds_list[0].aa_or_nt)
        rng = random.Random(1)
        for _ in range(30):
            operon = _mutate_operon(operon, rng)
            self.assertEqual(translate(operon.cds_list[0].aa_or_nt), original_protein)

    def test_mutation_respects_rbs_constraint_mask(self):
        host = OperonHost(pack=get_hostpack("ecoli"))
        cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCUAA", rbs_constraint="GGGGG" + "N" * 15)
        operon = Operon(
            promoter_dna="G" * 5, cds_list=[cds1], rbs_list=["GGGGG" + "AAAAAAAAAAAAAAA"],
            terminator_dna="", host=host, intergenic_policy="free",
        )
        rng = random.Random(2)
        for _ in range(30):
            operon = _mutate_operon(operon, rng)
            self.assertEqual(operon.rbs_list[0][:5], "GGGGG")


class DesignEndToEndTest(unittest.TestCase):
    def setUp(self):
        self.hostpack = get_hostpack("ecoli")
        self.calc = RBSCalculator(self.hostpack, temperature_c=37.0, use_vienna_if_available=False)
        self.operon = make_operon()

    def test_returns_a_nonempty_pareto_front(self):
        objectives = [tir_error_objective({"geneA": 1000.0}, self.calc), repeat_length_objective()]
        constraints = [synthesis_feasibility_constraint(), protein_sequence_constraint(self.operon)]
        rng = random.Random(42)
        results = design(self.operon, self.calc, objectives, constraints, population_size=16, n_generations=6, rng=rng)
        self.assertTrue(results)
        self.assertLessEqual(len(results), 20)
        for d in results:
            self.assertIsInstance(d, Design)
            self.assertEqual(len(d.objective_values), 2)

    def test_returned_front_is_actually_nondominated(self):
        objectives = [tir_error_objective({"geneA": 1000.0}, self.calc), repeat_length_objective()]
        constraints = [protein_sequence_constraint(self.operon)]
        rng = random.Random(7)
        results = design(self.operon, self.calc, objectives, constraints, population_size=16, n_generations=6, rng=rng)
        for i, a in enumerate(results):
            for j, b in enumerate(results):
                if i == j:
                    continue
                self.assertFalse(_dominates(a.objective_values, b.objective_values))

    def test_every_returned_design_still_satisfies_hard_constraints(self):
        objectives = [tir_error_objective({"geneA": 1000.0}, self.calc)]
        constraints = [protein_sequence_constraint(self.operon), synthesis_feasibility_constraint()]
        rng = random.Random(9)
        results = design(self.operon, self.calc, objectives, constraints, population_size=12, n_generations=5, rng=rng)
        for d in results:
            self.assertTrue(protein_sequence_constraint(self.operon)(d.operon, d.assembled))
            self.assertTrue(synthesis_feasibility_constraint()(d.operon, d.assembled))

    def test_results_are_deduplicated(self):
        objectives = [tir_error_objective({"geneA": 1000.0}, self.calc)]
        constraints = []
        rng = random.Random(11)
        results = design(self.operon, self.calc, objectives, constraints, population_size=10, n_generations=3, rng=rng)
        keys = [(tuple(c.aa_or_nt for c in d.operon.cds_list), tuple(d.operon.rbs_list)) for d in results]
        self.assertEqual(len(keys), len(set(keys)))

    def test_single_generation_still_returns_something(self):
        objectives = [tir_error_objective({"geneA": 1000.0}, self.calc)]
        rng = random.Random(13)
        results = design(self.operon, self.calc, objectives, [], population_size=6, n_generations=1, rng=rng)
        self.assertTrue(results)


if __name__ == "__main__":
    unittest.main()
