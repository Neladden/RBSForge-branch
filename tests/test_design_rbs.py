import random
import unittest

from operon.design.design_rbs import (
    DesignResult,
    TOL_KCAL_MOL,
    _has_start_codon,
    _matches_code,
    design_rbs,
    dg_target_from_tir,
)
from rbsforge.hostpack import get_hostpack

CDS = "AUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAUAA" + "CAA" * 10


class DgTargetTest(unittest.TestCase):
    def test_stronger_target_tir_gives_more_negative_dg_target(self):
        strong = dg_target_from_tir(50000.0, beta=0.45)
        weak = dg_target_from_tir(10.0, beta=0.45)
        self.assertLess(strong, weak)

    def test_matches_v1_style_rate_inverse(self):
        from rbsforge.scoring import v1_style_rate

        # V1_LOG_K=7.824 is the published, rounded constant (not full-
        # precision ln(2500)), so the round-trip is close but not exact.
        target = 2500.0
        dg = dg_target_from_tir(target, beta=0.45)
        self.assertAlmostEqual(v1_style_rate(dg, rt_eff=1.0 / 0.45), target, delta=1.0)


class IupacTest(unittest.TestCase):
    def test_n_matches_every_base(self):
        for base in "ACGU":
            self.assertTrue(_matches_code(base, "N"))

    def test_r_matches_only_purines(self):
        self.assertTrue(_matches_code("A", "R"))
        self.assertTrue(_matches_code("G", "R"))
        self.assertFalse(_matches_code("C", "R"))
        self.assertFalse(_matches_code("U", "R"))

    def test_literal_code_matches_only_itself(self):
        self.assertTrue(_matches_code("G", "G"))
        self.assertFalse(_matches_code("A", "G"))


class HasStartCodonTest(unittest.TestCase):
    def test_detects_aug(self):
        self.assertTrue(_has_start_codon("CCCCAUGCCCC"))

    def test_no_start_codon_in_clean_sequence(self):
        self.assertFalse(_has_start_codon("CCCCCCCCCCCC"))


class DesignRbsTest(unittest.TestCase):
    def setUp(self):
        self.hostpack = get_hostpack("ecoli")

    def test_hits_a_moderate_target_within_tolerance(self):
        rng = random.Random(42)
        result = design_rbs(CDS, target_tir=1000.0, hostpack=self.hostpack, max_iter=4000, rng=rng)
        self.assertIsInstance(result, DesignResult)
        self.assertTrue(result.hit_tol)
        self.assertLessEqual(result.objective, TOL_KCAL_MOL)

    def test_returned_rbs_has_no_internal_start_codon(self):
        rng = random.Random(1)
        result = design_rbs(CDS, target_tir=500.0, hostpack=self.hostpack, max_iter=3000, rng=rng)
        self.assertFalse(_has_start_codon(result.rbs))

    def test_returned_rbs_within_length_bounds(self):
        rng = random.Random(2)
        result = design_rbs(CDS, target_tir=500.0, hostpack=self.hostpack, rbs_len_range=(20, 35), max_iter=2000, rng=rng)
        self.assertGreaterEqual(len(result.rbs), 20)
        self.assertLessEqual(len(result.rbs), 35)

    def test_hit_tol_designs_span_a_wide_tir_range(self):
        for target in (10.0, 1000.0, 50000.0):
            rng = random.Random(7)
            result = design_rbs(CDS, target_tir=target, hostpack=self.hostpack, max_iter=4000, rng=rng)
            self.assertTrue(result.hit_tol, f"failed to converge for target_tir={target}")
            # within roughly an order of magnitude given the 0.25 kcal/mol tolerance
            self.assertGreater(result.predicted_tir, target / 3)
            self.assertLess(result.predicted_tir, target * 3)

    def test_fixed_length_iupac_mask_is_respected(self):
        mask = "N" * 20
        rng = random.Random(3)
        result = design_rbs(CDS, target_tir=5000.0, hostpack=self.hostpack, constraint_iupac=mask, max_iter=3000, rng=rng)
        self.assertEqual(len(result.rbs), len(mask))

    def test_literal_positions_in_mask_are_frozen(self):
        mask = "GGGGG" + "N" * 15
        rng = random.Random(4)
        result = design_rbs(CDS, target_tir=3000.0, hostpack=self.hostpack, constraint_iupac=mask, max_iter=3000, rng=rng)
        self.assertEqual(result.rbs[:5], "GGGGG")

    def test_constant_upstream_is_not_mutated_and_precedes_rbs(self):
        upstream = "AUGAAACCCGGGUAAUAAUAA"
        rng = random.Random(5)
        result = design_rbs(CDS, target_tir=2000.0, hostpack=self.hostpack, constant_upstream=upstream, max_iter=2000, rng=rng)
        # re-scoring with the same upstream + returned RBS should reproduce dg_total exactly
        from rbsforge.rbs.calculator import RBSCalculator

        calc = RBSCalculator(self.hostpack)
        mrna = upstream + result.rbs + CDS
        recheck = calc.predict_one(mrna, len(upstream) + len(result.rbs))
        self.assertAlmostEqual(recheck.delta_g_total, result.dg_total, places=9)

    def test_deterministic_with_seeded_rng(self):
        r1 = design_rbs(CDS, target_tir=1000.0, hostpack=self.hostpack, max_iter=1500, rng=random.Random(99))
        r2 = design_rbs(CDS, target_tir=1000.0, hostpack=self.hostpack, max_iter=1500, rng=random.Random(99))
        self.assertEqual(r1.rbs, r2.rbs)

    def test_max_iter_cap_still_returns_best_so_far(self):
        rng = random.Random(6)
        result = design_rbs(CDS, target_tir=1000.0, hostpack=self.hostpack, max_iter=1, rng=rng)
        self.assertIsInstance(result, DesignResult)
        self.assertGreaterEqual(result.n_evals, 1)


if __name__ == "__main__":
    unittest.main()
