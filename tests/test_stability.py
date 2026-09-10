import unittest

from operon.assembly import assemble
from operon.core import CDS, Operon, OperonHost
from operon.stability import (
    StabilityResult,
    has_5p_hairpin,
    mrna_stability,
    rnase_accessibility_score,
    ribosome_protection_dx,
)
from rbsforge.hostpack import get_hostpack


def make_operon(rbs, host):
    cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAUAA")
    operon = Operon(
        promoter_dna="G" * 5, cds_list=[cds1], rbs_list=[rbs], terminator_dna="",
        host=host, intergenic_policy="free",
    )
    return assemble(operon)


class RibosomeProtectionDxTest(unittest.TestCase):
    def setUp(self):
        self.host = OperonHost(pack=get_hostpack("ecoli"))

    def test_stronger_tir_gives_more_negative_dx(self):
        dx_strong = ribosome_protection_dx(50000.0, self.host)
        dx_weak = ribosome_protection_dx(10.0, self.host)
        self.assertLess(dx_strong, dx_weak)

    def test_very_weak_tir_can_give_positive_dx(self):
        dx = ribosome_protection_dx(0.5, self.host)
        self.assertGreater(dx, 0)

    def test_zero_tir_is_infinite_dx(self):
        self.assertEqual(ribosome_protection_dx(0.0, self.host), float("inf"))


class RnaseAccessibilityTest(unittest.TestCase):
    def test_unpaired_au_rich_region_scores_higher_than_gc_hairpin(self):
        au_rich = "A" * 30
        gc_hairpin = "GGGGG" + "AAAA" + "CCCCC" + "A" * 16
        self.assertGreater(rnase_accessibility_score(au_rich), rnase_accessibility_score(gc_hairpin))

    def test_score_is_nonnegative(self):
        self.assertGreaterEqual(rnase_accessibility_score("ACGUACGUACGU"), 0.0)


class Has5pHairpinTest(unittest.TestCase):
    def test_detects_a_real_hairpin_near_the_5p_end(self):
        utr = "GGGGG" + "AAAA" + "CCCCC" + "AAAAA"
        self.assertTrue(has_5p_hairpin(utr))

    def test_flat_utr_has_no_hairpin(self):
        self.assertFalse(has_5p_hairpin("A" * 20))

    def test_stem_shorter_than_min_is_not_reported(self):
        utr = "GG" + "AAAA" + "CC" + "AAAAA"  # 2 bp stem < default min_stem=4
        self.assertFalse(has_5p_hairpin(utr))


class MrnaStabilityTest(unittest.TestCase):
    def setUp(self):
        self.host = OperonHost(pack=get_hostpack("ecoli"))

    def test_returns_features_without_coefficients(self):
        a = make_operon("AAAAAAAAAAAAAAAAAAAA", self.host)
        result = mrna_stability(a, [50000.0], self.host)
        self.assertIsInstance(result, StabilityResult)
        self.assertIsNone(result.k_decay)
        self.assertIsNone(result.t_half)
        self.assertEqual(len(result.per_cistron_dx), 1)

    def test_returns_rate_only_when_coefficients_supplied(self):
        a = make_operon("AAAAAAAAAAAAAAAAAAAA", self.host)
        coeffs = dict(a0=0.0, a1=1.0, a2=1.0, a3=1.0, k0=0.01, b1=0.01, b2=-0.01)
        result = mrna_stability(a, [50000.0], self.host, coefficients=coeffs)
        self.assertIsNotNone(result.k_decay)
        self.assertIsNotNone(result.t_half)
        self.assertAlmostEqual(result.t_half, __import__("math").log(2) / result.k_decay)
        self.assertEqual(result.coefficients_used, coeffs)

        expected_log_mrna = (
            coeffs["a0"]
            + coeffs["a1"] * (-result.features.rnase_score)
            + coeffs["a2"] * (-result.features.dx)
            + coeffs["a3"] * (1.0 if result.features.has_5p_hairpin else 0.0)
        )
        self.assertAlmostEqual(result.log_mrna_level, expected_log_mrna)

    def test_multi_cistron_uses_minimum_protection(self):
        cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAUAA")
        cds2 = CDS(id="geneB", aa_or_nt="AUGCCCAAACCCAAACCCAAACCCAAACCCUAA")
        operon = Operon(
            promoter_dna="G" * 5, cds_list=[cds1, cds2], rbs_list=["AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAA"],
            terminator_dna="", host=self.host, intergenic_policy="free",
        )
        a = assemble(operon)
        result = mrna_stability(a, [50000.0, 5.0], self.host)
        self.assertEqual(result.features.dx, min(result.per_cistron_dx))

    def test_never_invents_coefficients_of_its_own(self):
        import inspect
        from operon import stability

        source = inspect.getsource(stability)
        # the module must not ship a default coefficients dict of its own
        self.assertNotIn("DEFAULT_COEFFICIENTS", source)


if __name__ == "__main__":
    unittest.main()
