import unittest

from rbsforge.rbs.calculator import RBSCalculator
from operon.htisc import DEFAULT_TAU_HTISC, scan_htisc


class ScanHtiscTest(unittest.TestCase):
    def setUp(self):
        self.calc = RBSCalculator.for_species("ecoli", temperature_c=37.0, use_vienna_if_available=False)
        annotated_utr = "AAAAAAAAAAAGGAGGAAAAA"
        annotated_cds = "AUG" + "AAACAAACAAA"
        internal_site = "AGGAGGAAAAA" + "AUG" + "CCCCCCCCCC"
        self.mrna = annotated_utr + annotated_cds + internal_site + "CAAACAAACAAA"
        self.annotated_start = self.mrna.index("AUG")

    def test_excludes_annotated_start(self):
        hits = scan_htisc(self.mrna, [self.annotated_start], self.calc, tau_htisc=100.0)
        self.assertTrue(all(h.index != self.annotated_start for h in hits))

    def test_finds_strong_internal_start(self):
        hits = scan_htisc(self.mrna, [self.annotated_start], self.calc, tau_htisc=100.0)
        internal_index = self.mrna.index("AUG", self.annotated_start + 3)
        self.assertTrue(any(h.index == internal_index for h in hits))

    def test_threshold_excludes_weak_hits(self):
        hits = scan_htisc(self.mrna, [self.annotated_start], self.calc, tau_htisc=10_000_000.0)
        self.assertEqual(hits, [])

    def test_leaderless_starts_are_never_hits(self):
        leaderless_mrna = "C" * 30 + "AUG" + "AAACCCCCCCCCCCCCCCCCCC"
        hits = scan_htisc(leaderless_mrna, [], self.calc, tau_htisc=0.0)
        self.assertEqual(hits, [])

    def test_default_tau_is_1000_au_on_v1_scale(self):
        self.assertEqual(DEFAULT_TAU_HTISC, 1000.0)

    def test_hits_ranked_strongest_first(self):
        hits = scan_htisc(self.mrna, [self.annotated_start], self.calc, tau_htisc=0.0)
        rates = [h.v1_style_rate for h in hits]
        self.assertEqual(rates, sorted(rates, reverse=True))


if __name__ == "__main__":
    unittest.main()
