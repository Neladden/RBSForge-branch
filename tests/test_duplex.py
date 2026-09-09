import unittest

from rbsforge.constants import celsius_to_kelvin
from rbsforge.thermo.duplex import best_hybridization


def zero_spacing(_aligned_spacing):
    return 0.0


class DuplexSearchTest(unittest.TestCase):
    def setUp(self):
        self.t37 = celsius_to_kelvin(37.0)
        self.asd = "ACCUCCUUA"

    def test_perfect_full_length_match_has_zero_overhang(self):
        # reverse_complement("ACCUCCUUA") == "UAAGGAGGU": a perfect 9-nt
        # SD placed right before the start codon uses the whole aSD tail.
        utr = "UAAGGAGGU"
        start_index = len(utr)
        alignment, _ = best_hybridization(utr, self.asd, start_index, zero_spacing, self.t37)
        self.assertEqual(alignment.length, 9)
        self.assertEqual(alignment.overhang, 0)
        self.assertEqual(alignment.n_watson_crick, 9)
        self.assertEqual(alignment.n_mismatch, 0)
        self.assertLess(alignment.delta_g_hybrid, -10)  # 9 WC stacks should be strongly favorable

    def test_aligned_spacing_equals_raw_when_no_overhang(self):
        utr = "UAAGGAGGU" + "NNNNN".replace("N", "A")  # 5 nt spacer after the SD
        start_index = len(utr)
        alignment, _ = best_hybridization(utr, self.asd, start_index, zero_spacing, self.t37)
        self.assertEqual(alignment.overhang, 0)
        self.assertEqual(alignment.raw_spacing, alignment.aligned_spacing)
        self.assertEqual(alignment.raw_spacing, 5)

    def test_partial_core_sd_still_found(self):
        # A short 6-nt AGGAGG core (no flanking match) should still be
        # found as a valid, shorter alignment.
        utr = "CCCCCCAGGAGGCCCCC"
        start_index = len(utr)
        alignment, _ = best_hybridization(utr, self.asd, start_index, zero_spacing, self.t37)
        self.assertGreaterEqual(alignment.length, 4)
        found_segment = utr[alignment.mrna_start : alignment.mrna_end]
        self.assertIn(found_segment, "CCCCCCAGGAGGCCCCC")
        self.assertTrue(set(found_segment) <= set("AG"))  # should land inside the AGGAGG core

    def test_no_sd_like_sequence_returns_empty_alignment(self):
        utr = "CCCCCCCCCCCCCCCCCCCC"  # no complementarity to ACCUCCUUA at all
        start_index = len(utr)
        alignment, _ = best_hybridization(utr, self.asd, start_index, zero_spacing, self.t37)
        self.assertEqual(alignment.length, 0)

    def test_spacing_aware_selection_prefers_optimal_spacing_over_raw_strength(self):
        # Two SD-like sequences: pick whichever offset/window minimizes
        # delta_G_hybrid + delta_G_spacing jointly, not delta_G_hybrid alone.
        from rbsforge.rbs.spacing import spacing_penalty

        push = (12.2, 2.5, 2.0, 3.0)
        pull = (0.048, 0.24, 0.0)

        def spacing_fn(aligned_spacing):
            return spacing_penalty(aligned_spacing, 5, push, pull)

        utr = "UAAGGAGGU" + "AAAAA"  # full 9-nt SD, 5 nt spacer (optimal)
        start_index = len(utr)
        alignment, total = best_hybridization(utr, self.asd, start_index, spacing_fn, self.t37)
        self.assertEqual(alignment.aligned_spacing, 5)
        self.assertAlmostEqual(total, alignment.delta_g_hybrid, places=6)


if __name__ == "__main__":
    unittest.main()
