import unittest

from rbsforge.constants import celsius_to_kelvin
from rbsforge.rbs.standby import standby_penalty
from rbsforge.thermo.fold import BuiltinFolder


class StandbyPenaltyTest(unittest.TestCase):
    def setUp(self):
        self.folder = BuiltinFolder(celsius_to_kelvin(37.0))

    def test_no_room_gives_zero_penalty(self):
        penalty, start, end = standby_penalty("AGGAGGAAAUG", sd_duplex_start=0, standby_site_nt=4, folder=self.folder)
        self.assertEqual(penalty, 0.0)
        self.assertEqual((start, end), (0, 0))

    def test_penalty_is_never_negative(self):
        # A window with a hairpin that overlaps the standby region.
        window = "GGGGGCCCCCAAAAAGGGGGCCCCC"
        penalty, start, end = standby_penalty(window, sd_duplex_start=15, standby_site_nt=4, folder=self.folder)
        self.assertGreaterEqual(penalty, 0.0)
        self.assertEqual(end - start, 4)

    def test_structured_standby_region_costs_more_than_unstructured(self):
        structured = "GGGGCCCCAAAAAAAAAAAAAAAA"  # strong hairpin right where standby will sit
        unstructured = "AAAACCCCAAAAAAAAAAAAAAAA"  # same length, no complementarity
        p_structured, _, _ = standby_penalty(structured, sd_duplex_start=8, standby_site_nt=4, folder=self.folder)
        p_unstructured, _, _ = standby_penalty(unstructured, sd_duplex_start=8, standby_site_nt=4, folder=self.folder)
        self.assertGreaterEqual(p_structured, p_unstructured)


if __name__ == "__main__":
    unittest.main()
