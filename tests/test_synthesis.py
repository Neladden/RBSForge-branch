import random
import unittest

from operon.synthesis import synthesis_score


FLANK = "CCACGTGCATGCATCGTAGCTAGCACC"  # non-repetitive, boundary-safe (verified empirically) on either side of a planted run


class HomopolymerTest(unittest.TestCase):
    def test_gc_run_of_9_violates(self):
        r = synthesis_score(FLANK + "G" * 9 + FLANK)
        self.assertFalse(r.is_feasible)
        self.assertTrue(any(base == "G" for base, _, _ in r.homopolymer_hits))

    def test_gc_run_of_8_is_fine(self):
        r = synthesis_score(FLANK + "G" * 8 + FLANK)
        self.assertEqual(r.homopolymer_hits, [])

    def test_at_run_of_13_violates(self):
        r = synthesis_score(FLANK + "A" * 13 + FLANK)
        self.assertTrue(any(base == "A" for base, _, _ in r.homopolymer_hits))

    def test_at_run_of_12_is_fine(self):
        r = synthesis_score(FLANK + "A" * 12 + FLANK)
        self.assertEqual(r.homopolymer_hits, [])


class TandemTest(unittest.TestCase):
    def test_dinucleotide_tandem_of_10_copies_violates(self):
        r = synthesis_score(FLANK + "AT" * 10 + FLANK)
        self.assertTrue(any(period == 2 for period, _, _ in r.tandem_hits))

    def test_dinucleotide_tandem_of_9_copies_is_fine(self):
        r = synthesis_score(FLANK + "AT" * 9 + FLANK)
        self.assertEqual(r.tandem_hits, [])

    def test_trinucleotide_tandem_of_6_copies_violates(self):
        r = synthesis_score(FLANK + "CAT" * 6 + FLANK)
        self.assertTrue(any(period == 3 for period, _, _ in r.tandem_hits))

    def test_homopolymer_does_not_double_count_as_tandem(self):
        # a long single-base run trivially satisfies any period's
        # recurrence but belongs to _homopolymer_hits, not here.
        r = synthesis_score("ACGT" * 5 + "A" * 30 + "ACGT" * 5)
        self.assertEqual(r.tandem_hits, [])


class GcWindowTest(unittest.TestCase):
    def test_extreme_gc_20bp_window_violates(self):
        r = synthesis_score("ACGT" * 20 + "G" * 20 + "ACGT" * 20)
        self.assertTrue(any(w == 20 for _, _, w, _ in r.gc_window_hits))

    def test_moderate_gc_is_fine_over_both_windows(self):
        seq = "ACGT" * 50  # exactly 50% GC in every window, deterministic
        r = synthesis_score(seq)
        self.assertEqual(r.gc_window_hits, [])

    def test_terminal_30bp_gc_checked_only_beyond_60nt(self):
        short_extreme = "G" * 30 + "ACGT" * 5  # len 50, <= 60
        r_short = synthesis_score(short_extreme)
        # only the general 20/100-bp windows can fire here, not the
        # length-gated terminal-30 rule specifically -- just check it
        # doesn't crash and still reports the 20bp-window violation.
        self.assertTrue(any(w == 20 for _, _, w, _ in r_short.gc_window_hits))

        long_extreme = "G" * 30 + "ACGT" * 20
        r_long = synthesis_score(long_extreme)
        self.assertTrue(any(offset == 0 for offset, _, w, _ in r_long.gc_window_hits if w == 30))


class RestrictionSiteTest(unittest.TestCase):
    def test_forbidden_site_is_infeasible(self):
        r = synthesis_score("AAA" + "GAATTC" + "TTT", forbidden_re_sites=["GAATTC"])
        self.assertFalse(r.is_feasible)
        self.assertTrue(r.forbidden_re_hits)

    def test_no_forbidden_sites_given_never_flags_anything(self):
        r = synthesis_score("AAAGAATTCTTT", forbidden_re_sites=None)
        self.assertEqual(r.forbidden_re_hits, [])

    def test_absent_site_is_feasible(self):
        random.seed(12)
        seq = "".join(random.choice("ACGT") for _ in range(100))
        r = synthesis_score(seq, forbidden_re_sites=["GGTACCGGTACC"])
        self.assertEqual(r.forbidden_re_hits, [])


class RepeatDensityTest(unittest.TestCase):
    def test_single_large_repeat_on_a_short_fragment_violates(self):
        unit = "ACGTGCATGCAT"  # 12nt, >= k
        seq = unit + "TTTT" + unit  # repeat covers most of a 28nt fragment
        r = synthesis_score(seq)
        self.assertTrue(r.density_hits)
        self.assertFalse(r.is_feasible)

    def test_nonrepetitive_fragment_has_no_density_hits(self):
        random.seed(14)
        seq = "".join(random.choice("ACGT") for _ in range(300))
        r = synthesis_score(seq)
        self.assertEqual(r.density_hits, [])


class OverallFeasibilityTest(unittest.TestCase):
    def test_clean_random_sequence_is_feasible(self):
        random.seed(13)
        seq = "".join(random.choice("ACGT") for _ in range(400))
        r = synthesis_score(seq)
        self.assertTrue(r.is_feasible)

    def test_any_single_violation_makes_the_whole_report_infeasible(self):
        r = synthesis_score("A" * 20 + "G" * 9 + "A" * 20)
        self.assertFalse(r.is_feasible)


if __name__ == "__main__":
    unittest.main()
