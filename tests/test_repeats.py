import unittest

from operon.repeats import _revcomp, find_repeats, scan_motifs


class DirectRepeatTest(unittest.TestCase):
    def test_plants_a_direct_12mer_twice(self):
        unit = "ACGTACGTACGT"
        seq = "A" * 5 + unit + "C" * 20 + unit + "A" * 5
        r = find_repeats(seq, k=12)
        hit = next(h for h in r.direct if h.unit == unit)
        self.assertEqual(hit.positions, [5, 5 + len(unit) + 20])

    def test_no_direct_repeat_in_nonrepetitive_sequence(self):
        seq = "ACGTGCATGCATGCTAGCTAGGCATCGATCGATGCTAGCATCG"
        r = find_repeats(seq, k=12)
        self.assertEqual(r.direct, [])

    def test_too_short_sequence_returns_empty_report(self):
        r = find_repeats("ACGT", k=12)
        self.assertEqual(r.all_hits(), [])


class InvertedRepeatTest(unittest.TestCase):
    def test_plants_a_non_palindromic_inverted_repeat(self):
        unit = "AAAAACCCCCTT"
        self.assertNotEqual(unit, _revcomp(unit))
        seq = unit + "G" * 30 + _revcomp(unit)
        r = find_repeats(seq, k=12)
        hit = next(h for h in r.inverted if h.unit == unit)
        self.assertEqual(hit.positions, [0])
        self.assertEqual(hit.rc_positions, [len(unit) + 30])

    def test_self_palindromic_kmer_is_not_reported_as_inverted(self):
        # a k-mer that equals its own reverse complement doesn't imply a
        # second, distinct locus by itself.
        palindrome = "AAACCCGGGTTT"
        self.assertEqual(palindrome, _revcomp(palindrome))
        r = find_repeats(palindrome + "G" * 20, k=12)
        self.assertEqual(r.inverted, [])

    def test_each_inverted_pair_reported_once_not_twice(self):
        unit = "AAAAACCCCCTT"
        seq = unit + "G" * 30 + _revcomp(unit)
        r = find_repeats(seq, k=12)
        units_seen = {h.unit for h in r.inverted}
        self.assertEqual(len(r.inverted), 1)
        self.assertTrue(units_seen == {unit} or units_seen == {_revcomp(unit)})


class TandemRepeatTest(unittest.TestCase):
    def test_trinucleotide_tandem_run_detected(self):
        # flank with a base absent from the CAG unit so the run can't be
        # phased into the flank by coincidence (e.g. a trailing G would
        # complete a "GCA" rotation of the same period-3 pattern).
        seq = "T" * 10 + "CAG" * 10 + "T" * 10
        r = find_repeats(seq, k=12)
        hit = next(h for h in r.tandem if h.period == 3)
        self.assertEqual(hit.positions, [10])
        self.assertEqual(hit.length, 30)
        self.assertGreaterEqual(hit.copies, 10)

    def test_short_tandem_below_k_is_not_reported_as_general_tandem(self):
        seq = "G" * 10 + "CAG" * 2 + "G" * 10  # span 6 < k=12
        r = find_repeats(seq, k=12)
        self.assertEqual(r.tandem, [])

    def test_homopolymer_is_a_period_1_tandem(self):
        seq = "C" * 20
        r = find_repeats(seq + "A" * 20, k=12)
        hit = next(h for h in r.tandem if h.period == 1)
        self.assertEqual(hit.unit, "C")
        self.assertEqual(hit.length, 20)


class TerminalRepeatTest(unittest.TestCase):
    def test_tandem_at_sequence_start_is_terminal(self):
        seq = "A" * 5 + "GATCGATCGATCGATC"
        r = find_repeats(seq, k=12)
        self.assertTrue(any(h.positions == [0] for h in r.terminal))

    def test_tandem_at_sequence_end_is_terminal(self):
        seq = "GATCGATCGATCGATC" + "A" * 5
        r = find_repeats(seq, k=12)
        hit = next(h for h in r.terminal)
        self.assertEqual(hit.positions[0] + hit.length, len(seq))

    def test_interior_tandem_below_terminal_span_is_not_terminal(self):
        nonrepetitive_flank = "ACGTGCATGCAT"
        seq = nonrepetitive_flank + "CA" * 2 + nonrepetitive_flank[::-1]  # span 4 < 5, interior
        r = find_repeats(seq, k=12)
        self.assertEqual(r.terminal, [])


class LongestAndAllHitsTest(unittest.TestCase):
    def test_longest_reflects_the_biggest_hit(self):
        seq = "G" * 10 + "CAG" * 20 + "G" * 10
        r = find_repeats(seq, k=12)
        self.assertEqual(r.longest(), max(h.length for h in r.all_hits()))

    def test_empty_report_longest_is_zero(self):
        r = find_repeats("ACGT", k=12)
        self.assertEqual(r.longest(), 0)


class ScanMotifsTest(unittest.TestCase):
    def test_finds_exact_forward_and_reverse_hit(self):
        motif = "GCTAGC"
        seq = "AAA" + motif + "TTT"
        hits = scan_motifs(seq, [motif])
        fwd = [h for h in hits if h["strand"] == 1]
        self.assertEqual(fwd, [{"motif": motif, "position": 3, "strand": 1, "mismatches": 0}])

    def test_no_hit_for_absent_motif(self):
        hits = scan_motifs("AAAAAAAAAAAAAAAAAAAA", ["GCTAGC"])
        self.assertEqual(hits, [])

    def test_one_mismatch_tolerance(self):
        motif = "GCTAGC"
        seq = "AAA" + "GCTACC" + "TTT"  # 1 mismatch from GCTAGC
        exact_hits = scan_motifs(seq, [motif], max_mismatches=0)
        fuzzy_hits = scan_motifs(seq, [motif], max_mismatches=1)
        self.assertEqual([h for h in exact_hits if h["strand"] == 1], [])
        self.assertTrue(any(h["strand"] == 1 and h["mismatches"] == 1 for h in fuzzy_hits))

    def test_motif_longer_than_sequence_is_skipped_not_an_error(self):
        hits = scan_motifs("ACGT", ["ACGTACGTACGTACGT"])
        self.assertEqual(hits, [])


if __name__ == "__main__":
    unittest.main()
