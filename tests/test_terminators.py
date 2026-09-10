import unittest

from operon.terminators import scan_intrinsic_terminators, scan_rho_terminators


def _make_intrinsic_terminator(flank_len=12):
    # GC stem, AAAA loop, poly-U tail -- a strong, unambiguous intrinsic terminator
    term = "GGGGGG" + "AAAA" + "CCCCCC" + "UUUUUUUU"
    return "A" * flank_len + term + "A" * flank_len, term


class IntrinsicTerminatorTest(unittest.TestCase):
    def test_finds_a_planted_strong_terminator_on_forward_strand(self):
        seq, _ = _make_intrinsic_terminator()
        hits = scan_intrinsic_terminators(seq)
        self.assertTrue(any(h.strand == 1 for h in hits))

    def test_hairpin_mfe_is_reported_and_stable(self):
        seq, _ = _make_intrinsic_terminator()
        hits = scan_intrinsic_terminators(seq)
        fwd = next(h for h in hits if h.strand == 1)
        self.assertLessEqual(fwd.score, -7.0)  # the crude MFE gate itself

    def test_no_hits_on_a_plain_alternating_sequence(self):
        plain = "ACGU" * 20
        self.assertEqual(scan_intrinsic_terminators(plain), [])

    def test_weak_stem_without_u_tract_is_not_a_hit_on_that_strand(self):
        # a real hairpin but no downstream U-tract on the forward strand
        # (note: A-only flanks would turn into a U-tract on the *reverse*
        # strand after reverse-complementing, so only assert about fwd).
        no_tail = "A" * 12 + "GGGGGG" + "AAAA" + "CCCCCC" + "A" * 12
        fwd_hits = [h for h in scan_intrinsic_terminators(no_tail) if h.strand == 1]
        self.assertEqual(fwd_hits, [])

    def test_annotated_span_excludes_the_intended_terminator(self):
        seq, _ = _make_intrinsic_terminator()
        all_hits = scan_intrinsic_terminators(seq)
        fwd = next(h for h in all_hits if h.strand == 1)
        filtered = scan_intrinsic_terminators(seq, annotated_span=(fwd.start, fwd.end), annotated_window=5)
        self.assertNotIn(fwd, filtered)
        self.assertLess(len(filtered), len(all_hits))

    def test_annotated_span_none_returns_everything(self):
        seq, _ = _make_intrinsic_terminator()
        self.assertEqual(
            len(scan_intrinsic_terminators(seq, annotated_span=None)),
            len(scan_intrinsic_terminators(seq)),
        )


class RhoTerminatorTest(unittest.TestCase):
    def _make_rut_region(self):
        window = list("A" * 78)
        for k in range(6):
            window[k * 13] = "C"
        window[5] = "G"
        return "".join(window)

    def test_finds_rut_candidate_with_downstream_pyrimidine_pause(self):
        rut = self._make_rut_region()
        seq = "A" * 20 + rut + "C" * 8 + "A" * 20
        hits = scan_rho_terminators(seq)
        self.assertTrue(any(h.strand == 1 for h in hits))

    def test_no_hit_without_a_downstream_pause(self):
        rut = self._make_rut_region()
        seq = "A" * 20 + rut + "A" * 120  # no pyrimidine tract, no hairpin, within search range
        hits = scan_rho_terminators(seq)
        self.assertEqual([h for h in hits if h.strand == 1 and h.start < 20 + len(rut)], [])

    def test_no_hit_when_g_equals_or_exceeds_c(self):
        # G's densely and evenly spread so *every* 78-nt sub-window (the
        # scanner doesn't only look at one fixed alignment) has at least
        # as many G's as the handful of C's -- C/G ratio never exceeds 1
        # anywhere, regardless of exact window placement.
        window = list("AG" * 60)  # 120 nt, alternating A/G: 60 G's throughout
        for k in range(6):
            window[k * 15] = "C"
        seq = "A" * 20 + "".join(window) + "C" * 8 + "A" * 20
        hits = scan_rho_terminators(seq)
        self.assertEqual(hits, [])

    def test_plain_sequence_has_no_rho_hits(self):
        plain = "ACGU" * 60
        self.assertEqual(scan_rho_terminators(plain), [])


if __name__ == "__main__":
    unittest.main()
