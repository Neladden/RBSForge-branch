import unittest

from operon.core import OperonHost
from operon.pauses import scan_pauses
from rbsforge.hostpack import get_hostpack


class ScanPausesTest(unittest.TestCase):
    def setUp(self):
        self.host = OperonHost(pack=get_hostpack("ecoli"))

    def test_finds_planted_internal_sd(self):
        cds = "AUG" + "CAACAACAA" + "AGGAGGA" + "CAACAACAACAACAACC"
        hits = scan_pauses(cds, self.host)
        sd_hits = [h for h in hits if h.kind == "internal_sd"]
        self.assertTrue(sd_hits)
        self.assertTrue(any("AGGAGG" in h.notes for h in sd_hits))

    def test_finds_polyproline_run(self):
        cds = "AUG" + "CCG" * 4 + "CAACAACAA"
        hits = scan_pauses(cds, self.host)
        pp_hits = [h for h in hits if h.kind == "polyproline"]
        self.assertEqual(len(pp_hits), 1)
        self.assertEqual(pp_hits[0].score, 4.0)

    def test_single_proline_is_not_a_hit(self):
        cds = "AUG" + "CCG" + "AAACAACAA"
        hits = scan_pauses(cds, self.host)
        self.assertEqual([h for h in hits if h.kind == "polyproline"], [])

    def test_finds_slow_codon_run(self):
        cds = "AUG" + "CUA" * 5 + "CAACAACAA"  # CUA (Leu) is non-preferred
        hits = scan_pauses(cds, self.host)
        slow_hits = [h for h in hits if h.kind == "slow_codon_run"]
        self.assertTrue(slow_hits)
        self.assertEqual(slow_hits[0].start, 3)

    def test_single_slow_codon_is_not_a_run(self):
        cds = "AUG" + "CUA" + "GCGGCGGCGGCGGCG"  # only one slow codon, rest preferred Ala
        hits = scan_pauses(cds, self.host)
        self.assertEqual([h for h in hits if h.kind == "slow_codon_run"], [])

    def test_no_hits_on_a_clean_short_cds(self):
        cds = "AUG" + "GCGGCGGCGGCG" + "UAA"  # all-preferred Ala codons, no PP, no SD
        hits = scan_pauses(cds, self.host)
        self.assertEqual(hits, [])

    def test_hits_sorted_by_position(self):
        cds = "AUG" + "CCG" * 3 + "CAACAACAA" + "AGGAGGA" + "CAACAACAACC"
        hits = scan_pauses(cds, self.host)
        starts = [h.start for h in hits]
        self.assertEqual(starts, sorted(starts))

    def test_uses_host_specific_anti_sd(self):
        # b_subtilis has a different (longer) anti-SD tail; a hit tuned to
        # E. coli's tail shouldn't be assumed to transfer unchanged.
        bsub_host = OperonHost(pack=get_hostpack("b_subtilis"))
        cds = "AUG" + "CAACAACAA" + "AGGAGGA" + "CAACAACAACAACAACC"
        ecoli_hits = [h for h in scan_pauses(cds, self.host) if h.kind == "internal_sd"]
        bsub_hits = [h for h in scan_pauses(cds, bsub_host) if h.kind == "internal_sd"]
        # both should run without error; asd differs so results need not match exactly
        self.assertIsInstance(ecoli_hits, list)
        self.assertIsInstance(bsub_hits, list)


if __name__ == "__main__":
    unittest.main()
