import unittest

from operon.promoter_calculator import predict, scan_promoters
from operon.promoter_calculator.calculator import _revcomp

# A consensus sigma70 promoter: -35 (TTGACA), an in-range 17 bp spacer,
# -10 (TATAAT), a 6 nt discriminator, and a 20 nt ITR, with flanking
# padding so a full UP element and downstream ITR both fit.
CONSENSUS_PROMOTER = (
    "A" * 30 + "TTGACA" + "A" * 17 + "TATAAT" + "C" * 6 + "G" * 20 + "A" * 10
)
# Same layout, hexamers replaced with sequence unrelated to either consensus motif.
WEAK_PROMOTER = (
    "A" * 30 + "GGGGGG" + "A" * 17 + "CCCCCC" + "C" * 6 + "G" * 20 + "A" * 10
)


class PromoterCalculatorEndToEndTest(unittest.TestCase):
    def test_consensus_hexamers_score_stronger_than_scrambled(self):
        """spec section 20, implementation-order row 10: a canonical
        TATAAT...17 bp...TTGACA construct should be flagged as a strong TSS."""
        strong = predict(CONSENSUS_PROMOTER, organism="ecoli").best
        weak = predict(WEAK_PROMOTER, organism="ecoli").best
        self.assertIsNotNone(strong)
        self.assertIsNotNone(weak)
        self.assertGreater(strong.tx_rate, weak.tx_rate)
        self.assertLess(strong.dg_total, weak.dg_total)

    def test_breakdown_terms_sum_to_dg_total(self):
        best = predict(CONSENSUS_PROMOTER, organism="ecoli").best
        self.assertAlmostEqual(sum(best.breakdown.terms().values()), best.dg_total, places=9)

    def test_element_spans_recover_the_hexamers(self):
        best = predict(CONSENSUS_PROMOTER, organism="ecoli").best
        seq = CONSENSUS_PROMOTER
        s, e = best.hex35_span
        self.assertEqual(seq[s:e], best.hex35)
        s, e = best.hex10_span
        self.assertEqual(seq[s:e], best.hex10)
        self.assertEqual(best.hex35, "TTGACA")
        self.assertEqual(best.hex10, "TATAAT")

    def test_reverse_strand_hit_has_identical_energy_to_forward(self):
        """Embedding a promoter as the reverse complement of a construct
        must score identically to scanning it on the forward strand --
        the coordinate-flip bookkeeping must not perturb the energy."""
        forward_best = predict(CONSENSUS_PROMOTER, organism="ecoli").best

        construct = "C" * 40 + _revcomp(CONSENSUS_PROMOTER) + "C" * 40
        result = predict(construct, organism="ecoli")
        reverse_best = result.best

        self.assertEqual(reverse_best.strand, -1)
        self.assertAlmostEqual(reverse_best.dg_total, forward_best.dg_total, places=9)
        self.assertAlmostEqual(reverse_best.tx_rate, forward_best.tx_rate, places=6)

    def test_scan_promoters_excludes_the_intended_tss(self):
        result = predict(CONSENSUS_PROMOTER, organism="ecoli")
        intended = result.best.tss

        all_hits = scan_promoters(CONSENSUS_PROMOTER, intended_tss=None)
        cryptic = scan_promoters(CONSENSUS_PROMOTER, intended_tss=intended, tau_tx=0.0)

        self.assertTrue(any(h.tss == intended for h in all_hits))
        self.assertFalse(any(h.tss == intended and h.strand == 1 for h in cryptic))

    def test_no_hits_when_sequence_too_short(self):
        result = predict("ACGT" * 5, organism="ecoli")
        self.assertIsNone(result.best)
        self.assertEqual(result.forward, {})
        self.assertEqual(result.reverse, {})


if __name__ == "__main__":
    unittest.main()
