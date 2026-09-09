import unittest

from rbsforge import predict
from rbsforge.rbs.calculator import RBSCalculator


class PredictEndToEndTest(unittest.TestCase):
    def setUp(self):
        # Strong SD ("AGGAGG"), 6-nt raw spacer, then AUG.
        self.strong_seq = "AGGAGGACAACTAAATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT"
        # No SD-like sequence at all before the AUG.
        self.leaderless_seq = "CCCCCCCCCCCCCCCCCCCCAUGAAACGCAUUAGCACCACC"
        # Weak/short spacer variant of the strong sequence (SD right next to start codon).
        self.compressed_seq = "AGGAGGAUGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT"

    def test_finds_start_codon_and_scores_it(self):
        result = predict(self.strong_seq, species="ecoli", temperature_c=37.0)
        self.assertIsNotNone(result.best)
        self.assertEqual(result.best.codon, "AUG")
        self.assertIsNotNone(result.best.delta_g_total)
        self.assertGreater(result.best.translation_initiation_rate, 0)

    def test_breakdown_terms_sum_to_total(self):
        result = predict(self.strong_seq, species="ecoli", temperature_c=37.0)
        best = result.best
        self.assertAlmostEqual(sum(best.breakdown.terms.values()), best.delta_g_total, places=6)
        self.assertEqual(
            set(best.breakdown.terms),
            {"standby", "mRNA:rRNA", "spacing", "start", "stacking", "mRNA"},
        )

    def test_leaderless_start_is_skipped_from_ranking(self):
        result = predict(self.leaderless_seq, species="ecoli", temperature_c=37.0)
        self.assertTrue(all(r.leaderless for r in result.results))
        self.assertIsNone(result.best)

    def test_compressed_spacing_scores_worse_than_optimal(self):
        strong = predict(self.strong_seq, species="ecoli", temperature_c=37.0).best
        compressed = predict(self.compressed_seq, species="ecoli", temperature_c=37.0).best
        self.assertGreater(
            compressed.breakdown.terms["spacing"], strong.breakdown.terms["spacing"]
        )

    def test_temperature_changes_prediction(self):
        low_t = predict(self.strong_seq, species="ecoli", temperature_c=25.0).best
        high_t = predict(self.strong_seq, species="ecoli", temperature_c=70.0).best
        self.assertNotAlmostEqual(low_t.delta_g_total, high_t.delta_g_total, places=2)
        # SD:aSD duplex should weaken (less negative) at the higher temperature.
        self.assertGreater(
            high_t.breakdown.terms["mRNA:rRNA"], low_t.breakdown.terms["mRNA:rRNA"]
        )

    def test_species_changes_prediction(self):
        ecoli = predict(self.strong_seq, species="ecoli", temperature_c=37.0).best
        bsub = predict(self.strong_seq, species="b_subtilis", temperature_c=37.0).best
        self.assertNotAlmostEqual(ecoli.delta_g_total, bsub.delta_g_total, places=2)

    def test_custom_hostpack_via_calculator(self):
        from rbsforge.hostpack import HostPack

        custom = HostPack(
            name="test organism",
            phylogeny="unknown",
            gram_stain="unknown",
            t_growth_c=45.0,
            asd="ACCUCCUUA",
        )
        calc = RBSCalculator(custom, temperature_c=45.0)
        result = calc.predict(self.strong_seq)
        self.assertIsNotNone(result.best)
        self.assertEqual(result.temperature_c, 45.0)


if __name__ == "__main__":
    unittest.main()
