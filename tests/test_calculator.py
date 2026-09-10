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


class PredictOneFootprintUnfoldingTest(unittest.TestCase):
    """Spec section 6 / Appendix B: predict_one's two-fold bound-state
    unfolding, replacing the old blanket -1*unconstrained-MFE mRNA term."""

    def setUp(self):
        self.calc = RBSCalculator.for_species("ecoli", temperature_c=37.0, use_vienna_if_available=False)
        self.utr = "CCCCCCCCCC" + "AGGAGG" + "AAAAA"
        self.no_structure_cds = "CAA" * 12
        # 12 nt hairpin (4 bp GC stem, 4 nt loop): entirely within
        # [start, start + 3 + footprint_cds) = offsets 0..15 relative to
        # start when placed immediately after the start codon (offset 3..14).
        self.footprint_hairpin = "GGGG" + "AAAA" + "CCCC"

    def test_no_structure_gives_zero_unfolding(self):
        mrna = self.utr + "AUG" + self.no_structure_cds
        r = self.calc.predict_one(mrna, mrna.index("AUG", 5))
        self.assertFalse(r.leaderless)
        self.assertAlmostEqual(r.breakdown.terms["mRNA"], 0.0, places=6)

    def test_hairpin_wholly_in_footprint_costs_its_full_stability(self):
        mrna = self.utr + "AUG" + self.footprint_hairpin + "CAAA" * 6
        r = self.calc.predict_one(mrna, mrna.index("AUG", 5))
        self.assertFalse(r.leaderless)
        # The ribosome's footprint covers exactly this hairpin, so it must
        # pay the hairpin's full stability to occupy that region -- a real,
        # positive cost, not zero and not negative.
        self.assertGreater(r.breakdown.terms["mRNA"], 0.0)

    def test_hairpin_outside_footprint_does_not_cost_via_this_term(self):
        # A hairpin far upstream of the SD, untouched by the footprint+start
        # occupied set, does not have to be unfolded for the ribosome to
        # bind at the start codon -- this term alone reports ~0 for it.
        # (A real SD-region hairpin's effect on initiation would have to
        # show up through the SD:aSD duplex search instead; that search is
        # sequence-only in this model and does not yet see structure -- see
        # rbsforge/docs/MODEL.md.)
        upstream_hairpin = "GGGG" + "UUUU" + "CCCC"
        mrna = upstream_hairpin + self.utr + "AUG" + self.no_structure_cds
        r = self.calc.predict_one(mrna, mrna.index("AUG", 5))
        self.assertFalse(r.leaderless)
        self.assertAlmostEqual(r.breakdown.terms["mRNA"], 0.0, places=6)

    def test_mixed_structure_window_hurts_less_than_old_blanket_penalty(self):
        # A window with a hairpin *outside* the footprint (upstream, near
        # the 5' end of the UTR) and a second hairpin *inside* it: the new
        # footprint-aware term must charge less than the old
        # -1*unconstrained-MFE approach would have, because the upstream
        # hairpin no longer has to be "paid for" at all.
        upstream_hairpin = "GGGG" + "UUUU" + "CCCC"
        mrna = (
            upstream_hairpin
            + self.utr
            + "AUG"
            + self.footprint_hairpin
            + "CAAA" * 6
        )
        start = mrna.index("AUG", 5)
        r = self.calc.predict_one(mrna, start)
        self.assertFalse(r.leaderless)

        hp = self.calc.hostpack
        pre_len = min(start, start, 200)
        window_start = start - pre_len
        window_end = min(len(mrna), start + 3 + hp.cutoff_post)
        window = mrna[window_start:window_end]
        old_blanket_term = -self.calc.folder.fold(window).delta_g

        self.assertLess(r.breakdown.terms["mRNA"], old_blanket_term)

    def test_footprint_span_reported(self):
        mrna = self.utr + "AUG" + self.no_structure_cds
        start = mrna.index("AUG", 5)
        r = self.calc.predict_one(mrna, start)
        hp = self.calc.hostpack
        self.assertEqual(r.footprint_span, (start, start + hp.footprint_cds))

    def test_extra_unpaired_global_forces_additional_positions_open(self):
        # The coupling hook: positions outside the default footprint can
        # still be forced unpaired via extra_unpaired_global, and doing so
        # can only raise (never lower) the bound-state energy, so the
        # reported unfolding cost can only grow.
        mrna = self.utr + "AUG" + self.footprint_hairpin + "CAAA" * 6
        start = mrna.index("AUG", 5)
        baseline = self.calc.predict_one(mrna, start)
        forced_extra = self.calc.predict_one(
            mrna, start, extra_unpaired_global=range(0, start)
        )
        self.assertFalse(baseline.leaderless)
        self.assertFalse(forced_extra.leaderless)
        self.assertGreaterEqual(
            forced_extra.breakdown.terms["mRNA"], baseline.breakdown.terms["mRNA"]
        )

    def test_predict_matches_predict_one_for_every_start(self):
        mrna = self.utr + "AUG" + self.footprint_hairpin + "CAAA" * 6
        result = self.calc.predict(mrna)
        for r in result.results:
            if r.leaderless:
                continue
            r2 = self.calc.predict_one(mrna, r.index)
            self.assertAlmostEqual(r.delta_g_total, r2.delta_g_total, places=9)


if __name__ == "__main__":
    unittest.main()
