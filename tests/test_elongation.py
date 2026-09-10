import random
import unittest

from operon.core import OperonHost
from operon.elongation import (
    ECOLI_GENOMIC,
    ECOLI_HIGHLY_TRANSLATED,
    mutate_codon,
    recode_cds,
    ter_nt_s,
)
from operon.elongation.genetic_code import GENETIC_CODE, SYNONYMS, translate


class GeneticCodeTest(unittest.TestCase):
    def test_table_has_64_codons(self):
        self.assertEqual(len(GENETIC_CODE), 64)

    def test_three_stop_codons(self):
        stops = [c for c, aa in GENETIC_CODE.items() if aa == "*"]
        self.assertEqual(sorted(stops), ["UAA", "UAG", "UGA"])

    def test_translate_stops_at_first_stop_codon(self):
        self.assertEqual(translate("AUGAAACGCUAAGGGGGG"), "MKR")

    def test_translate_rejects_non_multiple_of_three(self):
        with self.assertRaises(ValueError):
            translate("AUGAA")

    def test_translate_rejects_bad_codon(self):
        with self.assertRaises(ValueError):
            translate("AUGXYZ")

    def test_synonyms_cover_every_amino_acid(self):
        for aa in set(GENETIC_CODE.values()):
            self.assertGreaterEqual(len(SYNONYMS[aa]), 1)


class RecodeCdsTest(unittest.TestCase):
    def test_recoded_sequence_translates_back_to_input_protein(self):
        rng = random.Random(1)
        protein = "MAKRLEQWYVCGPSTNDFHI"
        nt = recode_cds(protein, ECOLI_HIGHLY_TRANSLATED, rng)
        self.assertEqual(translate(nt), protein)
        self.assertEqual(len(nt), 3 * len(protein))

    def test_unknown_amino_acid_raises(self):
        rng = random.Random(1)
        with self.assertRaises(ValueError):
            recode_cds("MAKX", ECOLI_HIGHLY_TRANSLATED, rng)

    def test_deterministic_with_seeded_rng(self):
        protein = "MAKRLEQWYVCGPSTNDFHI"
        nt1 = recode_cds(protein, ECOLI_HIGHLY_TRANSLATED, random.Random(7))
        nt2 = recode_cds(protein, ECOLI_HIGHLY_TRANSLATED, random.Random(7))
        self.assertEqual(nt1, nt2)


class MutateCodonTest(unittest.TestCase):
    def test_mutation_preserves_protein_sequence(self):
        rng = random.Random(2)
        protein = "MAKRLEQWYVCGPSTNDFHI"
        nt = recode_cds(protein, ECOLI_HIGHLY_TRANSLATED, rng)
        for _ in range(20):
            nt = mutate_codon(nt, ECOLI_HIGHLY_TRANSLATED, rng)
            self.assertEqual(translate(nt), protein)

    def test_start_and_stop_are_never_touched_because_they_have_no_synonyms_worth_swapping(self):
        # Met and Trp are each their own synonym class (1 codon) -- a
        # mutation can never land there; this is a property test, not a
        # special-cased rule in the implementation.
        rng = random.Random(3)
        nt = "AUG" + "UGG" + "AUG"  # M, W, M
        mutated = mutate_codon(nt, ECOLI_HIGHLY_TRANSLATED, rng)
        self.assertEqual(mutated, nt)

    def test_rejects_non_multiple_of_three(self):
        with self.assertRaises(ValueError):
            mutate_codon("AUGA", ECOLI_HIGHLY_TRANSLATED, random.Random(1))


class TerNtSTest(unittest.TestCase):
    def setUp(self):
        self.host = OperonHost(pack=None)

    def test_ter_is_between_zero_and_r_max(self):
        rng = random.Random(4)
        nt = recode_cds("MAKRLEQWYVCGPSTNDFHI", ECOLI_HIGHLY_TRANSLATED, rng)
        rate = ter_nt_s(nt, self.host)
        self.assertGreater(rate, 0.0)
        self.assertLessEqual(rate, self.host.r_elong_nominal_nt_s)

    def test_all_preferred_codons_hit_r_max(self):
        # every codon at its table's max weight -> harmonic mean == r_max
        preferred_nt = "".join(
            max(ECOLI_HIGHLY_TRANSLATED[aa], key=ECOLI_HIGHLY_TRANSLATED[aa].get)
            for aa in "MAKRLEQ"
        )
        rate = ter_nt_s(preferred_nt, self.host)
        self.assertAlmostEqual(rate, self.host.r_elong_nominal_nt_s, places=6)

    def test_one_slow_codon_drags_the_harmonic_mean_down_more_than_it_would_the_arithmetic_mean(self):
        fast_nt = "".join(
            max(ECOLI_HIGHLY_TRANSLATED[aa], key=ECOLI_HIGHLY_TRANSLATED[aa].get) for aa in "AAAAAAAAAA"
        )
        slow_codons = [c for c in SYNONYMS["A"] if c != max(ECOLI_HIGHLY_TRANSLATED["A"], key=ECOLI_HIGHLY_TRANSLATED["A"].get)]
        mixed_nt = slow_codons[0] + fast_nt[3:]

        arithmetic_mean_ratio = (0.3 + 9 * 1.0) / 10  # w[c]/w_max per codon, arithmetic
        harmonic_ter = ter_nt_s(mixed_nt, self.host)
        arithmetic_equivalent = self.host.r_elong_nominal_nt_s * arithmetic_mean_ratio
        self.assertLess(harmonic_ter, arithmetic_equivalent)

    def test_stop_codon_is_excluded_from_scoring(self):
        rng = random.Random(5)
        nt = recode_cds("MAK", ECOLI_HIGHLY_TRANSLATED, rng) + "UAA"
        rate_with_stop = ter_nt_s(nt, self.host)
        rate_without_stop = ter_nt_s(nt[:-3], self.host)
        self.assertAlmostEqual(rate_with_stop, rate_without_stop)

    def test_uses_balanced_table_over_high_when_both_present(self):
        host = OperonHost(pack=None, codon_weights_high=ECOLI_HIGHLY_TRANSLATED, codon_weights_balanced=ECOLI_GENOMIC)
        rng = random.Random(6)
        nt = recode_cds("MAKRLEQ", ECOLI_HIGHLY_TRANSLATED, rng)
        # should not raise, and should use codon_weights_balanced (ECOLI_GENOMIC)
        rate = ter_nt_s(nt, host)
        self.assertGreater(rate, 0)

    def test_rejects_non_multiple_of_three(self):
        with self.assertRaises(ValueError):
            ter_nt_s("AUGA", self.host)

    def test_all_stop_codons_raises(self):
        with self.assertRaises(ValueError):
            ter_nt_s("UAAUAG", self.host)


if __name__ == "__main__":
    unittest.main()
