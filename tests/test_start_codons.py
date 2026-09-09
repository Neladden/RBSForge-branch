import unittest

from rbsforge.rbs.start_codons import find_start_codons, start_codon_energy
from rbsforge.hostpack import DEFAULT_START_CODON_DG


class StartCodonTest(unittest.TestCase):
    def test_finds_all_frames(self):
        # AUG at position 1 (frame-shifted relative to a would-be position 0)
        hits = find_start_codons("CAUGCC", allowed_codons=("AUG",))
        self.assertEqual(hits, [(1, "AUG")])

    def test_finds_multiple_codon_types(self):
        hits = find_start_codons("AUGGUGUUG", allowed_codons=("AUG", "GUG", "UUG"))
        self.assertEqual(hits, [(0, "AUG"), (3, "GUG"), (6, "UUG")])

    def test_dna_input_normalized(self):
        hits = find_start_codons("CATGCC", allowed_codons=("AUG",))
        self.assertEqual(hits, [(1, "AUG")])

    def test_energy_table_matches_published_values(self):
        self.assertAlmostEqual(start_codon_energy("AUG", DEFAULT_START_CODON_DG), -1.194)
        self.assertAlmostEqual(start_codon_energy("GUG", DEFAULT_START_CODON_DG), -0.0748)
        self.assertAlmostEqual(start_codon_energy("UUG", DEFAULT_START_CODON_DG), -0.0435)
        self.assertAlmostEqual(start_codon_energy("CUG", DEFAULT_START_CODON_DG), -0.03406)

    def test_aug_is_most_favorable(self):
        aug = start_codon_energy("AUG", DEFAULT_START_CODON_DG)
        gug = start_codon_energy("GUG", DEFAULT_START_CODON_DG)
        uug = start_codon_energy("UUG", DEFAULT_START_CODON_DG)
        self.assertLess(aug, gug)
        self.assertLess(gug, uug)

    def test_unrecognized_codon_falls_back_conservatively(self):
        fallback = start_codon_energy("AAA", DEFAULT_START_CODON_DG)
        self.assertGreater(fallback, max(DEFAULT_START_CODON_DG.values()))


if __name__ == "__main__":
    unittest.main()
