import unittest

from rbsforge.hostpack import BUILTIN_HOSTPACKS, get_hostpack


class HostPackTest(unittest.TestCase):
    def test_ecoli_lookup_by_alias(self):
        for alias in ("ecoli", "E. coli", "E_COLI_K12", "e-coli"):
            hp = get_hostpack(alias)
            self.assertEqual(hp.name, "Escherichia coli K-12")

    def test_unknown_species_raises_with_available_list(self):
        with self.assertRaises(KeyError) as ctx:
            get_hostpack("nonexistent_organism")
        self.assertIn("ecoli", str(ctx.exception))

    def test_all_builtin_packs_have_normalized_rna_asd(self):
        for hp in BUILTIN_HOSTPACKS.values():
            self.assertNotIn("T", hp.asd)
            self.assertTrue(set(hp.asd) <= set("ACGU"))
            self.assertGreaterEqual(len(hp.asd), 9)

    def test_gram_positive_hosts_have_wider_spacing_optimum(self):
        ecoli = get_hostpack("ecoli")
        bsub = get_hostpack("b_subtilis")
        self.assertEqual(ecoli.gram_stain, "negative")
        self.assertEqual(bsub.gram_stain, "positive")
        self.assertGreater(bsub.s_opt, ecoli.s_opt)


if __name__ == "__main__":
    unittest.main()
