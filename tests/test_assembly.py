import unittest

from operon.assembly import assemble
from operon.core import CDS, Operon, OperonHost


def make_operon(cds_list, rbs_list, promoter="GGGGGGGGGG", terminator="TTTTTTTTTT",
                 policy="free", tss=None):
    return Operon(
        promoter_dna=promoter,
        cds_list=cds_list,
        rbs_list=rbs_list,
        terminator_dna=terminator,
        host=OperonHost(pack=None),
        intergenic_policy=policy,
        tss=tss,
    )


class SingleCdsAssemblyTest(unittest.TestCase):
    def test_layout_and_features(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAACGCTAA")
        operon = make_operon([cds1], ["AGGAGGAAAA"])
        a = assemble(operon)

        self.assertEqual(a.dna, "GGGGGGGGGG" + "AGGAGGAAAA" + "ATGAAACGCTAA" + "TTTTTTTTTT")
        kinds = {f["id"]: f for f in a.features}
        self.assertEqual(kinds["promoter"], {"id": "promoter", "type": "promoter", "start": 0, "end": 10, "strand": 1})
        self.assertEqual(kinds["RBS_geneA"]["type"], "rbs")
        self.assertEqual(kinds["geneA"]["type"], "CDS")
        self.assertEqual(kinds["terminator"]["type"], "terminator")
        self.assertEqual(a.junctions, [])

    def test_mrna_starts_at_default_tss_and_is_rna(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAACGCTAA")
        operon = make_operon([cds1], ["AGGAGGAAAA"])
        a = assemble(operon)
        # default TSS = end of promoter_dna
        self.assertEqual(a.mrna, "AGGAGGAAAA" + "AUGAAACGCUAA" + "UUUUUUUUUU")
        self.assertNotIn("T", a.mrna)

    def test_explicit_tss_inside_promoter(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAA")
        operon = make_operon([cds1], ["AGGAGG"], promoter="CCCCCGGGGG", tss=5)
        a = assemble(operon)
        self.assertEqual(a.mrna, a.dna[5:].upper().replace("T", "U"))

    def test_starts_and_cds_end_are_mrna_relative(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAACGCTAA")
        operon = make_operon([cds1], ["AGGAGGAAAA"])
        a = assemble(operon)
        tss = len("GGGGGGGGGG")
        expected_start = len("GGGGGGGGGG" + "AGGAGGAAAA") - tss
        self.assertEqual(a.starts, [expected_start])
        self.assertEqual(a.mrna[a.starts[0] : a.starts[0] + 3], "AUG")
        self.assertEqual(a.cds_end, [expected_start + len("ATGAAACGCTAA")])

    def test_empty_cds_list_rejected(self):
        operon = make_operon([], [])
        with self.assertRaises(ValueError):
            assemble(operon)

    def test_unknown_policy_rejected(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATG")
        operon = make_operon([cds1], [""], policy="bogus")
        with self.assertRaises(ValueError):
            assemble(operon)


class JunctionPolicyTest(unittest.TestCase):
    def test_abut_gives_d_zero(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAATAA")
        cds2 = CDS(id="geneB", aa_or_nt="ATGCCCTAA")
        operon = make_operon([cds1, cds2], ["AGGAGG", ""], policy="abut")
        a = assemble(operon)
        self.assertEqual(len(a.junctions), 1)
        self.assertEqual(a.junctions[0].d_nt, 0)
        self.assertEqual(a.junctions[0].between, ("geneA", "geneB"))
        # no gap: geneB's CDS sequence immediately follows geneA's
        self.assertIn("ATGAAATAA" + "ATGCCCTAA", a.dna)

    def test_spacer_n_gives_d_equal_n(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAATAA")
        cds2 = CDS(id="geneB", aa_or_nt="ATGCCCTAA")
        operon = make_operon([cds1, cds2], ["AGGAGG", "AAAAAAAAAA"], policy="spacer-10")
        a = assemble(operon)
        self.assertEqual(a.junctions[0].d_nt, 10)

    def test_spacer_wrong_length_rejected(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAATAA")
        cds2 = CDS(id="geneB", aa_or_nt="ATGCCCTAA")
        operon = make_operon([cds1, cds2], ["AGGAGG", "AAAA"], policy="spacer-10")
        with self.assertRaises(ValueError):
            assemble(operon)

    def test_free_policy_d_equals_rbs_length(self):
        cds1 = CDS(id="geneA", aa_or_nt="ATGAAATAA")
        cds2 = CDS(id="geneB", aa_or_nt="ATGCCCTAA")
        operon = make_operon([cds1, cds2], ["AGGAGG", "AAAAA"], policy="free")
        a = assemble(operon)
        self.assertEqual(a.junctions[0].d_nt, 5)

    def test_overlap_4_gives_d_minus_4_and_canonical_auga_at_the_junction(self):
        # spec section 3.1: canonical 4 nt overlap -> d = -4. CDS_i's own
        # leading 4 nt (here "AUGA...") become the physically shared span.
        cds1 = CDS(id="geneA", aa_or_nt="AAAAAAAAAA")
        cds2 = CDS(id="geneB", aa_or_nt="ATGACCCCCC")
        operon = make_operon([cds1, cds2], ["AGGAGG", ""], policy="overlap-4")
        a = assemble(operon)
        j = a.junctions[0]
        self.assertEqual(j.d_nt, -4)
        self.assertEqual(a.dna[j.start_i_plus_1 : j.start_i_plus_1 + 4], "ATGA")
        # the last 4 nt of geneA's own sequence were superseded by geneB's lead-in
        self.assertEqual(j.stop_i - j.start_i_plus_1, 4)

    def test_overlap_25_is_insulating_distance(self):
        cds1 = CDS(id="geneA", aa_or_nt="A" * 40)
        cds2 = CDS(id="geneB", aa_or_nt="ATG" + "C" * 30)
        operon = make_operon([cds1, cds2], ["AGGAGG", ""], policy="overlap-25")
        a = assemble(operon)
        self.assertEqual(a.junctions[0].d_nt, -25)

    def test_overlap_longer_than_downstream_cds_rejected(self):
        cds1 = CDS(id="geneA", aa_or_nt="A" * 40)
        cds2 = CDS(id="geneB", aa_or_nt="ATG")
        operon = make_operon([cds1, cds2], ["AGGAGG", ""], policy="overlap-25")
        with self.assertRaises(ValueError):
            assemble(operon)

    def test_three_cistron_operon_has_two_junctions(self):
        cds = [CDS(id=f"gene{i}", aa_or_nt="ATG" + "C" * 9 + "TAA") for i in range(3)]
        operon = make_operon(cds, ["AGGAGG", "", ""], policy="abut")
        a = assemble(operon)
        self.assertEqual(len(a.junctions), 2)
        self.assertEqual([j.between for j in a.junctions], [("gene0", "gene1"), ("gene1", "gene2")])
        self.assertEqual(len(a.starts), 3)
        self.assertEqual(len(a.cds_end), 3)


if __name__ == "__main__":
    unittest.main()
