import unittest
import warnings

from rbsforge.rbs.calculator import RBSCalculator
from operon.assembly import assemble
from operon.core import CDS, Operon, OperonHost
from operon.coupling import LeaderlessCouplingError, coupled_tirs, k_reinitiation, to_physical


def make_two_cistron_operon(policy, rbs2="AGGAGGAAAAA"):
    cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAUAA")
    cds2 = CDS(id="geneB", aa_or_nt="AUGCCCAAACCCAAACCCAAACCCAAACCCUAA")
    return Operon(
        promoter_dna="G" * 10,
        cds_list=[cds1, cds2],
        rbs_list=["AGGAGGAAAAA", rbs2],
        terminator_dna="T" * 10,
        host=OperonHost(pack=None),
        intergenic_policy=policy,
    )


class KReinitiationTest(unittest.TestCase):
    def test_published_fixed_points(self):
        self.assertAlmostEqual(k_reinitiation(-4), 0.022)
        self.assertAlmostEqual(k_reinitiation(0), 0.0072)
        self.assertAlmostEqual(k_reinitiation(25), 0.0072)
        self.assertAlmostEqual(k_reinitiation(-25), 0.0072 / 11.6)

    def test_beyond_25_decays(self):
        self.assertLess(k_reinitiation(125), k_reinitiation(25))

    def test_interpolated_region_is_monotonic_between_endpoints(self):
        lo = k_reinitiation(-25)
        hi = k_reinitiation(-5)
        mid = k_reinitiation(-15)
        self.assertLess(lo, mid)
        self.assertLess(mid, hi)


class ToPhysicalTest(unittest.TestCase):
    def test_is_documented_identity_seam(self):
        # No published au->s^-1 formula exists (spec section 7.2): this
        # must stay identity until a real reporter-library refit says
        # otherwise -- don't let a future edit silently invent one.
        self.assertEqual(to_physical(12345.0, OperonHost(pack=None)), 12345.0)


class CoupledTirsTest(unittest.TestCase):
    def setUp(self):
        self.calc = RBSCalculator.for_species("ecoli", temperature_c=37.0, use_vienna_if_available=False)
        self.host = OperonHost(pack=self.calc.hostpack)

    def test_refuses_non_vienna_folder_by_default(self):
        a = assemble(make_two_cistron_operon("overlap-4"))
        with self.assertRaises(RuntimeError):
            coupled_tirs(a, self.host, self.calc)

    def test_warns_and_proceeds_when_vienna_not_required(self):
        a = assemble(make_two_cistron_operon("overlap-4"))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            r = coupled_tirs(a, self.host, self.calc, require_vienna=False)
        self.assertEqual(len(caught), 1)
        self.assertEqual(len(r), 2)
        self.assertTrue(all(isinstance(x, float) for x in r))

    def test_first_cds_is_uncoupled_mono_cistronic_tir(self):
        a = assemble(make_two_cistron_operon("overlap-4"))
        r = coupled_tirs(a, self.host, self.calc, require_vienna=False)
        mono = self.calc.predict_one(a.mrna, a.starts[0])
        self.assertAlmostEqual(r[0], mono.v1_style_rate)

    def test_d_minus_4_gives_larger_reinitiation_contribution_than_insulation(self):
        a_auga = assemble(make_two_cistron_operon("overlap-4"))
        a_insulated = assemble(make_two_cistron_operon("overlap-25"))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r_auga = coupled_tirs(a_auga, self.host, self.calc, require_vienna=False)
            r_insulated = coupled_tirs(a_insulated, self.host, self.calc, require_vienna=False)

        r0 = r_auga[0]
        reinit_auga = self.host.k_p * k_reinitiation(a_auga.junctions[0].d_nt) * r0
        reinit_insulated = self.host.k_p * k_reinitiation(a_insulated.junctions[0].d_nt) * r0
        self.assertGreater(reinit_auga, reinit_insulated)
        # d=-4's much larger re-init leak should show up in the coupled TIR too.
        self.assertGreater(r_auga[1], r_insulated[1])

    def test_leaderless_first_cds_raises(self):
        leaderless_cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCAUUAGCACCACC")
        cds2 = CDS(id="geneB", aa_or_nt="AUGCCCCCCCCCUAA")
        operon = Operon(
            promoter_dna="C" * 30,  # no SD-like sequence anywhere -> leaderless
            cds_list=[leaderless_cds1, cds2],
            rbs_list=["", ""],
            terminator_dna="",
            host=OperonHost(pack=None),
            intergenic_policy="abut",
        )
        a = assemble(operon)
        with self.assertRaises(LeaderlessCouplingError):
            coupled_tirs(a, self.host, self.calc, require_vienna=False)

    def test_leaderless_downstream_cds_keeps_only_reinitiation(self):
        # CDS1 is long enough (>200 nt, no SD-like content of its own) that
        # RBS1's "AGGAGG" falls outside predict_one's 200-nt max search
        # window for CDS2's start -- CDS2 genuinely has no SD candidate.
        cds1 = CDS(id="geneA", aa_or_nt="AUG" + "CAA" * 70 + "UAA")
        cds2 = CDS(id="geneB", aa_or_nt="AUGCCCCCCCCCUAA")
        operon = Operon(
            promoter_dna="G" * 10,
            cds_list=[cds1, cds2],
            rbs_list=["AGGAGGAAAAA", ""],
            terminator_dna="T" * 10,
            host=OperonHost(pack=None),
            intergenic_policy="abut",
        )
        a = assemble(operon)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = coupled_tirs(a, self.host, self.calc, require_vienna=False)
        downstream = self.calc.predict_one(a.mrna, a.starts[1])
        self.assertTrue(downstream.leaderless)
        expected_reinit = self.host.k_p * k_reinitiation(a.junctions[0].d_nt) * r[0]
        self.assertAlmostEqual(r[1], expected_reinit)

    def test_single_cds_operon_has_no_coupling_terms(self):
        cds1 = CDS(id="geneA", aa_or_nt="AUGAAACGCAUUAGCACCACCAUUACCACCACCAUCACCAUUACCACAUAA")
        operon = Operon(
            promoter_dna="G" * 10, cds_list=[cds1], rbs_list=["AGGAGGAAAAA"],
            terminator_dna="T" * 10, host=OperonHost(pack=None), intergenic_policy="free",
        )
        a = assemble(operon)
        r = coupled_tirs(a, self.host, self.calc, require_vienna=False)
        self.assertEqual(len(r), 1)


if __name__ == "__main__":
    unittest.main()
