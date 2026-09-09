import unittest

from rbsforge.constants import celsius_to_kelvin
from rbsforge.thermo.nn_params import hairpin_loop_energy, watson_crick_stack


class TemperatureScalingTest(unittest.TestCase):
    def test_stack_energy_at_37_matches_table(self):
        # GC/CG delta_G_37 tabulated as -3.42 kcal/mol.
        dg = watson_crick_stack("GC/CG", celsius_to_kelvin(37.0))
        self.assertAlmostEqual(dg, -3.42, places=2)

    def test_duplex_weakens_at_higher_temperature(self):
        # model summary section 18: helices should become LESS stable
        # (less negative delta_G) as temperature rises.
        dg_37 = watson_crick_stack("GC/CG", celsius_to_kelvin(37.0))
        dg_60 = watson_crick_stack("GC/CG", celsius_to_kelvin(60.0))
        dg_80 = watson_crick_stack("GC/CG", celsius_to_kelvin(80.0))
        self.assertGreater(dg_60, dg_37)
        self.assertGreater(dg_80, dg_60)

    def test_au_rich_stack_weakens_faster_than_gc_rich(self):
        # A weaker (less enthalpy-dense) AU-rich stack should lose a
        # smaller absolute amount of stability than a GC-rich one when
        # heated the same amount, matching the qualitative claim in
        # model summary section 12.1/18 that GC-rich duplexes are more
        # heat-resistant -- but check it via the fractional (relative)
        # destabilization, which should be *larger* for the weaker duplex.
        au_37 = watson_crick_stack("AU/UA", celsius_to_kelvin(37.0))
        au_70 = watson_crick_stack("AU/UA", celsius_to_kelvin(70.0))
        gc_37 = watson_crick_stack("GC/CG", celsius_to_kelvin(37.0))
        gc_70 = watson_crick_stack("GC/CG", celsius_to_kelvin(70.0))
        au_fraction_remaining = au_70 / au_37
        gc_fraction_remaining = gc_70 / gc_37
        self.assertLess(au_fraction_remaining, gc_fraction_remaining)

    def test_hairpin_loop_extrapolation_beyond_table(self):
        dg9 = hairpin_loop_energy(9, celsius_to_kelvin(37.0))
        dg20 = hairpin_loop_energy(20, celsius_to_kelvin(37.0))
        self.assertGreater(dg20, dg9)  # larger loops cost more (less favorable)

    def test_hairpin_too_small_is_impossible(self):
        self.assertEqual(hairpin_loop_energy(2, celsius_to_kelvin(37.0)), float("inf"))


if __name__ == "__main__":
    unittest.main()
