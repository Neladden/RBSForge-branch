import unittest

from rbsforge.rbs.spacing import spacing_penalty

PUSH = (12.2, 2.5, 2.0, 3.0)
PULL = (0.048, 0.24, 0.0)


class SpacingPenaltyTest(unittest.TestCase):
    def test_zero_at_optimum(self):
        self.assertAlmostEqual(spacing_penalty(5, 5, PUSH, PULL), 0.0, places=6)

    def test_continuity_at_optimum(self):
        # push (ds -> 0-) and pull (ds == 0) should agree closely at the boundary
        just_below = spacing_penalty(4, 5, PUSH, PULL)
        at_opt = spacing_penalty(5, 5, PUSH, PULL)
        self.assertLess(just_below, 0.1)
        self.assertLess(at_opt, 0.1)

    def test_stretch_is_quadratic_pull(self):
        # ds = 3 -> 0.048*9 + 0.24*3 = 0.432 + 0.72 = 1.152
        self.assertAlmostEqual(spacing_penalty(8, 5, PUSH, PULL), 1.152, places=3)

    def test_compression_grows_as_spacing_shrinks(self):
        p_minus1 = spacing_penalty(4, 5, PUSH, PULL)
        p_minus3 = spacing_penalty(2, 5, PUSH, PULL)
        p_minus6 = spacing_penalty(-1, 5, PUSH, PULL)
        self.assertLess(p_minus1, p_minus3)
        self.assertLess(p_minus3, p_minus6)

    def test_compression_penalty_steeper_than_extension(self):
        # Equal-magnitude deviation: compression should cost more than extension
        # once away from the near-zero boundary region.
        compressed = spacing_penalty(5 - 6, 5, PUSH, PULL)
        extended = spacing_penalty(5 + 6, 5, PUSH, PULL)
        self.assertGreater(compressed, extended)

    def test_custom_s_opt(self):
        self.assertAlmostEqual(spacing_penalty(8, 8, PUSH, PULL), 0.0, places=6)


if __name__ == "__main__":
    unittest.main()
