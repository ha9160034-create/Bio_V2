"""Tests for ΔScore uncertainty handling in the comparison UI."""

import unittest

from docking.interactions import interpret_score_difference


class TestScoreDifferenceUncertainty(unittest.TestCase):
    def test_sign_only_without_uncertainty(self):
        info = interpret_score_difference(-0.8)
        self.assertEqual(info.direction, "more favorable")
        self.assertNotIn("within the", str(info))

    def test_small_difference_without_uncertainty_uses_floor(self):
        info = interpret_score_difference(-0.1)
        self.assertIn("protocol noise floor", str(info))
        self.assertEqual(info.direction, "unchanged")

    def test_within_noise_band_claims_no_direction(self):
        info = interpret_score_difference(-0.8, uncertainty_kcal=1.0)
        self.assertIn("within the estimated", str(info))
        self.assertEqual(info.direction, "unchanged")

    def test_outside_noise_band_keeps_direction(self):
        info = interpret_score_difference(-0.8, uncertainty_kcal=0.2)
        self.assertEqual(info.direction, "more favorable")

    def test_noise_text_remains_non_clinical(self):
        info = interpret_score_difference(0.3, uncertainty_kcal=0.9)
        self.assertIn("not binding affinity", str(info))
        self.assertIn("clinical guidance", str(info))


if __name__ == "__main__":
    unittest.main()
