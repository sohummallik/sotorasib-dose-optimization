"""Checks for planning behavior and invalid inputs in the scenario analysis."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from design_sensitivity import sample_size_per_arm, precision_half_width


class DesignSensitivityTests(unittest.TestCase):
    def test_reversing_arm_labels_preserves_required_size(self):
        self.assertEqual(sample_size_per_arm(.35, .25, .05, .8),
                         sample_size_per_arm(.25, .35, .05, .8))

    def test_more_power_and_smaller_differences_require_more_patients(self):
        base = sample_size_per_arm(.35, .25, .05, .8)
        self.assertGreater(sample_size_per_arm(.35, .25, .05, .9), base)
        self.assertGreater(sample_size_per_arm(.30, .25, .05, .8), base)

    def test_four_times_the_sample_halves_approximate_interval_width(self):
        self.assertAlmostEqual(precision_half_width(.35, .25, 100) / 2,
                               precision_half_width(.35, .25, 400), places=12)

    def test_invalid_planning_inputs_are_rejected(self):
        for values in [(.25, .25, .05, .8), (1.1, .25, .05, .8),
                       (.35, .25, .05, .04), (float("nan"), .25, .05, .8)]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                sample_size_per_arm(*values)
        with self.assertRaises(ValueError):
            precision_half_width(.35, .25, 0)


if __name__ == "__main__":
    unittest.main()
