"""Reject unsupported covariates and nonphysical PK parameter sets."""
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from model import individual_params, load_params, typical_subject


class PKInputValidationTests(unittest.TestCase):
    def setUp(self):
        self.parameters = load_params()

    def subject(self, **changes):
        return dict(typical_subject(), **changes)

    def test_unrecognized_categories_do_not_fall_back_to_reference(self):
        for change in [{"ecog": 99}, {"ecog": True}, {"sex": "female"},
                       {"tumor": "small"}, {"race": "unknown"}, {"race": None}]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                individual_params(self.parameters, self.subject(**change))

    def test_indicator_strings_and_numbers_are_not_booleans(self):
        for name in ["ppi", "highfat"]:
            for value in ["False", "True", 0, 1, None]:
                with self.subTest(name=name, value=value), self.assertRaisesRegex(ValueError, "boolean"):
                    individual_params(self.parameters, self.subject(**{name: value}))

    def test_albumin_is_finite_positive_and_cannot_make_clearance_negative(self):
        for value in [0, -1, np.nan, np.inf, True, "40"]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "alb_gL"):
                individual_params(self.parameters, self.subject(alb_gL=value))
        # A positive albumin value may still make the linear multiplier negative.
        with self.assertRaisesRegex(ValueError, "Calculated PK parameters"):
            individual_params(self.parameters, self.subject(alb_gL=1))

    def test_random_effect_shape_finiteness_and_overflow(self):
        for eta in [1, [0, 0], [[0, 0, 0]], [0, np.nan, 0], [np.inf, 0, 0]]:
            with self.subTest(eta=eta), self.assertRaisesRegex(ValueError, "eta"):
                individual_params(self.parameters, self.subject(), eta=eta)
        for value in [-1000, 1000]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "Calculated PK parameters"):
                individual_params(self.parameters, self.subject(), eta=[0, value, 0])

    def test_relative_bioavailability_overrides_are_positive_not_bounded_by_one(self):
        for name in ["f1bs_override", "f1ss_override"]:
            for value in [0, -1, np.nan, np.inf, True, "1"]:
                with self.subTest(name=name, value=value), self.assertRaisesRegex(ValueError, name):
                    individual_params(self.parameters, self.subject(), **{name: value})
        result = individual_params(self.parameters, self.subject(),
                                   f1bs_override=1.2, f1ss_override=1.1)
        self.assertEqual(result["F1BS"], 1.2)
        self.assertEqual(result["F1SS"], 1.1)

    def test_nonphysical_structural_inputs_are_rejected(self):
        for name, value in [("CLSS", -1), ("V2", 0), ("K24", np.inf), ("KIND_F", np.nan)]:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Calculated PK parameters"):
                individual_params(dict(self.parameters, **{name: value}), self.subject())

    def test_all_published_contrast_inputs_remain_accepted(self):
        contrasts = json.loads((ROOT / "data/covariate_targets.json").read_text())["contrasts"]
        for contrast in contrasts:
            with self.subTest(contrast=contrast["contrast"]):
                result = individual_params(self.parameters, self.subject(**contrast["overrides"]))
                self.assertTrue(all(np.isfinite(value) and value > 0 for value in result.values()))
        default = individual_params(self.parameters, self.subject())
        explicit = individual_params(self.parameters, self.subject(ppi=np.bool_(False)), eta=np.zeros(3))
        self.assertEqual(default, explicit)


if __name__ == "__main__":
    unittest.main()
