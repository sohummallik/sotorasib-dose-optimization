"""Reject malformed aggregate counts before they can produce plausible estimates."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from grouped_exposure_response import fit_grouped


class GroupedInputValidationTests(unittest.TestCase):
    def setUp(self):
        self.data = pd.read_csv(ROOT / 'data/fda_figure26_quartiles.csv')
        self.coordinates = self.data.auc_label.to_numpy(dtype=float)

    def test_impossible_response_counts_are_rejected(self):
        # These three errors previously returned numerical odds ratios.
        for value in [-1, 25.5, 58, np.nan, np.inf]:
            with self.subTest(responders=value):
                data = self.data.astype({'responders': float})
                data.loc[0, 'responders'] = value
                with self.assertRaises(ValueError):
                    fit_grouped(data, self.coordinates)

    def test_invalid_totals_are_rejected(self):
        for value in [0, -1, 57.5, np.nan, np.inf]:
            with self.subTest(total=value):
                data = self.data.astype({'total': float})
                data.loc[0, 'total'] = value
                with self.assertRaises(ValueError):
                    fit_grouped(data, self.coordinates)

    def test_nonfinite_nonpositive_or_misaligned_coordinates_are_rejected(self):
        for value in [0, -1, np.nan, np.inf]:
            coordinates = self.coordinates.copy()
            coordinates[0] = value
            with self.subTest(coordinate=value), self.assertRaises(ValueError):
                fit_grouped(self.data, coordinates)
        for coordinates in [self.coordinates[:-1], self.coordinates.reshape(4, 1)]:
            with self.subTest(shape=coordinates.shape), self.assertRaises(ValueError):
                fit_grouped(self.data, coordinates)

    def test_slope_requires_distinct_coordinates(self):
        with self.assertRaisesRegex(ValueError, 'distinct'):
            fit_grouped(self.data, np.repeat(13000, 4))

    def test_missing_count_columns_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'columns'):
            fit_grouped(self.data.drop(columns='total'), self.coordinates)

    def test_complete_separation_cannot_report_a_finite_odds_ratio(self):
        data = self.data.assign(responders=[0, 0, 57, 57])
        with self.assertRaisesRegex(ValueError, 'separated'):
            fit_grouped(data, self.coordinates)


if __name__ == '__main__':
    unittest.main()
