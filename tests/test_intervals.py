from __future__ import annotations

import unittest

from workflow.lib.intervals import FeatureInterval, PointIntervalIndex


class IntervalTests(unittest.TestCase):
    def test_strand_and_closed_coordinates(self):
        index = PointIntervalIndex(
            [
                FeatureInterval("2L", 10, 20, "+", "a"),
                FeatureInterval("2L", 20, 30, "+", "b"),
                FeatureInterval("2L", 10, 20, "-", "c"),
            ],
            bin_size=8,
        )
        self.assertEqual({x.feature_id for x in index.query("2L", 20, "+")}, {"a", "b"})
        self.assertEqual({x.feature_id for x in index.query("2L", 10, "-")}, {"c"})
        self.assertEqual(index.query("2L", 9, "+"), [])


if __name__ == "__main__":
    unittest.main()

