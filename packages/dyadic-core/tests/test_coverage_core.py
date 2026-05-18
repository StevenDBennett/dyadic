"""Additional core tests for coverage."""

import unittest

from dyadic_core import DualNumber, TwoAdicProcessor


class TestDualNumberCoords(unittest.TestCase):
    def test_coords_returns_triple(self):
        d = DualNumber(42, k=16)
        coords = d.coords()
        self.assertEqual(len(coords), 3)
        v, alpha, e = coords
        self.assertEqual(v, d.v)
        self.assertEqual(alpha, d.alpha)
        self.assertEqual(e, d.e)


class TestProcessorDlog(unittest.TestCase):
    def test_dlog_returns_triple(self):
        proc = TwoAdicProcessor(16)
        d = DualNumber(42, 16)
        result = proc.dlog(d)
        self.assertEqual(len(result), 3)
        self.assertEqual(result, d.coords())


if __name__ == "__main__":
    unittest.main()
