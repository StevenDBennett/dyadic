"""Additional core tests for coverage."""

import unittest

from dyadic_core import DualNumber, TwoAdicProcessor, dlog_bootstrap


class TestDualNumberCoords(unittest.TestCase):
    def test_coords_returns_triple(self):
        d = DualNumber(42, k=16)
        coords = d.coords()
        self.assertEqual(len(coords), 3)
        v, alpha, e = coords
        self.assertEqual(v, d.v)
        self.assertEqual(alpha, d.alpha)
        self.assertEqual(e, d.e)


class TestFromCoords(unittest.TestCase):
    def test_from_coords_roundtrip(self):
        d = DualNumber.from_coords(v=0, alpha=0, e=42, k=16)
        self.assertTrue(d.verify())
        self.assertEqual(d.value, pow(5, 42, 1 << 16))

    def test_from_coords_zero(self):
        d = DualNumber.from_coords(v=None, alpha=0, e=0, k=16)
        self.assertTrue(d.is_zero)

    def test_from_coords_v_above_k(self):
        d = DualNumber.from_coords(v=20, alpha=0, e=0, k=16)
        self.assertTrue(d.is_zero)

    def test_from_coords_alpha_normalised(self):
        d = DualNumber.from_coords(v=0, alpha=3, e=0, k=16)
        self.assertEqual(d.alpha, 1)

    def test_from_coords_e_negative(self):
        d = DualNumber.from_coords(v=0, alpha=0, e=-1, k=16)
        self.assertTrue(d.verify())


class TestProcessorPow(unittest.TestCase):
    def test_pow_zero(self):
        proc = TwoAdicProcessor(16)
        d = DualNumber(3, k=16)
        result = proc.pow(d, 0)
        self.assertEqual(result.value, 1)

    def test_pow_one(self):
        proc = TwoAdicProcessor(16)
        d = DualNumber(3, k=16)
        result = proc.pow(d, 1)
        self.assertEqual(result.value, 3)

    def test_pow_negative(self):
        proc = TwoAdicProcessor(16)
        d = DualNumber(3, k=16)
        inv = proc.inv(d)
        result = proc.pow(d, -1)
        self.assertEqual(result.value, inv.value)


class TestDlogBootstrap(unittest.TestCase):
    def test_bootstrap_identity(self):
        for k in [4, 8, 12]:
            e = dlog_bootstrap(5, k)
            self.assertEqual(pow(5, e, 1 << k), 5 % (1 << k))

    def test_bootstrap_known(self):
        e = dlog_bootstrap(5, 16)
        self.assertEqual(e, 1)

    def test_bootstrap_small_k(self):
        self.assertEqual(dlog_bootstrap(1, 2), 0)


class TestProcessorDlog(unittest.TestCase):
    def test_dlog_returns_triple(self):
        proc = TwoAdicProcessor(16)
        d = DualNumber(42, 16)
        result = proc.dlog(d)
        self.assertEqual(len(result), 3)
        self.assertEqual(result, d.coords())


if __name__ == "__main__":
    unittest.main()
