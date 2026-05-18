"""Tests for dyadic_core.exponent — ExponentSpace."""
import unittest

from dyadic_core.exponent import ExponentSpace
from dyadic_core import bitmask


class TestExponentSpace(unittest.TestCase):
    def test_lift_zero(self):
        es = ExponentSpace(5, 8)
        self.assertEqual(es.lift(0), 1)

    def test_lift_one(self):
        es = ExponentSpace(5, 8)
        self.assertEqual(es.lift(1), 5)

    def test_lift_matches_pow(self):
        es = ExponentSpace(5, 12)
        for e in range(8):
            self.assertEqual(es.lift(e), pow(5, e, 1 << 12))

    def test_periodicity(self):
        es = ExponentSpace(5, 8)
        N = 1 << 6
        self.assertEqual(es.lift(0), es.lift(N))

    def test_difference_operator(self):
        es = ExponentSpace(5, 8)
        f = lambda e: pow(5, e, 1 << 8)
        for e in range(4):
            diff = es.difference(f, e)
            expected = (f(e + 1) - f(e)) & bitmask(8)
            self.assertEqual(diff, expected)

    def test_integral(self):
        es = ExponentSpace(5, 8)
        f = lambda e: 1
        integral = es.integrate(f)
        self.assertEqual(integral, es.N)


if __name__ == "__main__":
    unittest.main()
