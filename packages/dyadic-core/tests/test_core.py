"""Tests for dyadic_core — DualNumber, modinv, dlog, TwoAdicProcessor."""

import unittest

from dyadic_core import (
    DualNumber,
    TwoAdicProcessor,
    bitmask,
    dlog_residual_tracking,
    dual_add,
    modinv_newton,
    two_adic_dlog,
    two_adic_log5,
    valuation,
)


class TestMask(unittest.TestCase):
    def testbitmask_identity(self):
        for k in [1, 2, 3, 8, 16, 32]:
            self.assertEqual(bitmask(k), (1 << k) - 1)

    def testbitmask_range(self):
        for k in range(1, 33):
            m = bitmask(k)
            self.assertEqual(m.bit_length(), k)


class TestValuation(unittest.TestCase):
    def test_powers_of_two(self):
        for n in range(10):
            self.assertEqual(valuation(1 << n), n)

    def test_zero(self):
        self.assertIsNone(valuation(0))

    def test_odd(self):
        for n in [1, 3, 5, 7, 9, 255]:
            self.assertEqual(valuation(n), 0)

    def test_mixed(self):
        cases = {12: 2, 24: 3, 100: 2, 1024: 10, 1025: 0}
        for n, expected in cases.items():
            self.assertEqual(valuation(n), expected)


class TestModinvNewton(unittest.TestCase):
    def test_small_k(self):
        for k in [8, 16, 32]:
            for a in [1, 3, 5, 7, 9]:
                inv = modinv_newton(a, k)
                self.assertEqual((a * inv) & bitmask(k), 1)

    def test_random(self):
        import random

        for k in [8, 12, 16, 24]:
            for _ in range(20):
                a = random.randrange(1, 1 << k, 2)
                inv = modinv_newton(a, k)
                self.assertEqual((a * inv) & bitmask(k), 1)

    def test_even_raises(self):
        with self.assertRaises(ValueError):
            modinv_newton(2, 8)

    def test_k_zero_raises(self):
        with self.assertRaises(ValueError):
            modinv_newton(3, 0)


class TestTwoAdicLog5(unittest.TestCase):
    def test_suffix_stability(self):
        L8 = two_adic_log5(8)
        L16 = two_adic_log5(16)
        self.assertEqual(L16 & bitmask(8), L8)

    def test_caching(self):
        L1 = two_adic_log5(12)
        L2 = two_adic_log5(12)
        self.assertIs(L1, L2)


class TestTwoAdicDlog(unittest.TestCase):
    def test_known_values(self):
        for e_true in [0, 1, 2, 3, 7, 15, 31, 63]:
            k = 16
            a = pow(5, e_true, 1 << k)
            result = two_adic_dlog(a, k)
            self.assertIsNotNone(result)
            alpha, e = result
            if a & 3 == 1:
                self.assertEqual(e, e_true)
                self.assertEqual(alpha, 0)

    def test_odd_returns_tuple(self):
        for a in [1, 3, 5, 7, 255]:
            result = two_adic_dlog(a, 16)
            self.assertIsNotNone(result)
            alpha, e = result
            self.assertIn(alpha, (0, 1))

    def test_even_returns_none(self):
        for a in [0, 2, 4, 6, 8]:
            self.assertIsNone(two_adic_dlog(a, 16))

    def test_verify_roundtrip(self):
        for a in [1, 3, 5, 7, 255, 1023, 32767]:
            k = 16
            result = two_adic_dlog(a, k)
            if result is not None:
                alpha, e = result
                expected = a
                recomputed = pow(5, e, 1 << k)
                if alpha:
                    recomputed = (-recomputed) & bitmask(k)
                self.assertEqual(recomputed, expected)


class TestDualNumber(unittest.TestCase):
    def test_zero(self):
        d = DualNumber(0, 16)
        self.assertTrue(d.is_zero)
        self.assertEqual(d.value, 0)

    def test_one(self):
        d = DualNumber(1, 16)
        self.assertFalse(d.is_zero)
        self.assertEqual(d.value, 1)
        self.assertTrue(d.verify())

    def test_negative_one(self):
        d = DualNumber(-1 & bitmask(16), 16)
        self.assertTrue(d.verify())

    def test_power_of_two(self):
        d = DualNumber(8, 16)
        self.assertTrue(d.verify())
        self.assertEqual(d.v, 3)

    def test_overflow(self):
        d = DualNumber(1 << 17, 16)
        self.assertTrue(d.is_zero)

    def test_from_coords_roundtrip(self):
        for v in [0, 1, 2]:
            for alpha in [0, 1]:
                for e in [0, 1, 5, 10]:
                    k = 12
                    d = DualNumber.from_coords(v, alpha, e, k)
                    self.assertTrue(d.verify())
                    self.assertEqual(d.v, v)
                    self.assertEqual(d.alpha, alpha)
                    if not d.is_zero:
                        self.assertEqual(d.e, e % (1 << (k - 2)))

    def test_k_too_small_raises(self):
        with self.assertRaises(ValueError):
            DualNumber(1, 2)


class TestTwoAdicProcessor(unittest.TestCase):
    def setUp(self):
        self.proc = TwoAdicProcessor(16)

    def test_mul(self):
        a = DualNumber(3, 16)
        b = DualNumber(5, 16)
        c = self.proc.mul(a, b)
        self.assertEqual(c.value, 15)

    def test_mul_commutative(self):
        a = DualNumber(7, 16)
        b = DualNumber(11, 16)
        ab = self.proc.mul(a, b)
        ba = self.proc.mul(b, a)
        self.assertEqual(ab.value, ba.value)

    def test_mul_associative(self):
        a = DualNumber(3, 16)
        b = DualNumber(5, 16)
        c = DualNumber(7, 16)
        ab_c = self.proc.mul(self.proc.mul(a, b), c)
        a_bc = self.proc.mul(a, self.proc.mul(b, c))
        self.assertEqual(ab_c.value, a_bc.value)

    def test_inv(self):
        a = DualNumber(3, 16)
        inv_a = self.proc.inv(a)
        prod = self.proc.mul(a, inv_a)
        self.assertEqual(prod.value, 1)

    def test_inv_non_unit_raises(self):
        a = DualNumber(2, 16)
        with self.assertRaises(ValueError):
            self.proc.inv(a)

    def test_pow(self):
        a = DualNumber(3, 16)
        a4 = self.proc.pow(a, 4)
        self.assertEqual(a4.value, 81)

    def test_pow_zero(self):
        a = DualNumber(3, 16)
        a0 = self.proc.pow(a, 0)
        self.assertEqual(a0.value, 1)

    def test_pow_negative(self):
        a = DualNumber(3, 16)
        a_neg = self.proc.pow(a, -1)
        prod = self.proc.mul(a, a_neg)
        self.assertEqual(prod.value, 1)

    def test_mul_zero(self):
        a = DualNumber(3, 16)
        z = DualNumber(0, 16)
        c = self.proc.mul(a, z)
        self.assertTrue(c.is_zero)

    def test_overflow(self):
        a = DualNumber(1 << 14, 16)
        b = DualNumber(4, 16)
        c = self.proc.mul(a, b)
        self.assertTrue(c.is_zero)


class TestSynthesisFindings(unittest.TestCase):
    def test_exponential_isometry(self):
        k = 12
        for e in range(1, 32):
            diff = (pow(5, e, 1 << k) - 1) & bitmask(k)
            v2_diff = valuation(diff) if diff != 0 else k
            v2_e = valuation(e)
            self.assertEqual(v2_diff, v2_e + 2)

    def test_lut_bootstrap_correctness(self):
        k = 16
        for e_true in range(64):
            a = pow(5, e_true, 1 << k)
            result = two_adic_dlog(a, k)
            if result is not None:
                alpha, e = result
                if a & 3 == 1:
                    self.assertEqual(e, e_true)


class TestDlogResidualTracking(unittest.TestCase):
    def test_known_exponent(self):
        k = 10
        for e_true in [0, 1, 2, 3, 7, 15]:
            a = pow(5, e_true, 1 << k)
            if a & 3 == 1:
                e, history = dlog_residual_tracking(a, k)
                self.assertEqual(e, e_true)

    def test_history_nonempty(self):
        k = 10
        a = pow(5, 7, 1 << k)
        _, history = dlog_residual_tracking(a, k)
        self.assertGreater(len(history), 0)

    def test_history_keys(self):
        k = 10
        a = pow(5, 3, 1 << k)
        _, history = dlog_residual_tracking(a, k)
        required = {
            "bits",
            "eprec_before",
            "eprec_after",
            "tau_before",
            "v2_before",
            "tau_after",
            "v2_after",
            "delta",
        }
        for entry in history:
            self.assertTrue(required.issubset(entry.keys()))

    def test_v2_doubles(self):
        k = 12
        a = pow(5, 5, 1 << k)
        _, history = dlog_residual_tracking(a, k)
        for entry in history:
            if entry["v2_before"] is not None:
                self.assertGreaterEqual(
                    entry["v2_after"], 2 * entry["v2_before"], msg=f"v2 {entry}"
                )

    def test_small_k_returns_empty(self):
        e, history = dlog_residual_tracking(1, 2)
        self.assertEqual(e, 0)
        self.assertEqual(history, [])

    def test_invalid_a_raises(self):
        with self.assertRaises(ValueError):
            dlog_residual_tracking(3, 8)


class TestDualAdd(unittest.TestCase):
    """dual_add: each structural case verified against integer round-trip."""

    def _check(self, v_a, alpha_a, e_a, v_b, alpha_b, e_b, k):
        mask = bitmask(k)
        a = pow(5, e_a, 1 << k)
        if alpha_a:
            a = (-a) & mask
        a = (a << v_a) & mask
        b = pow(5, e_b, 1 << k)
        if alpha_b:
            b = (-b) & mask
        b = (b << v_b) & mask
        s = (a + b) & mask

        result = dual_add(v_a, alpha_a, e_a, v_b, alpha_b, e_b, k)

        if s == 0:
            self.assertIsNone(result[0], f"Expected cancellation, got {result}")
            return

        d = DualNumber(s, k)
        self.assertEqual(result, (d.v, d.alpha, d.e), f"{a} + {b} = {s}")

    def test_diff_valuation(self):
        for v in range(1, 5):
            self._check(0, 0, 0, v, 0, 0, 8)

    def test_same_sign_same_e(self):
        for e in range(4):
            self._check(0, 0, e, 0, 0, e, 8)

    def test_same_sign_same_e_alpha1(self):
        for e in range(4):
            self._check(0, 1, e, 0, 1, e, 8)

    def test_same_sign_diff_e(self):
        for e_a, e_b in [(0, 1), (1, 3), (0, 5), (2, 7)]:
            self._check(0, 0, e_a, 0, 0, e_b, 8)

    def test_opp_sign_diff_e(self):
        for e_a, e_b in [(0, 1), (1, 3), (0, 5), (2, 7)]:
            self._check(0, 0, e_a, 0, 1, e_b, 8)

    def test_exact_cancellation(self):
        for e in range(4):
            result = dual_add(0, 0, e, 0, 1, e, 8)
            self.assertIsNone(result[0])

    def test_different_v_with_signs(self):
        self._check(0, 1, 1, 2, 0, 3, 8)

    def test_k_edge(self):
        for k in [3, 4, 8, 16]:
            self._check(0, 0, 1, 0, 0, 3, k)
            self._check(0, 0, 1, 0, 1, 3, k)
            self._check(1, 0, 0, 3, 0, 0, k)


class TestDLogNewtonStep(unittest.TestCase):
    """Tests for dlog_newton_step and DLogNewtonStep."""

    def test_newton_step_refines_dlog(self):
        from dyadic_core import dlog_newton_step

        k = 16
        e_seed = 3
        a = pow(5, 7, 1 << k)
        log5_unit = 69  # two_adic_log5(16) >> 2  (pre-computed)
        new_e, _ = dlog_newton_step(
            e_seed,
            a,
            log5_unit,
            bits=10,
            mask=(1 << 10) - 1,
            emask=(1 << 8) - 1,
            new_eprec=8,
        )
        self.assertIsInstance(new_e, int)
        self.assertGreater(new_e, 0)


if __name__ == "__main__":
    unittest.main()
