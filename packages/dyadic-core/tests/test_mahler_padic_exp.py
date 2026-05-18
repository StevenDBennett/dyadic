"""
Tests for padic_exp, padic_log, g0 (core.py) and MahlerCalculus (mahler.py).
"""

import unittest

from dyadic_core import (
    MahlerCalculus,
    bitmask,
    g0,
    padic_exp,
    padic_log,
    two_adic_log5,
    valuation,
)


class TestPadicExp(unittest.TestCase):
    def test_exp_zero(self):
        self.assertEqual(padic_exp(0, 8), 1)

    def test_exp_minus4(self):
        for k in [8, 16, 32]:
            result = padic_exp(-4, k)
            self.assertEqual(result % 4, 1, f"exp(−4) should be ≡ 1 mod 4 at k={k}")

    def test_exp_log_roundtrip(self):
        for k in [8, 16]:
            mod = 1 << k
            for g in range(1, min(mod, 256), 4):
                log_g = padic_log(g, k)
                exp_log = padic_exp(log_g, k)
                self.assertEqual(exp_log, g % mod, f"exp(log({g})) ≠ {g} mod 2^{k}")

    def test_log_exp_roundtrip(self):
        for k in [8, 16]:
            mod = 1 << k
            for x in [4, 8, 12, 16, 20, 24, 28, 32]:
                if x >= mod:
                    continue
                exp_x = padic_exp(x, k)
                log_exp = padic_log(exp_x, k)
                self.assertEqual(log_exp, x % mod, f"log(exp({x})) ≠ {x} mod 2^{k}")

    def test_exp_requires_v2_ge_2(self):
        with self.assertRaises(ValueError):
            padic_exp(2, 8)
        with self.assertRaises(ValueError):
            padic_exp(1, 8)

    def test_log_requires_v2_ge_2(self):
        with self.assertRaises(ValueError):
            padic_log(3, 8)
        with self.assertRaises(ValueError):
            padic_log(2, 8)

    def test_exp_additive_property(self):
        k = 16
        mod = 1 << k
        x, y = 4, 8
        exp_x = padic_exp(x, k)
        exp_y = padic_exp(y, k)
        exp_xy = padic_exp(x + y, k)
        self.assertEqual(exp_xy, (exp_x * exp_y) % mod)

    def test_exp_at_k32(self):
        g0_16 = padic_exp(-4, 16)
        g0_32 = padic_exp(-4, 32)
        self.assertEqual(g0_32 & bitmask(16), g0_16)

    def test_two_adic_log5_ground_truth(self):
        # Regression test: padic_log(5, 10) fixed after duplicate-exp2_neg4 bug.
        # The correct value 636 was verified via series computation.
        from dyadic_core import two_adic_log5

        val = two_adic_log5(10)
        self.assertEqual(val, 636)
        # Suffix stability: higher-k values truncate to the same lower bits
        for k in [12, 16, 20]:
            higher = two_adic_log5(k)
            self.assertEqual(higher & bitmask(10), val)


class TestG0(unittest.TestCase):
    def test_g0_mod_4_is_1(self):
        for k in [8, 16, 32]:
            self.assertEqual(g0(k) % 4, 1)

    def test_g0_log_is_minus4(self):
        for k in [8, 16, 32]:
            mod = 1 << k
            g = g0(k)
            log_g = padic_log(g, k)
            expected = (-4) % mod
            self.assertEqual(log_g, expected, f"log(g₀) ≠ −4 mod 2^{k}")

    def test_g0_agrees_with_neg123_to_13_bits(self):
        k = 32
        g = g0(k)
        neg123 = (-123) % (1 << k)
        diff = abs(g - neg123) % (1 << k)
        self.assertEqual(valuation(diff), 13)


class TestMahlerPolynomial(unittest.TestCase):
    def test_C_x_0_is_1(self):
        self.assertEqual(MahlerCalculus.mahler_polynomial(0, 5), 1)

    def test_C_x_1_is_x(self):
        for x in range(10):
            self.assertEqual(MahlerCalculus.mahler_polynomial(1, x), x)

    def test_C_x_2_is_x_x_minus_1_over_2(self):
        for x in range(10):
            expected = x * (x - 1) // 2
            self.assertEqual(MahlerCalculus.mahler_polynomial(2, x), expected)

    def test_negative_n_is_zero(self):
        self.assertEqual(MahlerCalculus.mahler_polynomial(-1, 5), 0)


class TestMahlerTransform(unittest.TestCase):
    def test_roundtrip_constant_function(self):
        f = [7, 7, 7, 7, 7]
        coeffs = MahlerCalculus.to_mahler(f, 4)
        for x in range(5):
            self.assertEqual(MahlerCalculus.from_mahler(coeffs, x), f[x])

    def test_roundtrip_linear_function(self):
        f = [0, 1, 2, 3, 4, 5]
        coeffs = MahlerCalculus.to_mahler(f, 5)
        for x in range(6):
            self.assertEqual(MahlerCalculus.from_mahler(coeffs, x), f[x])

    def test_roundtrip_quadratic_function(self):
        f = [x * x for x in range(8)]
        coeffs = MahlerCalculus.to_mahler(f, 7)
        for x in range(8):
            self.assertEqual(MahlerCalculus.from_mahler(coeffs, x), f[x])

    def test_constant_function_has_only_a0(self):
        f = [42, 42, 42, 42]
        coeffs = MahlerCalculus.to_mahler(f, 3)
        self.assertEqual(coeffs, [42, 0, 0, 0])


class TestDiracOperator(unittest.TestCase):
    def test_D_e0_is_empty(self):
        self.assertEqual(MahlerCalculus.dirac_operator([1]), [])

    def test_D_e1_is_minus_e0(self):
        self.assertEqual(MahlerCalculus.dirac_operator([0, 1]), [-1])

    def test_D_e2_is_minus_e1(self):
        self.assertEqual(MahlerCalculus.dirac_operator([0, 0, 1]), [0, -1])

    def test_D_empty(self):
        self.assertEqual(MahlerCalculus.dirac_operator([]), [])


class TestVolterraOperator(unittest.TestCase):
    def test_T_empty(self):
        self.assertEqual(MahlerCalculus.volterra_operator([]), [])

    def test_T_raises_on_nonzero_a0(self):
        with self.assertRaises(ValueError):
            MahlerCalculus.volterra_operator([1])
        with self.assertRaises(ValueError):
            MahlerCalculus.volterra_operator([5, 3, 2])

    def test_T_e1_is_minus_e2(self):
        self.assertEqual(MahlerCalculus.volterra_operator([0, 1]), [0, 0, -1])

    def test_T_e2_is_minus_e3(self):
        self.assertEqual(MahlerCalculus.volterra_operator([0, 0, 1]), [0, 0, 0, -1])


class TestDtTdAsymmetry(unittest.TestCase):
    def test_D_e0_is_zero(self):
        self.assertEqual(MahlerCalculus.dirac_operator([1]), [])

    def test_D_e1_exits_ker_epsilon(self):
        De1 = MahlerCalculus.dirac_operator([0, 1])
        self.assertEqual(De1, [-1])
        self.assertNotEqual(De1[0], 0)

    def test_T_raises_outside_ker_epsilon(self):
        with self.assertRaises(ValueError):
            MahlerCalculus.volterra_operator([-1])
        with self.assertRaises(ValueError):
            MahlerCalculus.volterra_operator([1])

    def test_TD_e0_is_zero_not_e0(self):
        De0 = MahlerCalculus.dirac_operator([1])
        TDe0 = MahlerCalculus.volterra_operator(De0)
        self.assertEqual(TDe0, [])
        self.assertNotEqual(TDe0, [1])

    def test_TD_e1_raises(self):
        De1 = MahlerCalculus.dirac_operator([0, 1])
        with self.assertRaises(ValueError):
            MahlerCalculus.volterra_operator(De1)

    def test_DT_e1_is_e1(self):
        e1 = [0, 1]
        Te1 = MahlerCalculus.volterra_operator(e1)
        self.assertEqual(Te1, [0, 0, -1])
        DTe1 = MahlerCalculus.dirac_operator(Te1)
        self.assertEqual(DTe1, e1)

    def test_DT_identity_on_all_ker_epsilon(self):
        for n in range(1, 9):
            e_n = [0] * (n + 1)
            e_n[n] = 1
            Te_n = MahlerCalculus.volterra_operator(e_n)
            DTe_n = MahlerCalculus.dirac_operator(Te_n)
            self.assertEqual(DTe_n, e_n, f"D∘T(e_{n}) = {DTe_n}, expected {e_n}")

    def test_DT_identity_on_general_ker_epsilon_element(self):
        f = [0, 3, 2, -1, 7, 0, -5]
        Tf = MahlerCalculus.volterra_operator(f)
        DTf = MahlerCalculus.dirac_operator(Tf)
        self.assertEqual(DTf, f)

    def test_TD_identity_on_span_e2_and_above(self):
        for n in range(2, 9):
            e_n = [0] * (n + 1)
            e_n[n] = 1
            De_n = MahlerCalculus.dirac_operator(e_n)
            TDe_n = MahlerCalculus.volterra_operator(De_n)
            self.assertEqual(TDe_n, e_n, f"T∘D(e_{n}) = {TDe_n}, expected {e_n}")

    def test_T_is_right_inverse_not_left_at_e1(self):
        e1 = [0, 1]
        DTe1 = MahlerCalculus.dirac_operator(MahlerCalculus.volterra_operator(e1))
        self.assertEqual(DTe1, e1)
        with self.assertRaises(ValueError):
            MahlerCalculus.volterra_operator(MahlerCalculus.dirac_operator(e1))


class TestTruncate(unittest.TestCase):
    def test_truncate_mod_2k(self):
        coeffs = [1, 2, 3, 257]
        truncated = MahlerCalculus.truncate(coeffs, 8)
        self.assertEqual(truncated, [1, 2, 3, 1])

    def test_truncate_reduces_length(self):
        coeffs = list(range(10))
        truncated = MahlerCalculus.truncate(coeffs, 4)
        self.assertEqual(len(truncated), 4)


class TestCliffDensityTheory(unittest.TestCase):
    def test_cliff_distribution_exact(self):
        k = 14
        total = 1 << (k - 2)
        counts = {}
        for L_val in range(total):
            L_signed = L_val if L_val <= total // 2 else L_val - total
            lp4 = L_signed + 4
            v_lp4 = valuation(abs(lp4))
            c = max(0, v_lp4 - 2) if v_lp4 is not None else k
            counts[c] = counts.get(c, 0) + 1

        theory = {0: 7 / 8}
        for j in range(1, 10):
            theory[j] = 2 ** (-(j + 3))

        for c, expected_prob in theory.items():
            actual_prob = counts.get(c, 0) / total
            self.assertAlmostEqual(
                actual_prob,
                expected_prob,
                places=4,
                msg=f"c={c}: prob={actual_prob:.6f}, theory={expected_prob:.6f}",
            )

    def test_expected_cliff_cost_is_quarter(self):
        k = 14
        total = 1 << (k - 2)
        total_c = sum(
            max(0, valuation(abs((L - total if total // 2 < L else L) + 4)) - 2)
            if (L - total if total // 2 < L else L) + 4 != 0
            else k
            for L in range(total)
        )
        expected_c = total_c / total
        self.assertAlmostEqual(
            expected_c, 0.25, delta=0.01, msg=f"E[c] = {expected_c:.4f}, expected 0.25"
        )

    def test_high_cliff_rare(self):
        prob_high_cliff = sum(2 ** (-(j + 3)) for j in range(4, 100))
        self.assertAlmostEqual(prob_high_cliff, 2 ** (-6), places=6)


class TestDivisorOptimality(unittest.TestCase):
    def _step(self, e: int, target: int, k: int, d: int):
        mod = 1 << k
        mod_e = 1 << (k - 2)
        L = two_adic_log5(k) >> 2
        g_e = pow(5, e, mod)
        delta = (g_e - target) % mod
        if valuation(delta) < d:
            return None
        denom = (g_e * L) % mod_e
        return (e - (delta >> d) * pow(denom, -1, mod_e)) % mod_e

    def test_all_divisors_all_j(self):
        k = 24
        mod_e = 1 << (k - 2)
        target_e = 777 % mod_e
        target = pow(5, target_e, 1 << k)
        for d in range(6):
            for j in range(2, 7):
                e = (target_e + (1 << j)) % mod_e
                e_new = self._step(e, target, k, d)
                if e_new is None:
                    continue
                diff = (e_new - target_e) % mod_e
                v_new = valuation(diff) if diff != 0 else k
                if d < 2:
                    pred = j
                elif d == 2:
                    pred = min(2 * j + 1, k - 1)
                else:
                    pred = max(0, j - (d - 2))
                self.assertEqual(v_new, pred, f"d={d}, j={j}: v_new={v_new}, predicted={pred}")

    def test_d2_is_strictly_best(self):
        k = 20
        mod_e = 1 << (k - 2)
        target_e = 333 % mod_e
        target = pow(5, target_e, 1 << k)
        for j in range(2, 7):
            e = (target_e + (1 << j)) % mod_e
            v_d2 = valuation((self._step(e, target, k, 2) - target_e) % mod_e)
            for d in [0, 1, 3, 4]:
                e_new = self._step(e, target, k, d)
                if e_new is None:
                    continue
                diff = (e_new - target_e) % mod_e
                v_d = valuation(diff) if diff != 0 else k
                self.assertGreater(
                    v_d2, v_d, f"d=2 (v={v_d2}) should beat d={d} (v={v_d}) at j={j}"
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
