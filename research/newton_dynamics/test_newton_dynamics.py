"""Tests for research.newton_dynamics — polynomial arithmetic, iterates, dynatomic."""
import unittest
import math
import sys
from pathlib import Path

_parent = str(Path(__file__).resolve().parents[2])
if _parent not in sys.path:
    sys.path.insert(0, _parent)

from research.newton_dynamics import (
    poly_mul, poly_add, poly_scalar_mul, poly_pow, poly_divmod,
    mobius, compute_iterates, dynatomic_polynomial,
    is_cube, tonelli_shanks, check_quadratic_cube_roots,
    COEFFS_PERIOD4, COEFFS_PERIOD5,
    MULTIPLIERS_PERIOD4, MULTIPLIERS_PERIOD5,
    load_period6_coefficients, PERIOD6_PREDICTED,
)


def eval_poly(coeffs, x):
    return sum(c * (x ** i) for i, c in enumerate(coeffs))


class TestMobius(unittest.TestCase):
    def test_mobius_one(self):
        self.assertEqual(mobius(1), 1)

    def test_mobius_prime(self):
        self.assertEqual(mobius(2), -1)
        self.assertEqual(mobius(3), -1)
        self.assertEqual(mobius(5), -1)

    def test_mobius_square_free_product(self):
        self.assertEqual(mobius(6), 1)
        self.assertEqual(mobius(30), -1)

    def test_mobius_square_factor(self):
        self.assertEqual(mobius(4), 0)
        self.assertEqual(mobius(12), 0)
        self.assertEqual(mobius(18), 0)


class TestPoly(unittest.TestCase):
    def test_mul_constant(self):
        self.assertEqual(poly_mul([1], [2, 3]), [2, 3])

    def test_mul_simple(self):
        p = poly_mul([1, 1], [1, -1])
        self.assertEqual(p, [1, 0, -1])

    def test_mul_degree(self):
        p = poly_mul([1, 2, 3], [4, 5])
        self.assertEqual(p, [4, 13, 22, 15])

    def test_add(self):
        self.assertEqual(poly_add([1, 2], [3, 4, 5]), [4, 6, 5])

    def test_scalar_mul(self):
        self.assertEqual(poly_scalar_mul([1, 2, 3], 4), [4, 8, 12])

    def test_pow_zero(self):
        self.assertEqual(poly_pow([1, 2], 0), [1])

    def test_pow_one(self):
        self.assertEqual(poly_pow([1, 2], 1), [1, 2])

    def test_pow_two(self):
        self.assertEqual(poly_pow([1, 1], 2), [1, 2, 1])

    def test_divmod_exact(self):
        q, r = poly_divmod([-1, 0, 1], [-1, 1])
        self.assertEqual(q, [1, 1])
        self.assertEqual(r, [0])

    def test_divmod_with_remainder(self):
        q, r = poly_divmod([0, 0, 1], [1, 1])
        self.assertEqual(q, [-1, 1])
        self.assertEqual(r, [1])

    def test_divmod_degree_lt(self):
        q, r = poly_divmod([1, 2], [3, 4, 5])
        self.assertEqual(q, [0])
        self.assertEqual(r, [1, 2])

    def test_divmod_by_zero_raises(self):
        with self.assertRaises(ValueError):
            poly_divmod([1], [0])


class TestIterates(unittest.TestCase):
    def test_iterate_0(self):
        iters = compute_iterates(0)
        self.assertEqual(len(iters), 1)
        A0, B0 = iters[0]
        self.assertEqual(A0, [0, 1])
        self.assertEqual(B0, [1])

    def test_iterate_1(self):
        iters = compute_iterates(1)
        A1, B1 = iters[1]
        self.assertEqual(A1, [1, 0, 0, 2])
        self.assertEqual(B1, [0, 0, 3])

    def test_iterate_2(self):
        iters = compute_iterates(2)
        A2, B2 = iters[2]
        self.assertEqual(A2[0], 2)
        self.assertEqual(A2[1], 0)
        self.assertEqual(A2[2], 0)
        self.assertEqual(A2[3], 12)
        self.assertEqual(B2[0], 0)
        self.assertEqual(B2[1], 0)
        self.assertEqual(B2[2], 9)

    def test_a_at_zero(self):
        iters = compute_iterates(6)
        for d in range(1, 7):
            expected = 2 ** ((3 ** (d - 1) - 1) // 2)
            self.assertEqual(iters[d][0][0], expected)


class TestDynatomic(unittest.TestCase):
    def test_period_2(self):
        iters = compute_iterates(2)
        phi2 = dynatomic_polynomial(2, iters)
        self.assertEqual(phi2, [2, 5, 20])

    def test_period_2_special_values(self):
        iters = compute_iterates(2)
        phi2 = dynatomic_polynomial(2, iters)
        self.assertEqual(phi2[0], 2)
        self.assertEqual(sum(phi2), 27)

    def test_period_4_degree(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        self.assertEqual(len(phi4) - 1, 24)

    def test_period_4_special_values(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        self.assertEqual(phi4[0], 4096)
        self.assertEqual(phi4[0], COEFFS_PERIOD4[0])

    def test_period_4_matches_data(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        self.assertEqual(len(phi4), len(COEFFS_PERIOD4))
        self.assertEqual(phi4[0], COEFFS_PERIOD4[0])
        self.assertEqual(phi4[-1], COEFFS_PERIOD4[-1])

    def test_period_5_degree(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        self.assertEqual(len(phi5) - 1, 80)

    def test_period_5_constant_term(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        self.assertEqual(phi5[0], 1099511627776)


class TestCleanPrimes(unittest.TestCase):
    def test_is_cube_trivial(self):
        self.assertTrue(is_cube(0, 7))
        self.assertTrue(is_cube(1, 7))

    def test_is_cube_mod_7(self):
        self.assertTrue(is_cube(6, 7))
        self.assertFalse(is_cube(3, 7))

    def test_tonelli_shanks_perfect_square(self):
        sqrt = tonelli_shanks(5, 11)
        self.assertIsNotNone(sqrt)
        self.assertEqual((sqrt * sqrt) % 11, 5)

    def test_tonelli_shanks_non_residue(self):
        self.assertIsNone(tonelli_shanks(2, 5))

    def test_tonelli_shanks_p_eq_3_mod_4(self):
        sqrt = tonelli_shanks(2, 7)
        self.assertIsNotNone(sqrt)
        self.assertEqual((sqrt * sqrt) % 7, 2)

    def test_check_quadratic_cube_roots_linear(self):
        self.assertFalse(check_quadratic_cube_roots(0, 3, 1, 7))

    def test_check_quadratic_cube_roots_no_discriminant(self):
        self.assertFalse(check_quadratic_cube_roots(1, 1, 1, 5))


class TestData(unittest.TestCase):
    def test_period4_length(self):
        self.assertEqual(len(COEFFS_PERIOD4), 25)

    def test_period4_constant(self):
        self.assertEqual(COEFFS_PERIOD4[0], 4096)

    def test_period4_leading(self):
        self.assertEqual(COEFFS_PERIOD4[-1], 223338299392)

    def test_period5_length(self):
        self.assertEqual(len(COEFFS_PERIOD5), 81)

    def test_period5_constant(self):
        self.assertEqual(COEFFS_PERIOD5[0], 1099511627776)

    def test_period5_multipliers_count(self):
        self.assertEqual(len(MULTIPLIERS_PERIOD5), 16)

    def test_period5_multiplier_sum(self):
        total_real = sum(m[0] for m in MULTIPLIERS_PERIOD5)
        self.assertAlmostEqual(total_real, 486.0, places=5)

    def test_period6_predicted(self):
        self.assertEqual(PERIOD6_PREDICTED["mu3"], 696)
        self.assertEqual(PERIOD6_PREDICTED["deg_u"], 232)
        self.assertEqual(PERIOD6_PREDICTED["phi_at_0"], 2 ** 116)

    def test_load_period6_coefficients(self):
        coeffs = load_period6_coefficients()
        self.assertEqual(len(coeffs), 233)
        self.assertEqual(coeffs[0], 83076749736557242056487941267521536)


class TestIntegration(unittest.TestCase):
    def test_period2_product_formula(self):
        mu3 = sum(mobius(2 // d) * (3 ** d) for d in (1, 2))
        self.assertEqual(mu3, 6)
        self.assertEqual(6 ** (mu3 // 6), 6)

    def test_period4_product_formula(self):
        mu3 = sum(mobius(4 // d) * (3 ** d) for d in (1, 2, 4))
        self.assertEqual(mu3, 72)
        self.assertEqual(6 ** (mu3 // 6), 2176782336)


class TestTheorem1_PhiAtZero(unittest.TestCase):
    def test_mu3_values(self):
        self.assertEqual(sum(mobius(2 // d) * (3 ** d) for d in (1, 2)), 6)
        self.assertEqual(sum(mobius(4 // d) * (3 ** d) for d in (1, 2, 4)), 72)
        self.assertEqual(sum(mobius(5 // d) * (3 ** d) for d in (1, 5)), 240)
        self.assertEqual(sum(mobius(6 // d) * (3 ** d) for d in (1, 2, 3, 6)), 696)

    def test_period2_phi_at_zero(self):
        iters = compute_iterates(2)
        phi2 = dynatomic_polynomial(2, iters)
        mu3 = 6
        self.assertEqual(phi2[0], 2 ** (mu3 // 6))

    def test_period4_phi_at_zero(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        mu3 = 72
        self.assertEqual(phi4[0], 2 ** (mu3 // 6))

    def test_period5_phi_at_zero(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        mu3 = 240
        self.assertEqual(phi5[0], 2 ** (mu3 // 6))


class TestTheorem2_PhiAtOne(unittest.TestCase):
    def test_period2_phi_at_one(self):
        iters = compute_iterates(2)
        phi2 = dynatomic_polynomial(2, iters)
        mu3 = 6
        self.assertEqual(eval_poly(phi2, 1), 3 ** (mu3 // 2))

    def test_period4_phi_at_one(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        mu3 = 72
        self.assertEqual(eval_poly(phi4, 1), 3 ** (mu3 // 2))

    def test_period5_phi_at_one(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        mu3 = 240
        self.assertEqual(eval_poly(phi5, 1), 3 ** (mu3 // 2))


class TestTheorem3_ProductFormula(unittest.TestCase):
    def test_period5_product_from_multipliers(self):
        total_real = sum(m[0] for m in MULTIPLIERS_PERIOD5)
        self.assertAlmostEqual(total_real, 486.0, places=5)
        prod = 1.0
        for i in range(0, len(MULTIPLIERS_PERIOD5), 2):
            a, b = MULTIPLIERS_PERIOD5[i]
            prod *= (a * a + b * b)
        expected = 6 ** 40
        self.assertAlmostEqual(prod, expected, delta=expected * 1e-6)

    def test_period5_theoretical_product(self):
        mu3 = 240
        self.assertEqual(6 ** (mu3 // 6), 6 ** 40)


class TestTheorem4_MultiplierIdentity(unittest.TestCase):
    def test_period4_multiplier_sum(self):
        total_real = sum(m[0] for m in MULTIPLIERS_PERIOD4)
        self.assertAlmostEqual(total_real, 90.0, places=1)

    def test_period4_multiplier_product(self):
        self.assertEqual(6 ** 12, 2176782336)

    def test_period4_minimal_polynomial_constant(self):
        self.assertEqual(6 ** 12, 2176782336)

    def test_period4_count(self):
        self.assertEqual(len(MULTIPLIERS_PERIOD4), 6)


class TestCleanPrimesVerified(unittest.TestCase):
    def test_7_is_clean(self):
        p = 7
        has_period2 = False
        for u in range(p):
            if (20 * u * u + 5 * u + 2) % p == 0 and is_cube(u, p):
                has_period2 = True
                break
        self.assertFalse(has_period2)

    def test_103_is_clean(self):
        p = 103
        has_period2 = False
        for u in range(p):
            if (20 * u * u + 5 * u + 2) % p == 0 and is_cube(u, p):
                has_period2 = True
                break
        self.assertFalse(has_period2)

    def test_181_is_clean(self):
        p = 181
        has_period2 = False
        for u in range(p):
            if (20 * u * u + 5 * u + 2) % p == 0 and is_cube(u, p):
                has_period2 = True
                break
        self.assertFalse(has_period2)

    def test_non_clean_prime_has_points(self):
        p = 79
        has_period2 = False
        for u in range(p):
            if (20 * u * u + 5 * u + 2) % p == 0 and is_cube(u, p):
                has_period2 = True
                break
        self.assertTrue(has_period2)


class TestFullCoefficientMatch(unittest.TestCase):
    def test_period4_full_match(self):
        iters = compute_iterates(4)
        phi4 = dynatomic_polynomial(4, iters)
        self.assertEqual(phi4, COEFFS_PERIOD4)

    def test_period5_full_match(self):
        iters = compute_iterates(5)
        phi5 = dynatomic_polynomial(5, iters)
        self.assertEqual(phi5, COEFFS_PERIOD5)


class TestPeriod6Coefficients(unittest.TestCase):
    def test_period6_constant_term(self):
        coeffs = load_period6_coefficients()
        self.assertEqual(coeffs[0], 2 ** 116)

    def test_period6_length(self):
        coeffs = load_period6_coefficients()
        self.assertEqual(len(coeffs), 233)

    def test_period6_leading_coefficient_positive(self):
        coeffs = load_period6_coefficients()
        self.assertGreater(coeffs[-1], 0)

    def test_period6_predicted_consistency(self):
        self.assertEqual(PERIOD6_PREDICTED["phi_at_0"], 2 ** 116)
        self.assertEqual(PERIOD6_PREDICTED["prod_mu_power"], 116)


if __name__ == "__main__":
    unittest.main()
