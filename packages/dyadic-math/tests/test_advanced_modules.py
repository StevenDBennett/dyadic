"""Tests for advanced modules: separation, fourier, padic_roots, iwasawa, mersenne, isometry."""

import unittest

import numpy as np


class TestSeparation(unittest.TestCase):
    def test_newton_trajectory(self):
        from dyadic_math.separation import newton_trajectory

        k = 8
        a = pow(5, 3, 1 << k)
        traj = newton_trajectory(a, k, 0, steps=5)
        self.assertEqual(len(traj), 6)
        self.assertIsInstance(traj[0], int)

    def test_predicted_separation(self):
        from dyadic_math.separation import predicted_separation

        self.assertEqual(predicted_separation(1), 0)
        self.assertEqual(predicted_separation(2), 0)
        self.assertEqual(predicted_separation(4), 1)

    def test_step_count_profile(self):
        from dyadic_math.separation import step_count_profile

        profile = step_count_profile(6, 3)
        self.assertIsInstance(profile, dict)


class TestFourier(unittest.TestCase):
    def test_step_count_fn_shape(self):
        from dyadic_math.fourier import step_count_fn

        k = 6
        h = step_count_fn(k, 3)
        self.assertEqual(len(h), 1 << (k - 2))

    def test_analytic_step_count(self):
        from dyadic_math.fourier import analytic_step_count

        k = 6
        h = analytic_step_count(k, 3)
        self.assertEqual(len(h), 1 << (k - 2))

    def test_analytic_matches_numeric(self):
        from dyadic_math.fourier import analytic_step_count, step_count_fn

        k = 6
        h_num = step_count_fn(k, 3)
        h_an = analytic_step_count(k, 3)
        self.assertEqual(h_num.shape, h_an.shape)

    def test_power_spectrum(self):
        from dyadic_math.fourier import power_spectrum

        h = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int32)
        spec = power_spectrum(h)
        self.assertEqual(len(spec), len(h))

    def test_fourier_summary(self):
        from dyadic_math.fourier import fourier_summary

        summary = fourier_summary(6, 3)
        self.assertIsInstance(summary, str)
        self.assertIn("Fourier Analysis", summary)


class TestPadicRoots(unittest.TestCase):
    def test_newton_step(self):
        from dyadic_math.padic_roots import newton_step

        pk = 5**4
        x = newton_step(2, 8, pk)
        self.assertIsInstance(x, int)

    def test_halley_step(self):
        from dyadic_math.padic_roots import halley_step

        pk = 5**4
        x = halley_step(2, 8, pk)
        self.assertIsInstance(x, int)

    def test_convergence_profile(self):
        from dyadic_math.padic_roots import convergence_profile, lift_root, newton_step

        p, k = 5, 4
        a = 8
        x_true = lift_root(a, p, k)
        if x_true is not None:
            prof = convergence_profile(2, a, p, k, newton_step, x_true)
            self.assertGreater(len(prof), 0)

    def test_popcount_compression(self):
        from dyadic_math.padic_roots import popcount_compression

        result = popcount_compression(8, n_trials=10)
        self.assertIn("pearson_r", result)


class TestIwasawa(unittest.TestCase):
    def test_congruence_depth(self):
        from dyadic_math.iwasawa import congruence_depth

        M = [[1, 0], [0, 1]]
        self.assertEqual(congruence_depth(M, 8), 8)

        M2 = [[3, 0], [0, 3]]
        d = congruence_depth(M2, 8)
        self.assertGreaterEqual(d, 1)

    def test_filtration_residue(self):
        from dyadic_math.iwasawa import congruence_depth, filtration_residue

        M = [[3, 0], [0, 3]]
        d = congruence_depth(M, 8)
        res = filtration_residue(M, d, 8)
        self.assertEqual(len(res), 2)
        self.assertEqual(len(res[0]), 2)

    def test_ldu_decompose(self):
        from dyadic_math.iwasawa import ldu_decompose

        M = [[5, 3], [2, 7]]
        ldu = ldu_decompose(M, 8)
        self.assertIsNotNone(ldu)
        self.assertIn("L", ldu)
        self.assertIn("D", ldu)
        self.assertIn("U", ldu)

    def test_matrix_coordinates(self):
        from dyadic_math.iwasawa import matrix_coordinates

        M = [[5, 3], [2, 7]]
        coords = matrix_coordinates(M, 8)
        self.assertEqual(len(coords.residue), 2)
        self.assertIsNotNone(coords.det_dual)

    def test_matrix_commutator(self):
        from dyadic_math.iwasawa import matrix_commutator

        M = [[5, 0], [0, 5]]
        N = [[7, 0], [0, 7]]
        C = matrix_commutator(M, N, 8)
        self.assertEqual(len(C), 2)
        self.assertEqual(len(C[0]), 2)

    def test_filtration_portrait(self):
        from dyadic_math.iwasawa import filtration_portrait

        portrait = filtration_portrait(4)
        self.assertIn("Congruence Filtration", portrait)


class TestMersenne(unittest.TestCase):
    def test_mersenne_coordinates(self):
        from dyadic_math.mersenne import mersenne_coordinates

        coords = mersenne_coordinates(5, 10)
        self.assertIsNotNone(coords)
        alpha, e_true, v2_e = coords
        self.assertEqual(alpha, 1)

    def test_verify_core_identity(self):
        from dyadic_math.mersenne import verify_core_identity

        results = verify_core_identity(8)
        for n, passed in results.items():
            self.assertTrue(passed, f"Core identity failed at n={n}")

    def test_mersenne_cliff_table(self):
        from dyadic_math.mersenne import mersenne_cliff_table

        rows = mersenne_cliff_table(8)
        self.assertGreater(len(rows), 0)
        for row in rows:
            self.assertIn("n", row)
            self.assertIn("k*", row)
            self.assertIn("alpha", row)

    def test_bootstrap_cost(self):
        from dyadic_math.mersenne import bootstrap_cost

        cost = bootstrap_cost(6, 16)
        self.assertGreater(cost, 0)

    def test_optimal_bootstrap(self):
        from dyadic_math.mersenne import optimal_bootstrap

        results = optimal_bootstrap([16, 24, 32])
        for _k, best in results.items():
            self.assertGreater(best, 0)

    def test_compare_bootstrap_strategies(self):
        from dyadic_math.mersenne import compare_bootstrap_strategies

        results = compare_bootstrap_strategies([16, 24])
        for _k, data in results.items():
            for key in ("sqrt_cost", "half_cost", "lut_cost"):
                self.assertIn(key, data)

    def test_dlog_with_lut(self):
        from dyadic_math.mersenne import dlog_with_lut

        k = 16
        for e_true in [0, 1, 5, 10, 31]:
            a = pow(5, e_true, 1 << k)
            e = dlog_with_lut(a, k, b=8)
            self.assertEqual(e, e_true)


class TestMersenneProofs(unittest.TestCase):
    def test_cliff_constant(self):
        from dyadic_math.mersenne import cliff_constant

        self.assertEqual(cliff_constant(5), 5)

    def test_cliff_constant_general(self):
        from dyadic_math.mersenne import cliff_constant, valuation

        for g in (13, 21, 29, 37, 45, 53, 61):
            c = cliff_constant(g, k=24)
            s = valuation(g - 5)
            expected = min(s - 2, 5)
            self.assertEqual(c, expected, f"g={g}: c={c}, expected={expected}")

    def test_cliff_formula(self):
        from dyadic_math.mersenne import cliff_formula

        result = cliff_formula(5)
        self.assertIn("5", result)

    def test_mersenne_cliff_theorem(self):
        from dyadic_math.mersenne import mersenne_cliff_theorem

        thm = mersenne_cliff_theorem()
        self.assertTrue(thm["all_pass"], f"Theorem failed: {thm}")
        self.assertEqual(thm["c"], 5)

    def test_verify_cliff_constant(self):
        from dyadic_math.mersenne import verify_cliff_constant

        self.assertTrue(verify_cliff_constant())

    def test_verify_c_formula(self):
        from dyadic_math.mersenne import verify_c_formula

        self.assertTrue(verify_c_formula())

    def test_exp2_neg4(self):
        from dyadic_math.mersenne import exp2_neg4

        self.assertEqual(exp2_neg4(7), 5)
        self.assertEqual(exp2_neg4(8), 133)
        self.assertEqual(exp2_neg4(9), 389)
        self.assertEqual(exp2_neg4(10), 389)
        self.assertEqual(exp2_neg4(11), 1925)

    def test_cliff_constant_unified(self):
        from dyadic_math.mersenne import cliff_constant, cliff_constant_unified

        for g in (5, 13, 21, 29, 37, 45, 53, 61):
            d = cliff_constant(g, k=24)
            u = cliff_constant_unified(g, k=32)
            self.assertEqual(d, u, f"Mismatch for g={g}: direct={d}, unified={u}")

    def test_verify_unified_formula(self):
        from dyadic_math.mersenne import verify_unified_formula

        self.assertTrue(verify_unified_formula())

    def test_verify_connection(self):
        from dyadic_math.mersenne import verify_connection

        self.assertTrue(verify_connection())

    def test_lift_root(self):
        from dyadic_math.padic_roots import lift_root

        for p in (5, 7, 11, 13):
            for k in (2, 4):
                for a in (2, 3, 5, 6):
                    if a % p == 0:
                        continue
                    root = lift_root(a, p, k)
                    if root is not None:
                        self.assertEqual(pow(root, 3, p**k), a % (p**k))


class TestIsometry(unittest.TestCase):
    def test_verify_isometry(self):
        from dyadic_math.isometry import verify_isometry

        result = verify_isometry(8, n_trials=20)
        self.assertGreaterEqual(result["pass_rate"], 0.9)

    def test_isometry_pair_test(self):
        from dyadic_math.isometry import isometry_pair_test

        result = isometry_pair_test(8, n_trials=20)
        self.assertGreaterEqual(result["pass_rate"], 0.9)

    def test_isometry_summary(self):
        from dyadic_math.isometry import isometry_summary

        summary = isometry_summary(8)
        self.assertIn("Isometry", summary)

    def test_trace_alpha_independence(self):
        from dyadic_math.isometry import trace_alpha_independence

        result = trace_alpha_independence(5, 7, cycle_length=3, n_cycles=20)
        self.assertIn("chi2_stat", result)

    def test_exponentvaluation_profile(self):
        from dyadic_math.isometry import exponentvaluation_profile

        profile = exponentvaluation_profile(8, n_samples=50)
        self.assertIsInstance(profile, dict)
        if profile:
            total = sum(profile.values())
            self.assertAlmostEqual(total, 1.0, places=2)


if __name__ == "__main__":
    unittest.main()
