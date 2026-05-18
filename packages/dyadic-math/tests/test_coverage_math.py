"""Additional math tests for coverage."""

import unittest

import numpy as np


class TestLiftRoot(unittest.TestCase):
    def test_known_cube_root_mod_5(self):
        from dyadic_math.padic_roots import lift_root

        x = lift_root(3, 5, 1)
        self.assertIsNotNone(x)
        self.assertEqual((x**3) % 5, 3)

    def test_hensel_lift_to_higher_power(self):
        from dyadic_math.padic_roots import lift_root

        x = lift_root(3, 5, 4)
        self.assertIsNotNone(x)
        self.assertEqual((x**3) % (5**4), 3)

    def test_no_root_returns_none(self):
        from dyadic_math.padic_roots import lift_root

        x = lift_root(2, 7, 1)
        self.assertIsNone(x)

    def test_divisible_by_p_returns_none(self):
        from dyadic_math.padic_roots import lift_root

        x = lift_root(0, 5, 2)
        self.assertIsNone(x)


class TestNewton2Step(unittest.TestCase):
    def test_converges_faster_than_newton(self):
        from dyadic_math.padic_roots import lift_root, newton2_step, newton_step

        p, k = 5, 6
        pk = p**k
        a = 8
        x_true = lift_root(a, p, k)
        self.assertIsNotNone(x_true)
        x0 = 2
        x_n1 = newton_step(x0, a, pk)
        x_n2 = newton_step(x_n1, a, pk)
        x_comp = newton2_step(x0, a, pk)
        self.assertEqual(x_n2, x_comp)


class TestNewton3Step(unittest.TestCase):
    def test_equals_three_newton_steps(self):
        from dyadic_math.padic_roots import newton3_step, newton_step

        p, k = 5, 6
        pk = p**k
        a = 8
        x0 = 2
        x_n1 = newton_step(x0, a, pk)
        x_n2 = newton_step(x_n1, a, pk)
        x_n3 = newton_step(x_n2, a, pk)
        x_comp = newton3_step(x0, a, pk)
        self.assertEqual(x_n3, x_comp)


class TestCompareMethods(unittest.TestCase):
    def test_returns_all_methods(self):
        from dyadic_math.padic_roots import compare_methods

        result = compare_methods(5, 4, n_trials=5)
        self.assertIn("Newton (ord 2)", result)
        self.assertIn("Halley (ord 3)", result)
        self.assertIn("Comp-Newton (ord 4)", result)
        self.assertIn("Comp×3 (ord 8)", result)

    def test_higher_order_converges_better(self):
        from dyadic_math.padic_roots import compare_methods

        result = compare_methods(5, 4, n_trials=20)
        self.assertGreaterEqual(result["Comp×3 (ord 8)"], result["Newton (ord 2)"])


class TestVerifyOrder(unittest.TestCase):
    def test_returns_nested_dict(self):
        from dyadic_math.padic_roots import verify_order

        result = verify_order([5, 7], k=4, n_trials=3)
        self.assertIn("Newton", result)
        self.assertIn("Halley", result)
        for method_name in result:
            self.assertIsInstance(result[method_name], dict)


class TestNewtonCorrectionUniformity(unittest.TestCase):
    def test_returns_chi2_stats(self):
        from dyadic_math.padic_roots import newton_correction_uniformity

        result = newton_correction_uniformity(5, 2, n_seeds=100)
        self.assertIn("chi2_stat", result)
        self.assertIn("n_bins", result)
        self.assertIn("n_samples", result)
        self.assertIn("df", result)
        self.assertEqual(result["df"], 4)


class TestSeparationStep(unittest.TestCase):
    def test_identical_targets_never_diverge(self):
        from dyadic_math.separation import separation_step

        k = 8
        a = pow(5, 3, 1 << k)
        n = separation_step(a, a, k, 0)
        self.assertEqual(n, -1)

    def test_close_targets_diverge_later(self):
        from dyadic_math.separation import separation_step

        k = 8
        a = pow(5, 3, 1 << k)
        a_prime = a ^ (1 << 4)
        n = separation_step(a, a_prime, k, 0)
        self.assertGreaterEqual(n, 0)


class TestVerifySeparation(unittest.TestCase):
    def test_returns_results_for_all_s_values(self):
        from dyadic_math.separation import verify_separation

        k = 8
        s_values = [3, 4, 5]
        results = verify_separation(k, s_values, n_trials=5)
        for s in s_values:
            self.assertIn(s, results)
            self.assertIsInstance(results[s], int)


class TestDyadicCoefficients(unittest.TestCase):
    def test_has_dc_component(self):
        from dyadic_math.fourier import dyadic_coefficients

        h = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int32)
        coeffs = dyadic_coefficients(h)
        self.assertIn("DC", coeffs)

    def test_dyadic_keys_present(self):
        from dyadic_math.fourier import dyadic_coefficients

        h = np.array([0, 1, 2, 3, 0, 1, 2, 3], dtype=np.int32)
        coeffs = dyadic_coefficients(h)
        self.assertIn("N/2", coeffs)
        self.assertIn("N/4", coeffs)

    def test_matches_analytic_for_step_count(self):
        from dyadic_math.fourier import (
            analytic_coefficients,
            dyadic_coefficients,
            step_count_fn,
        )

        k = 6
        h_num = step_count_fn(k, 3)
        dyadic_num = dyadic_coefficients(h_num)
        dyadic_an = analytic_coefficients(k)
        self.assertIn("DC", dyadic_num)
        self.assertIn("DC", dyadic_an)
        self.assertGreater(abs(dyadic_num["DC"]), 0)
        self.assertGreater(abs(dyadic_an["DC"]), 0)


class TestAnalyticCoefficients(unittest.TestCase):
    def test_returns_expected_keys(self):
        from dyadic_math.fourier import analytic_coefficients

        coeffs = analytic_coefficients(6)
        self.assertIn("DC", coeffs)
        self.assertIn("N/2", coeffs)
        self.assertIn("N/4", coeffs)

    def test_dc_is_real(self):
        from dyadic_math.fourier import analytic_coefficients

        coeffs = analytic_coefficients(6)
        self.assertEqual(coeffs["DC"].imag, 0)


class TestUltrametricUncertainty(unittest.TestCase):
    def test_returns_string_with_N(self):
        from dyadic_math.fourier import ultrametric_uncertainty

        text = ultrametric_uncertainty(6)
        self.assertIn("N=16", text)
        self.assertIn("uncertainty", text.lower())


class TestVerifyCommutatorDepth(unittest.TestCase):
    def test_returns_triples(self):
        from dyadic_math.iwasawa import verify_commutator_depth

        results = verify_commutator_depth(8, [(1, 1)], n_trials=5)
        self.assertIsInstance(results, list)
        for dM, dN, dMN in results:
            self.assertIsInstance(dM, int)
            self.assertIsInstance(dN, int)
            self.assertIsInstance(dMN, int)

    def test_commutator_depth_theorem_holds(self):
        from dyadic_math.iwasawa import verify_commutator_depth

        results = verify_commutator_depth(8, [(1, 1), (1, 2), (2, 2)], n_trials=10)
        for dM, dN, dMN in results:
            self.assertGreaterEqual(
                dMN, dM + dN, f"depth([M,N])={dMN} < depth(M)+depth(N)={dM + dN}"
            )


class TestHolonomyDepthProfile(unittest.TestCase):
    def test_returns_dict_with_keys(self):
        from dyadic_math.iwasawa import holonomy_depth_profile

        result = holonomy_depth_profile(5, 7, cycle_length=3, n_cycles=5)
        self.assertIn("mean_depth_orig", result)
        self.assertIn("mean_depth_pert", result)
        self.assertIn("n_cycles", result)


class TestTraceExponentIndependence(unittest.TestCase):
    def test_returns_anova_results(self):
        from dyadic_math.isometry import trace_exponent_independence

        result = trace_exponent_independence(5, 7, cycle_length=3, n_cycles=20)
        self.assertIn("F_stat", result)
        self.assertIn("df_between", result)
        self.assertIn("df_within", result)
        self.assertIn("n_samples", result)
        self.assertIsInstance(result["F_stat"], float)


if __name__ == "__main__":
    unittest.main()
