"""Tests for dyadic_math.weight_stability — WeightStabilityDiagnostics."""

import unittest

import numpy as np
from dyadic_math.weight_stability import WeightStabilityDiagnostics


class TestWeightStabilityDiagnostics(unittest.TestCase):
    def test_analyze_returns_dict(self):
        W = np.array([1, 3, 5, 7, 255], dtype=np.int64)
        st = WeightStabilityDiagnostics(k=8)
        stats = st.analyse(W)
        self.assertIsInstance(stats, dict)
        self.assertIn("alpha_fraction", stats)
        self.assertIn("mean_v2_e", stats)
        self.assertIn("cliff_risk", stats)

    def test_powers_of_five_have_no_cliff(self):
        W = np.array([pow(5, e, 256) for e in range(8)], dtype=np.int64)
        st = WeightStabilityDiagnostics(k=8)
        stats = st.analyse(W)
        self.assertLess(stats["cliff_risk"], 0.5)

    def test_even_weights_handled(self):
        W = np.array([2, 4, 6, 8, 10], dtype=np.int64)
        st = WeightStabilityDiagnostics(k=8)
        stats = st.analyse(W)
        self.assertIsInstance(stats, dict)

    def test_summary_keys(self):
        W = np.random.randint(0, 256, size=(10,), dtype=np.int64)
        st = WeightStabilityDiagnostics.from_precision_sweep(W, range(4, 10), k=8)
        st.compute()
        s = st.summary()
        for key in ("n_weights", "n_stable", "n_ghost", "ghost_fraction"):
            self.assertIn(key, s)

    def test_report_output(self):
        W = np.random.randint(0, 256, size=(10,), dtype=np.int64)
        st = WeightStabilityDiagnostics.from_precision_sweep(W, range(4, 10), k=8)
        st.compute()
        report = st.report()
        self.assertIsInstance(report, str)
        self.assertIn("WeightStabilityDiagnostics Report", report)

    def test_stable_weights_with_max_k(self):
        W = np.random.randint(0, 256, size=(20,), dtype=np.int64)
        st = WeightStabilityDiagnostics.from_precision_sweep(W, range(4, 10), k=8)
        st.compute()
        stable = st.stable_weights(max_k=6)
        ghost = st.ghost_weights(max_k=6)
        self.assertIsInstance(stable, list)
        self.assertIsInstance(ghost, list)
        self.assertEqual(len(stable) + len(ghost), len(st.cliffs))

    def test_stable_weights_no_arg(self):
        W = np.random.randint(0, 256, size=(10,), dtype=np.int64)
        st = WeightStabilityDiagnostics.from_precision_sweep(W, range(4, 10), k=8)
        st.compute()
        stable = st.stable_weights()
        ghost = st.ghost_weights()
        self.assertEqual(len(stable) + len(ghost), len(st.cliffs))

    def test_even_returns_none_cliff(self):
        st = WeightStabilityDiagnostics(k=16)
        self.assertIsNone(st.mersenne_cliff_score(4))

    def test_mersenne_number_has_expected_cliff(self):
        st = WeightStabilityDiagnostics(k=16)
        score = st.mersenne_cliff_score(31)
        self.assertIsNotNone(score)
        self.assertEqual(score, 5)

    def test_power_of_five_has_high_cliff(self):
        st = WeightStabilityDiagnostics(k=16)
        score = st.mersenne_cliff_score(5)
        self.assertIsNotNone(score)
        self.assertGreaterEqual(score, 2)

    def test_compare_to_random_returns_z_scores(self):
        st = WeightStabilityDiagnostics(k=8)
        W = np.array([1, 3, 5, 7, 255, 127, 63], dtype=np.int64)
        result = st.compare_to_random(W, n_samples=50)
        self.assertIn("z_alpha", result)
        self.assertIn("z_v2_e", result)
        self.assertIn("z_cliff_risk", result)
        for key in result:
            self.assertIsInstance(result[key], float)

    def test_cliff_histogram_returns_dict_of_ints(self):
        W = np.array([1, 3, 5, 7, 255, 127, 63, 31], dtype=np.int64)
        st = WeightStabilityDiagnostics.from_precision_sweep(W, range(4, 10), k=8)
        st.compute()
        hist = st.cliff_histogram()
        self.assertIsInstance(hist, dict)
        for k_val, count in hist.items():
            self.assertIsInstance(k_val, int)
            self.assertIsInstance(count, int)
            self.assertGreater(count, 0)


if __name__ == "__main__":
    unittest.main()
