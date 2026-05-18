"""Tests for dyadic_math.basin module."""

import contextlib
import io
import unittest

import numpy as np
from dyadic_math.basin import BasinExplorer, GhostHunt, LayerGhostDiagnosticV2, precision_sweep


class TestBasinExplorer(unittest.TestCase):
    def test_newton_step_converges_to_root(self):
        k = 8
        target_e = 3
        a = pow(5, target_e, 1 << k)
        explorer = BasinExplorer(k, 5, a)
        fate, val, _ = explorer.classify(target_e)
        self.assertEqual(fate, "converged")

    def test_positive_sector_all_converge(self):
        k = 8
        for e in range(8):
            a = pow(5, e, 1 << k)
            explorer = BasinExplorer(k, 5, a)
            port = explorer.portrait()
            self.assertGreater(len(port["converged"]), 0)

    def test_ghost_detection(self):
        k = 8
        a = pow(5, 3, 1 << k)
        explorer = BasinExplorer(k, 5, a)
        port = explorer.portrait()
        self.assertIsInstance(port, dict)
        self.assertIn("converged", port)
        self.assertIn("cycle", port)
        self.assertIn("diverged", port)

    def test_portrait_matrix_length(self):
        k = 6
        a = pow(5, 3, 1 << k)
        explorer = BasinExplorer(k, 5, a)
        fv = explorer.fate_vector()
        self.assertEqual(len(fv), explorer.N)

    def test_full_portrait_structure(self):
        k = 6
        a = pow(5, 3, 1 << k)
        explorer = BasinExplorer(k, 5, a)
        port = explorer.full_portrait()
        self.assertIn("converged", port)
        self.assertIn("cycles", port)
        self.assertIn("diverged", port)
        self.assertIsInstance(port["cycles"], dict)
        total = len(port["converged"]) + len(port["diverged"])
        for _period, entries in port["cycles"].items():
            total += len(entries)
        self.assertEqual(total, explorer.N)

    def test_full_portrait_cycles_have_periods(self):
        k = 6
        a = pow(5, 3, 1 << k)
        explorer = BasinExplorer(k, 5, a)
        port = explorer.full_portrait()
        for period, entries in port["cycles"].items():
            self.assertIsInstance(period, int)
            self.assertGreater(period, 0)
            for _seed, cycle_path in entries:
                self.assertEqual(len(cycle_path), period)

    def test_portrait_matrix_encoding(self):
        k = 6
        a = pow(5, 3, 1 << k)
        explorer = BasinExplorer(k, 5, a)
        pm = explorer.portrait_matrix()
        self.assertEqual(len(pm), explorer.N)
        for val in pm:
            self.assertLessEqual(val, 0)

    def test_invalid_g_accepts_any_odd(self):
        k = 8
        explorer = BasinExplorer(k, 3, 1)
        self.assertIsInstance(explorer, BasinExplorer)

    def test_invalid_k_raises(self):
        with self.assertRaises(ValueError):
            BasinExplorer(2, 5, 1)


class TestPrecisionSweep(unittest.TestCase):
    def test_sweep_returns_results(self):
        results = precision_sweep(4, 8, 5, 3)
        self.assertGreater(len(results), 0)
        for k, frac in results:
            self.assertIsInstance(k, int)
            self.assertIsInstance(frac, float)


class TestLayerGhostDiagnosticV2(unittest.TestCase):
    def test_diagnostic_matrix_shape(self):
        W = np.array([[1, 3], [5, 7]], dtype=np.int64)
        diag = LayerGhostDiagnosticV2(k=8)
        fate, conv, ghost, mean_e, v2_e = diag.diagnostic_matrix(W)
        self.assertEqual(fate.shape, W.shape)
        self.assertIsInstance(conv, float)
        self.assertIsInstance(ghost, float)

    def test_all_even_returns_negative_one(self):
        W = np.array([2, 4, 6], dtype=np.int64)
        diag = LayerGhostDiagnosticV2(k=8)
        fate, conv, ghost, _, _ = diag.diagnostic_matrix(W)
        self.assertTrue(np.all(fate == -1))
        self.assertEqual(conv, 0.0)
        self.assertEqual(ghost, 0.0)

    def test_powers_of_five_are_stable(self):
        W = np.array([pow(5, e, 256) for e in range(4)], dtype=np.int64)
        diag = LayerGhostDiagnosticV2(k=8)
        fate, conv, ghost, _, _ = diag.diagnostic_matrix(W)
        self.assertGreater(conv, 0.5)
        self.assertLess(ghost, 0.5)


class TestGhostHunt(unittest.TestCase):
    def test_precision_threshold_sweep(self):
        hunter = GhostHunt()
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            hunter.precision_threshold_sweep(4, 8, 3)
        output = f.getvalue()
        self.assertIn("Precision Sweep", output)

    def test_quantization_cliff(self):
        W = np.array([1, 3, 5, 7], dtype=np.int64)
        hunter = GhostHunt()
        results = hunter.quantization_cliff(W, k_min=4, k_max=10)
        self.assertIsInstance(results, dict)
        for k in range(4, 11):
            self.assertIn(k, results)
            self.assertGreaterEqual(results[k], 0.0)
            self.assertLessEqual(results[k], 1.0)


if __name__ == "__main__":
    unittest.main()
