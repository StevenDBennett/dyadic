"""Tests for dyadic_math.crt and dyadic_math.nonabelian."""
import unittest

from dyadic_math.crt import (
    CRTDualNumber, CRTDualProcessor, combined_stability,
    _primitive_root, _prime_dlog,
)
from dyadic_math.nonabelian import (
    NonAbelianCRTDual, phase_alignment_experiment,
)
from dyadic_core import mat_mul, mat_det


class TestCRTHelpers(unittest.TestCase):
    def test_primitive_root_known(self):
        roots = {3: 2, 5: 2, 7: 3, 11: 2, 13: 2, 17: 3}
        for p, expected in roots.items():
            self.assertEqual(_primitive_root(p), expected)

    def test_prime_dlog_known(self):
        self.assertEqual(_prime_dlog(1, 7, 3), 0)
        self.assertEqual(_prime_dlog(3, 7, 3), 1)


class TestCRTDualNumber(unittest.TestCase):
    def test_initialization(self):
        n = CRTDualNumber(42, k=6, p=7, g_p=3)
        self.assertIsNotNone(n.component_2)
        self.assertEqual(n.residue_p, 42 % 7)

    def test_verify(self):
        n = CRTDualNumber(42, k=6, p=7, g_p=3)
        self.assertTrue(n.verify())

    def test_zero_component(self):
        n = CRTDualNumber(0, k=6, p=7, g_p=3)
        self.assertTrue(n.component_2.is_zero)


class TestCRTDualProcessor(unittest.TestCase):
    def setUp(self):
        self.proc = CRTDualProcessor(k=6, p=7)

    def test_crt_reconstruct(self):
        r = self.proc.crt_reconstruct(5, 3)
        self.assertEqual(r % (1 << 6), 5)
        self.assertEqual(r % 7, 3)

    def test_mul(self):
        a = CRTDualNumber(3, k=6, p=7, g_p=self.proc.g_p)
        b = CRTDualNumber(5, k=6, p=7, g_p=self.proc.g_p)
        c = self.proc.mul(a, b)
        self.assertTrue(c.verify())
        self.assertEqual(c.component_2.value, (3 * 5) & ((1 << 6) - 1))
        self.assertEqual(c.residue_p, (3 * 5) % 7)

    def test_product(self):
        P = self.proc.product([3, 5, 7])
        self.assertTrue(P.verify())
        expected_2 = (3 * 5 * 7) & ((1 << 6) - 1)
        expected_p = (3 * 5 * 7) % 7
        self.assertEqual(P.component_2.value, expected_2)
        self.assertEqual(P.residue_p, expected_p)

    def test_cycle_product(self):
        nums = [CRTDualNumber(i, k=6, p=7, g_p=self.proc.g_p) for i in [3, 5, 7]]
        P = self.proc.cycle_product(nums)
        self.assertTrue(P.verify())

    def test_convergence_ratio(self):
        n = CRTDualNumber(5, k=6, p=7, g_p=self.proc.g_p)
        r = self.proc.convergence_ratio_2adic(n)
        self.assertGreaterEqual(r, 0.0)
        self.assertLessEqual(r, 1.0)


class TestCombinedStability(unittest.TestCase):
    def test_returns_dict(self):
        result = combined_stability(k=6, p=7, num_cycles=10)
        self.assertIn("pearson_r", result)
        self.assertIn("n_samples", result)
        self.assertIn("mean_v2", result)


class TestNonAbelianMatrixOps(unittest.TestCase):
    def testmat_mul(self):
        A = [[1, 2], [3, 4]]
        B = [[5, 6], [7, 8]]
        C = mat_mul(A, B, 100)
        self.assertEqual(C[0][0], (1*5 + 2*7) % 100)
        self.assertEqual(C[0][1], (1*6 + 2*8) % 100)

    def testmat_det(self):
        M = [[1, 2], [3, 4]]
        det = mat_det(M, 100)
        self.assertEqual(det, (1*4 - 2*3) % 100)


class TestNonAbelianCRTDual(unittest.TestCase):
    def setUp(self):
        self.nc = NonAbelianCRTDual(k=6, p=7)

    def test_holonomy_identity(self):
        mats = [[[1, 0], [0, 1]] for _ in range(3)]
        H = self.nc.holonomy(mats)
        self.assertEqual(H, [[1, 0], [0, 1]])

    def test_invariants_keys(self):
        mats = [[[2, 1], [1, 2]]]
        inv = self.nc.invariants(mats)
        for key in ("det_mod2k", "alpha_det", "trace_modp"):
            self.assertIn(key, inv)

    def test_convergence_ratio(self):
        mats = [[[2, 1], [1, 2]]]
        r = self.nc.convergence_ratio_full(mats)
        self.assertGreaterEqual(r, 0.0)


class TestPhaseAlignmentExperiment(unittest.TestCase):
    def test_phase_alignment_experiment(self):
        result = phase_alignment_experiment(k=5, p=7, n_cycles=10)
        self.assertIn("alignment", result)
        self.assertIn("n_trials", result)


if __name__ == "__main__":
    unittest.main()
