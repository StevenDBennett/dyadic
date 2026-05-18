"""
nonabelian.py
-------------
Non-Abelian CRT-dual system with GL(2) holonomy.

GL(2) holonomy over a cycle in Z/(2^k · p)Z with CRT-based
crossed invariants.  The holonomy H = Π M_i
lives in GL(2), and its determinant (mod 2^k) and trace (mod p)
provide crossed invariants that combine via CRT.

Note on ramp_break_strength
---------------------------
The original eps_crit metric in ramp_break_strength is degenerate:
for a shear matrix [[1+ε, 1], [0, 1]], det = 1+ε, and the smallest
ε giving an even determinant is always 1.  The genuinely interesting
signal is whether the perturbation flips the α-sector of the
determinant (phase alignment).  Use phase_alignment_experiment()
for this purpose.
"""

from __future__ import annotations

import numpy as np
from dyadic_core import mat_det, mat_mul, two_adic_dlog

from .crt import CRTDualProcessor


def _perturbation_matrix(d_target: int, mod_full: int) -> list[list[int]]:
    """Create a shear matrix [[d, 1], [0, 1]] mod mod_full."""
    return [[d_target % mod_full, 1], [0, 1]]


class NonAbelianCRTDual:
    """
    Non-Abelian (matrix-valued) CRT-dual system.

    Computes holonomy for a cycle of 2×2 matrices over
    Z/(2^k · p)Z, extracts joint invariants, and tests
    their correlation with 2-adic ghost stability.
    """

    def __init__(self, k: int, p: int) -> None:
        self.k = k
        self.p = p
        self.mod2 = 1 << k
        self.mod_full = self.mod2 * p
        self.crt = CRTDualProcessor(k, p)

    def holonomy(self, matrices: list[list[list[int]]]) -> list[list[int]]:
        """Product of all matrices in the cycle modulo mod_full."""
        result = [[1, 0], [0, 1]]
        for matrix in matrices:
            result = mat_mul(result, matrix, self.mod_full)
        return result

    def invariants(self, matrices: list[list[list[int]]]) -> dict[str, int | None]:
        """
        Extract invariants: det mod 2^k, det dual view, trace mod p,
        and CRT merge.

        For even determinants v2(det) >= 1, the 2-adic discrete log is
        undefined and ``alpha_det``/``e_det`` are set to ``None``.
        Check ``det_even`` to detect this case.
        """
        holonomy = self.holonomy(matrices)
        det_mod2k = mat_det(holonomy, self.mod2)
        det_even_flag = (det_mod2k & 1) == 0
        trace_modp = (holonomy[0][0] + holonomy[1][1]) % self.p

        # Dual-view decomposition (only defined for odd inputs)
        if not det_even_flag:
            dlog_result = two_adic_dlog(det_mod2k, self.k)
            if dlog_result is not None:
                alpha_det, e_det = dlog_result
            else:
                alpha_det, e_det = None, None
        else:
            alpha_det, e_det = None, None

        # CRT merge
        inv_2k = self.crt.crt_reconstruct(det_mod2k, trace_modp)
        inv_p = self.crt.crt_reconstruct(trace_modp, det_mod2k)

        return {
            "det_mod2k": det_mod2k,
            "det_even": det_even_flag,
            "alpha_det": alpha_det,
            "e_det": e_det,
            "trace_modp": trace_modp,
            "crt_2k_view": inv_2k,
            "crt_p_view": inv_p,
        }

    def convergence_ratio_full(self, matrices: list[list[list[int]]]) -> float:
        """Convergence ratio of the holonomy determinant 2-adic component."""
        from .basin import BasinExplorer

        holonomy = self.holonomy(matrices)
        det = mat_det(holonomy, self.mod2)
        if det & 1 == 0:
            return 0.0
        try:
            explorer = BasinExplorer(self.k, 5, det)
            portrait = explorer.portrait()
            n_converged = len(portrait["converged"])
            total = n_converged + len(portrait["cycle"])
            return n_converged / total if total > 0 else 0.0
        except (ValueError, ArithmeticError):
            return 0.0


def phase_alignment_experiment(
    k: int,
    p: int,
    cycle_length: int = 4,
    n_cycles: int = 30,
    seed: int | None = None,
) -> dict[str, float]:
    """
    Test whether a single-bit perturbation flips the α-sector of the
    holonomy determinant.  This is the replacement metric for the
    degenerate eps_crit in ramp_break_strength.

    Returns dict with 'alignment' (fraction of trials where a flip
    occurred) and 'n_trials'.
    """
    if seed is not None:
        np.random.seed(seed)
    nc = NonAbelianCRTDual(k, p)
    flips = 0
    valid = 0

    for _ in range(n_cycles):
        mats = [
            [
                [np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)],
                [np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)],
            ]
            for _ in range(cycle_length)
        ]

        inv = nc.invariants(mats)
        if inv["det_even"]:
            continue
        alpha0 = inv["alpha_det"]
        if alpha0 is None:
            continue

        det_m = inv["det_mod2k"]
        assert det_m is not None
        eps_mat = _perturbation_matrix(det_m + 1, nc.mod_full)
        perturbed = mats[:]
        perturbed[0] = mat_mul(perturbed[0], eps_mat, nc.mod_full)
        inv_p = nc.invariants(perturbed)
        if inv_p["det_even"]:
            continue
        alpha1 = inv_p["alpha_det"]
        if alpha1 is None:
            continue

        valid += 1
        if alpha1 != alpha0:
            flips += 1

    return {
        "alignment": flips / valid if valid > 0 else 0.0,
        "n_trials": valid,
    }
