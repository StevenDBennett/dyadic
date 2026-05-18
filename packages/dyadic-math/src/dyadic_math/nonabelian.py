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

from dyadic_core import bitmask, valuation, modinv_newton, two_adic_dlog, mat_mul, mat_det
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

    def holonomy(self, mats: list[list[list[int]]]) -> list[list[int]]:
        """Product of all matrices in the cycle modulo mod_full."""
        H = [[1, 0], [0, 1]]
        for M in mats:
            H = mat_mul(H, M, self.mod_full)
        return H

    def invariants(self, mats: list[list[list[int]]]) -> dict[str, int]:
        """
        Extract invariants: det mod 2^k, det dual view, trace mod p,
        and CRT merge.
        """
        H = self.holonomy(mats)
        det_mod2k = mat_det(H, self.mod2)
        trace_modp = (H[0][0] + H[1][1]) % self.p

        # Dual-view decomposition of determinant
        det_odd = det_mod2k if (det_mod2k & 1) else det_mod2k + 1
        dlog_result = two_adic_dlog(det_odd, self.k)
        if dlog_result is not None:
            alpha_det, e_det = dlog_result
        else:
            alpha_det, e_det = 0, 0

        # CRT merge
        inv_2k = self.crt.crt_reconstruct(det_mod2k, trace_modp)
        inv_p = self.crt.crt_reconstruct(trace_modp, det_mod2k)

        return {
            "det_mod2k": det_mod2k,
            "alpha_det": alpha_det,
            "e_det": e_det,
            "trace_modp": trace_modp,
            "crt_2k_view": inv_2k,
            "crt_p_view": inv_p,
        }

    def convergence_ratio_full(self, mats: list[list[list[int]]]) -> float:
        """Convergence ratio of the holonomy determinant 2-adic component."""
        from .basin import BasinExplorer

        H = self.holonomy(mats)
        det = mat_det(H, self.mod2)
        if det & 1 == 0:
            return 0.0
        try:
            explorer = BasinExplorer(self.k, 5, det)
            portrait = explorer.portrait()
            n_converged = len(portrait['converged'])
            total = n_converged + len(portrait['cycle'])
            return n_converged / total if total > 0 else 0.0
        except Exception:
            return 0.0


def phase_alignment_experiment(
    k: int, p: int, N_cycle: int = 4, n_cycles: int = 30,
) -> dict[str, float]:
    """
    Test whether a single-bit perturbation flips the α-sector of the
    holonomy determinant.  This is the replacement metric for the
    degenerate eps_crit in ramp_break_strength.

    Returns dict with 'alignment' (fraction of trials where a flip
    occurred) and 'n_trials'.
    """
    nc = NonAbelianCRTDual(k, p)
    flips = 0

    for _ in range(n_cycles):
        mats = [
            [[np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)],
             [np.random.randint(0, nc.mod_full), np.random.randint(0, nc.mod_full)]]
            for _ in range(N_cycle)
        ]

        inv = nc.invariants(mats)
        alpha0 = inv["alpha_det"]

        eps_mat = _perturbation_matrix(inv["det_mod2k"] + 1, nc.mod_full)
        perturbed = mats[:]
        perturbed[0] = mat_mul(perturbed[0], eps_mat, nc.mod_full)
        inv_p = nc.invariants(perturbed)

        if inv_p["alpha_det"] != alpha0:
            flips += 1

    return {
        "alignment": flips / n_cycles if n_cycles > 0 else 0.0,
        "n_trials": n_cycles,
    }
