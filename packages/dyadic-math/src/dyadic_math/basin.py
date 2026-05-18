"""
basin.py
--------
Newton basin analysis for the discrete-logarithm map.

Classifies all seeds e ∈ Z/2^(k-2) by their Newton convergence fate
(converged, ghost cycle, or diverged) for a given target a mod 2^k.

Critical note on the α=1 fix
----------------------------
For targets with α = 1 (i.e. a ≡ 3 mod 4), the Newton iteration
solves 5^e ≡ -a (mod 2^k) not 5^e ≡ a.  Without this fix the
iteration converges to a ghost fixed point at e* = dlog(a+2, k).
This fix was identified in the original 2-Adic-Newton-Dynamics
codebase and is essential for correctness.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
from dyadic_core import bitmask, modinv_newton, two_adic_dlog, two_adic_log5, valuation


@lru_cache(maxsize=128)
def _cached_portrait(k: int, g: int, target_a: int) -> dict[str, list[int]]:
    """LRU-cached basin portrait keyed by (k, g, target_a)."""
    explorer = BasinExplorer(k, g, target_a)
    result: dict[str, list[int]] = {"converged": [], "cycle": [], "diverged": []}
    for e0 in range(explorer.N):
        fate, val, _ = explorer.classify(e0)
        result[fate].append(e0)
    return result


class BasinExplorer:
    """
    Explore the Newton-basin landscape for target a modulo 2^k.

    Parameters
    ----------
    k : int
        Bit precision (≥ 3).
    g : int
        Generator; must satisfy g ≡ 5 (mod 8).
    target_a : int
        The target integer for which to solve 5^e ≡ target_a (mod 2^k).
    """

    def __init__(self, k: int, g: int, target_a: int) -> None:
        if k < 3:
            raise ValueError("k must be ≥ 3")
        self.k = k
        self.g = g
        self.mask = bitmask(k)
        self.N = 1 << (k - 2)

        self._target_a = target_a
        target = target_a & self.mask
        alpha = (target >> 1) & 1

        # α=1 fix: for targets ≡ 3 mod 4, solve 5^e ≡ -a instead of 5^e ≡ a
        if alpha:
            self.a = (-target) & self.mask
            self.alpha = 1
        else:
            self.a = target
            self.alpha = 0

        self.L = two_adic_log5(k) >> 2

    def newton_step(self, e: int) -> int:
        """Single Newton iteration for the map e ↦ e - f(e)/f'(e)."""
        mask = self.mask
        exp_mask = self.N - 1

        g_val = pow(5, e, 1 << self.k)
        f_val = (g_val - self.a) & mask
        df_val = (g_val * self.L) & exp_mask
        df_inv = modinv_newton(df_val, self.k - 2)
        delta = ((f_val >> 2) * df_inv) & exp_mask
        return (e - delta) & exp_mask

    def _trajectory(
        self, e0: int, max_steps: int = 64, track_period: bool = False
    ) -> tuple[str, int, list[int]]:
        """
        Core classification: run Newton iteration from seed e0.

        Parameters
        ----------
        track_period : bool
            If True, return cycle period instead of cycle point.
        """
        seen: dict[int, int] = {}
        e = e0
        path = [e]

        for step in range(max_steps):
            if track_period:
                if e in seen:
                    period = step - seen[e]
                    return ("cycle", period, path)
            else:
                if e in seen:
                    return ("cycle", e, path)
            seen[e] = step

            e_next = self.newton_step(e)
            path.append(e_next)

            if e_next == e:
                if pow(5, e, 1 << self.k) == self.a:
                    return ("converged", step + 1 if track_period else e, path)
                else:
                    period_or_point = 1 if track_period else e
                    return ("cycle", period_or_point, path)

            e = e_next

        return ("diverged", max_steps if track_period else e, path)

    def classify(self, e0: int, max_steps: int = 64) -> tuple[str, int, list[int]]:
        """
        Classify seed e0 by Newton convergence fate.

        Returns
        -------
        (fate, value, path) where fate is one of:
            'converged' — converges to the true root
            'cycle'     — caught in a ghost cycle (fixed point ≠ true root)
            'diverged'  — escaped the exponent domain
        """
        return self._trajectory(e0, max_steps, track_period=False)

    def portrait(self) -> dict[str, list[int]]:
        """Full basin portrait over all seeds in the exponent domain (cached)."""
        return _cached_portrait(self.k, self.g, self._target_a)

    def fate_vector(self) -> list[int]:
        """Compact encoding: 0=converged, 1=cycle, 2=diverged."""
        vec = []
        for e0 in range(self.N):
            fate, _, _ = self.classify(e0)
            vec.append(0 if fate == "converged" else (1 if fate == "cycle" else 2))
        return vec

    def full_portrait(self) -> dict[str, Any]:
        """
        Full basin portrait with cycle period details.

        Returns
        -------
        {
            'converged': [(seed, steps_to_converge), ...],
            'cycles': {period: [(seed, cycle_path), ...], ...},
            'diverged': [seed, ...],
        }
        """
        result: dict[str, Any] = {
            "converged": [],
            "cycles": {},
            "diverged": [],
        }
        for e0 in range(self.N):
            fate, val, path = self._classify_with_period(e0)
            if fate == "converged":
                result["converged"].append((e0, val))
            elif fate == "cycle":
                period = val
                cycle_path = path[-period:] if period > 0 else path[-1:]
                result["cycles"].setdefault(period, []).append((e0, cycle_path))
            else:
                result["diverged"].append(e0)
        return result

    def portrait_matrix(self) -> list[int]:
        """
        Fate encoding per seed (matches 2-Adic-Newton-Dynamics API).

        Returns
        -------
        list[int] where:
            0       → converged
            -period → cycle of that period
            -1      → diverged
        """
        fates = []
        for e0 in range(self.N):
            fate, val, _ = self._classify_with_period(e0)
            if fate == "converged":
                fates.append(0)
            elif fate == "cycle":
                fates.append(-val)
            else:
                fates.append(-1)
        return fates

    def _classify_with_period(self, e0: int, max_steps: int = 64) -> tuple[str, int, list[int]]:
        """Classify seed and return cycle period (not cycle point).

        Delegates to the shared _trajectory method.
        """
        return self._trajectory(e0, max_steps, track_period=True)


def precision_sweep(
    k_min: int, k_max: int, g: int = 5, target_e: int = 2
) -> list[tuple[int, float]]:
    """
    Sweep precision k to find where ghost seeds first appear.

    Parameters
    ----------
    k_min, k_max : int — range of bit-precisions
    g : int — generator (default 5)
    target_e : int — exponent of target (a = g^{target_e} mod 2^k) (default 2)

    Returns
    -------
    List of (k, ghost_fraction) pairs.
    """
    results = []
    first_ghost_k = None
    for k in range(k_min, k_max + 1):
        a = pow(g, target_e, 1 << k)
        explorer = BasinExplorer(k, g, a)
        portrait = explorer.portrait()
        n_ghosts = len(portrait["cycle"])
        frac = n_ghosts / explorer.N if explorer.N > 0 else 0.0
        results.append((k, frac))
        if frac > 0 and first_ghost_k is None:
            first_ghost_k = k
    if first_ghost_k is not None:
        print(f"  Critical precision: k = {first_ghost_k}  (ghost fraction > 0)")
    return results


class LayerGhostDiagnosticV2:
    """
    Per-layer 2-adic ghost diagnostic for a weight matrix.

    Classifies each weight by its α-sector and computes v₂(e)
    statistics for quick assessment of ghost vulnerability.

    Parameters
    ----------
    k : int — bit precision (≥ 3)
    g : int — generator (default 5)
    max_iter : int — max Newton iterations (default 100)
    """

    def __init__(self, k: int, g: int = 5, max_iter: int = 100) -> None:
        self.k = k
        self.g = g
        self.max_iter = max_iter

    def diagnostic_matrix(
        self,
        weights: np.ndarray,
    ) -> tuple[np.ndarray, float, float, float, float]:
        """
        Analyse each weight in weights and return per-layer ghost stats.

        Returns
        -------
        (fate, conv_ratio, ghost_ratio, mean_e, v2_e)
        """
        flat = weights.ravel()
        fate = np.zeros(len(flat), dtype=np.int32)
        conv_count = 0
        ghost_count = 0
        e_vals: list[int] = []
        v2_vals: list[int] = []

        for i, w in enumerate(flat):
            w_int = int(w) & ((1 << self.k) - 1)
            if w_int & 1 == 0:
                fate[i] = -1  # even, skip
                continue
            result = two_adic_dlog(w_int, self.k)
            if result is None:
                fate[i] = -1
                continue
            alpha, e_true = result
            if alpha == 0:
                fate[i] = 0  # stable sector
                conv_count += 1
            else:
                fate[i] = 1  # ghost sector
                ghost_count += 1
            e_vals.append(e_true)
            if e_true != 0:
                v2_val = valuation(e_true)
                v2_vals.append(self.k if v2_val is None else v2_val)
            else:
                v2_vals.append(self.k)

        total = conv_count + ghost_count
        conv_ratio = conv_count / total if total > 0 else 0.0
        ghost_ratio = ghost_count / total if total > 0 else 0.0
        mean_e = float(np.mean(e_vals)) if e_vals else 0.0
        mean_v2_e = float(np.mean(v2_vals)) if v2_vals else 0.0

        return fate.reshape(weights.shape), conv_ratio, ghost_ratio, mean_e, mean_v2_e


class GhostHunt:
    """
    Sweep precision to find critical thresholds where ghosts appear.

    Provides both single-target precision sweeps and full weight-matrix
    quantization-cliff analysis.
    """

    def __init__(self, g: int = 5, max_iter: int = 100) -> None:
        self.g = g
        self.max_iter = max_iter

    def precision_threshold_sweep(self, k_min: int, k_max: int, target_e: int) -> None:
        """
        Print precision sweep table for a single target.

        Parameters
        ----------
        k_min, k_max : int — range of bit-precisions
        target_e : int — exponent of target (a = g^{target_e} mod 2^k)
        """
        print(f"Precision Sweep  (g={self.g}, target=5^{target_e})")
        print(f"{'k':>4}  {'N':>6}  {'Ghosts':>8}  {'Frac':>8}  {'Status':>10}")
        print("-" * 42)

        for k in range(k_min, k_max + 1):
            a = pow(self.g, target_e, 1 << k)
            explorer = BasinExplorer(k, self.g, a)
            portrait = explorer.portrait()
            n_ghosts = len(portrait["cycle"])
            frac = n_ghosts / explorer.N if explorer.N > 0 else 0.0
            status = "GHOST!" if frac > 0 else "OK"
            print(f"{k:>4}  {explorer.N:>6}  {n_ghosts:>8}  {frac:>8.4f}  {status:>10}")

    def quantization_cliff(
        self, weights: np.ndarray, k_min: int = 4, k_max: int = 16
    ) -> dict[int, float]:
        """
        Measure ghost density vs precision k for a weight matrix.

        Returns dict mapping k → ghost_fraction.
        """
        results: dict[int, float] = {}
        for k in range(k_min, k_max + 1):
            diag = LayerGhostDiagnosticV2(k, self.g, self.max_iter)
            _, _, ghost_ratio, _, _ = diag.diagnostic_matrix(weights)
            results[k] = ghost_ratio
        return results
