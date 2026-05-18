"""
thermodynamics.py
-----------------
Graded 2-adic weight stability diagnostics for collections of integers.

The primary quantity is v₂(e_true) — the 2-adic valuation of the
true discrete-logarithm exponent.  This is the genuinely graded
stability measure: higher v₂(e_true) means the weight lies deeper
in the 2-adic congruence filtration, hence more stable under
quantization.

Mersenne numbers w = 2ⁿ - 1 are maximally fragile with
v₂(e_true) = n - 2 and cliff at k* = n + 2.

References
----------
- Mersenne Ghost Theorem (mersenne.py)
- Trajectory Separation Theorem (separation.py)
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
from dyadic_core import bitmask, two_adic_dlog, valuation

from dyadic_math.basin import BasinExplorer


class SeedThermodynamics:
    """
    Analysis of 2-adic weight thermodynamics for a collection of integers.

    Two API styles are supported:

    1. **Coordinate analysis**: given a weight matrix,
        compute comprehensive 2-adic statistics via ``analyse(W)``.

    2. **Precision-sweep analysis**: use the classmethod
       ``from_precision_sweep(weights, k_range)`` to create an instance
       configured for ghost-cliff detection, then call ``compute()``.

    Parameters
    ----------
    k : int
        Bit precision for coordinate analysis (default 16).
    g : int
        Generator; must satisfy g ≡ 5 (mod 8).  Default 5.
    """

    def __init__(self, k: int = 16, g: int = 5) -> None:
        self.k = k
        self.g = g
        self.mask = bitmask(k)
        self.N = 1 << (k - 2)
        self._weights: np.ndarray | None = None
        self._k_range: range | None = None
        self._profiles: dict[int, dict[int, float]] | None = None
        self._cliffs: dict[int, int | None] | None = None

    @classmethod
    def from_precision_sweep(
        cls,
        weights: np.ndarray,
        k_range: range,
        k: int = 16,
        g: int = 5,
    ) -> SeedThermodynamics:
        """
        Create an instance configured for precision-sweep analysis.

        Parameters
        ----------
        weights : np.ndarray
            Integer weight matrix to analyse.
        k_range : range
            Range of bit-precisions to sweep.
        k : int
            Bit precision for the underlying DualNumber arithmetic.
        g : int
            Generator; must satisfy g ≡ 5 (mod 8).

        Returns
        -------
        SeedThermodynamics
            Instance ready for ``compute()``.
        """
        instance = cls(k=k, g=g)
        instance._weights = weights
        instance._k_range = k_range
        return instance

    # ── Coordinate analysis (v1-style) ──────────────────────────────────────

    def weight_coordinates(self, w: int) -> tuple[float, int, int]:
        """
        2-adic dual coordinates (v, alpha, e) for a single weight.

        Returns (v, alpha, e) where v can be float('inf') for zero.
        """
        wbm = w & self.mask
        if wbm == 0:
            return (float("inf"), 0, 0)
        v = valuation(wbm)
        if v is None or v >= self.k:
            return (float("inf"), 0, 0)
        odd = (wbm >> v) & self.mask
        result = two_adic_dlog(odd, self.k)
        if result is None:
            return (float(v), 0, 0)
        alpha, e = result
        return (float(v), alpha, e)

    def analyse(self, weights: np.ndarray) -> dict[str, float]:
        """
        Compute comprehensive 2-adic statistics for weight matrix weights.

        Returns a dict with keys:
            alpha_fraction, mean_v2_e, std_v2_e, mean_e_norm,
            n_mersenne_zone, cliff_risk
        """
        alpha_frac, v2_e_vals, e_norms = 0.0, list[int](), []
        n_odd = 0

        for w in weights.flat:
            w_int = int(w)
            if w_int & 1 == 0:
                continue
            n_odd += 1
            _, alpha, _ = self.weight_coordinates(w_int)
            alpha_frac += alpha
            odd = abs(w_int) & self.mask
            result = two_adic_dlog(odd, self.k)
            if result is not None:
                _, e_true = result
                v2_e: int
                if e_true != 0:
                    v2_val = valuation(e_true)
                    v2_e = self.k if v2_val is None else v2_val
                else:
                    v2_e = self.k
                v2_e_vals.append(v2_e)
                e_norms.append(e_true / self.N)

        alpha_frac = alpha_frac / n_odd if n_odd > 0 else 0.0
        mean_v2_e = float(np.mean(v2_e_vals)) if v2_e_vals else 0.0
        std_v2_e = float(np.std(v2_e_vals)) if v2_e_vals else 0.0
        mean_e_norm = float(np.mean(e_norms)) if e_norms else 0.0
        n_mersenne_zone = sum(1 for v in v2_e_vals if v is not None and v >= self.k // 2)

        return {
            "alpha_fraction": alpha_frac,
            "mean_v2_e": mean_v2_e,
            "std_v2_e": std_v2_e,
            "mean_e_norm": mean_e_norm,
            "n_mersenne_zone": n_mersenne_zone,
            "cliff_risk": n_mersenne_zone / len(v2_e_vals) if v2_e_vals else 0.0,
        }

    def mersenne_cliff_score(self, w: int) -> int | None:
        """
        Minimum precision k* needed for stable quantization of weight w.

        Returns k* = v2(e_true) + 2, or None if w is even.
        """
        if w & 1 == 0:
            return None
        odd = abs(w) & self.mask
        result = two_adic_dlog(odd, self.k)
        if result is None:
            return None
        _, e_true = result
        if e_true == 0:
            return self.k
        v2_e = valuation(e_true)
        if v2_e is None:
            return self.k
        return v2_e + 2

    def compare_to_random(self, weights: np.ndarray, n_samples: int = 1000) -> dict[str, float]:
        """
        Z-score comparison of statistics against random baseline.

        Returns dict with z_alpha, z_v2_e, z_cliff_risk.
        """
        stats = self.analyse(weights)
        alphas, v2s, cliffs = [], [], []

        for _ in range(n_samples):
            random_weights = np.random.randint(0, 1 << self.k, size=weights.shape, dtype=np.int64)
            rstats = self.analyse(random_weights)
            alphas.append(rstats["alpha_fraction"])
            v2s.append(rstats["mean_v2_e"])
            cliffs.append(rstats["cliff_risk"])

        mu_a, std_a = float(np.mean(alphas)), float(np.std(alphas)) or 1.0
        mu_v, std_v = float(np.mean(v2s)), float(np.std(v2s)) or 1.0
        mu_c, std_c = float(np.mean(cliffs)), float(np.std(cliffs)) or 1.0

        return {
            "z_alpha": (stats["alpha_fraction"] - mu_a) / std_a,
            "z_v2_e": (stats["mean_v2_e"] - mu_v) / std_v,
            "z_cliff_risk": (stats["cliff_risk"] - mu_c) / std_c,
        }

    # ── Precision-sweep analysis (v2-style) ─────────────────────────────────

    def __call__(self, weights: np.ndarray, k_range: range) -> SeedThermodynamics:
        """Deprecated: use ``SeedThermodynamics.from_precision_sweep()`` instead."""
        warnings.warn(
            "Calling SeedThermodynamics(weights, k_range) is deprecated. "
            "Use SeedThermodynamics.from_precision_sweep(weights, k_range) instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.from_precision_sweep(weights, k_range, k=self.k, g=self.g)

    def compute(self) -> None:
        """Lazy computation of ghost density profiles per weight."""
        if self._profiles is not None:
            return
        assert self._weights is not None, "call from_precision_sweep first"
        assert self._k_range is not None, "call from_precision_sweep first"
        self._profiles = {}
        self._cliffs = {}
        weights_flat = self._weights.ravel()
        for idx in range(len(weights_flat)):
            w = int(weights_flat[idx])
            profile: dict[int, float] = {}
            cliff: int | None = None
            for k in self._k_range:
                if w & 1 == 0:
                    profile[k] = 0.0
                    continue
                try:
                    be = BasinExplorer(k, self.g, w)
                except Exception:
                    profile[k] = 0.0
                    continue
                p = be.portrait()
                frac = len(p["cycle"]) / be.N if be.N > 0 else 0.0
                profile[k] = frac
                if frac > 0 and cliff is None:
                    cliff = k
            self._profiles[idx] = profile
            self._cliffs[idx] = cliff

    @property
    def profiles(self) -> dict[int, dict[int, float]]:
        if self._profiles is None:
            self.compute()
        assert self._profiles is not None
        return self._profiles

    @property
    def cliffs(self) -> dict[int, int | None]:
        if self._cliffs is None:
            self.compute()
        assert self._cliffs is not None
        return self._cliffs

    def cliff_histogram(self) -> dict[int, int]:
        """Distribution of cliff precisions across weights."""
        hist: dict[int, int] = {}
        for c in self.cliffs.values():
            if c is not None:
                hist[c] = hist.get(c, 0) + 1
        return hist

    def stable_weights(self, max_k: int | None = None) -> list[int]:
        """
        Indices of weights with no ghost cliff (or cliff >= max_k).
        """
        if max_k is not None:
            return [idx for idx, c in self.cliffs.items() if c is None or c >= max_k]
        return [idx for idx, c in self.cliffs.items() if c is None]

    def ghost_weights(self, max_k: int | None = None) -> list[int]:
        """
        Indices of weights that have a ghost cliff (below max_k).
        """
        if max_k is not None:
            return [idx for idx, c in self.cliffs.items() if c is not None and c < max_k]
        return [idx for idx, c in self.cliffs.items() if c is not None]

    def summary(self) -> dict[str, float]:
        """Aggregate statistics across all weights."""
        vals = [c for c in self.cliffs.values() if c is not None]
        return {
            "n_weights": len(self.cliffs),
            "n_stable": len(self.stable_weights()),
            "n_ghost": len(self.ghost_weights()),
            "ghost_fraction": len(self.ghost_weights()) / max(len(self.cliffs), 1),
            "mean_cliff": float(np.mean(vals)) if vals else 0.0,
            "min_cliff": min(vals) if vals else 0.0,
            "max_cliff": max(vals) if vals else 0.0,
        }

    def report(self) -> str:
        """ASCII report of ghost cliff statistics."""
        s = self.summary()
        lines = [
            "SeedThermodynamics Report",
            f"  Weights:        {s['n_weights']}",
            f"  Stable:         {s['n_stable']}",
            f"  Ghost:          {s['n_ghost']} ({s['ghost_fraction']:.1%})",
            f"  Mean cliff:     k = {s['mean_cliff']:.1f}",
            f"  Min cliff:      k = {s['min_cliff']:.0f}",
            f"  Max cliff:      k = {s['max_cliff']:.0f}",
        ]

        hist = self.cliff_histogram()
        if hist:
            bar_max = max(hist.values())
            bar_width = 30
            lines.append("  Cliff histogram:")
            for k_val in sorted(hist):
                count = hist[k_val]
                bar_len = max(1, int(count / bar_max * bar_width))
                bar = "\u2588" * bar_len
                lines.append(f"    k={k_val:2d}  {bar} {count}")

        return "\n".join(lines)
