"""
isometry.py
-----------
Three capstone theorems tying together the Tidal Coordinates framework.

T6-a  Exponential Map Isometry
    v₂(5^e - 1) = v₂(e) + 2  for all e ≠ 0.
    The exponential map e ↦ 5^e scales 2-adic distances by exactly
    factor 4 (via the 2-adic Taylor series).

T6-c  Trace-mod-p Independence
    In GL(2) holonomies over Z/(2^k·p)Z, the 2-adic α-sector of
    det(H) and Tr(H) mod p are statistically independent.
"""

from __future__ import annotations

import random
from collections import Counter

import numpy as np
from dyadic_core import bitmask, two_adic_dlog, valuation

from .nonabelian import NonAbelianCRTDual

# ── T6-a: Exponential Map Isometry ──────────────────────────────────────────


def verify_isometry(k: int, n_trials: int = 200) -> dict[str, float]:
    """
    Verify v₂(5^e - 1) = v₂(e) + 2 for non-zero e.

    Returns pass rate and number of trials.
    """
    passed = 0
    for _ in range(n_trials):
        e = random.randrange(1, 1 << (k - 2))
        diff = (pow(5, e, 1 << k) - 1) & bitmask(k)
        v2_diff: int
        if diff == 0:
            v2_diff = k
        else:
            v2_diff_val = valuation(diff)
            v2_diff = k if v2_diff_val is None else v2_diff_val
        assert e != 0
        v2_e_val = valuation(e)
        assert v2_e_val is not None  # e > 0
        expected = v2_e_val + 2
        if v2_diff == expected or v2_diff >= k:
            passed += 1
    return {"pass_rate": passed / n_trials, "n_trials": n_trials}


def isometry_pair_test(k: int, n_trials: int = 200) -> dict[str, float]:
    """
    Verify v₂(5^{e1} - 5^{e2}) = v₂(e1 - e2) + 2 for e1 ≠ e2.
    """
    passed = 0
    collected = 0
    while collected < n_trials:
        e1 = random.randrange(0, 1 << (k - 2))
        e2 = random.randrange(0, 1 << (k - 2))
        if e1 == e2:
            continue
        collected += 1
        g1 = pow(5, e1, 1 << k)
        g2 = pow(5, e2, 1 << k)
        diff = (g1 - g2) & bitmask(k)
        v2_diff: int
        if diff != 0:
            v2_diff_val = valuation(diff)
            v2_diff = k if v2_diff_val is None else v2_diff_val
        else:
            v2_diff = k
        e_diff = (e1 - e2) & ((1 << (k - 2)) - 1)
        v2_e: int
        if e_diff != 0:
            v2_e_val = valuation(e_diff)
            v2_e = k if v2_e_val is None else v2_e_val
        else:
            v2_e = k
        expected = v2_e + 2
        if v2_diff == expected or v2_diff >= k:
            passed += 1
    return {"pass_rate": passed / n_trials, "n_trials": n_trials}


def isometry_summary(k: int) -> str:
    """Print the full bidirectional conditioning picture."""
    n_cases = 1 << (k - 2)
    lines = [
        f"Exponential Map Isometry Summary (k={k}, N={n_cases})",
        "  v₂(5^e - 1)     = v₂(e) + 2     (T6-a: isometry)",
        "  n*(s)           = ceil(log₂(s)) - 1  (T2: separation)",
        "",
        "  Together these imply the bidirectional Lipschitz condition:",
        "    v₂(e1 - e2) + 2 ≤ v₂(5^e1 - 5^e2) ≤ k",
        "  The Newton map is a 2-adic contraction with factor 1/2.",
    ]
    return "\n".join(lines)


# ── T6-c: Trace-mod-p Independence ──────────────────────────────────────────


def trace_alpha_independence(
    k: int, p: int, cycle_length: int = 4, n_cycles: int = 100
) -> dict[str, float]:
    """
    Chi-square test: α(det(H)) and Tr(H) mod p are independent.

    Returns chi2_stat, p_value, and degrees of freedom.
    """
    nc = NonAbelianCRTDual(k, p)
    observations: list[tuple[int, int]] = []

    for _ in range(n_cycles):
        mats = [
            [[random.randrange(0, nc.mod_full) for _ in range(2)] for _ in range(2)]
            for _ in range(cycle_length)
        ]
        inv = nc.invariants(mats)
        if inv["det_even"]:
            continue
        alpha = inv["alpha_det"]
        if alpha is None:
            continue
        observations.append((alpha, inv["trace_modp"]))

    # Contingency table: alpha (0,1) × trace_modp (0..p-1)
    table = Counter(observations)
    n = len(observations)
    alpha_counts = Counter(a for a, _ in observations)
    trace_counts = Counter(t for _, t in observations)

    chi2 = 0.0
    for alpha in (0, 1):
        for trace in range(p):
            observed = table.get((alpha, trace), 0)
            expected = alpha_counts[alpha] * trace_counts[trace] / n if n > 0 else 0
            if expected > 0:
                chi2 += (observed - expected) ** 2 / expected

    df = (2 - 1) * (p - 1)
    return {
        "chi2_stat": chi2,
        "df": df,
        "n_samples": n,
    }


def trace_exponent_independence(
    k: int, p: int, cycle_length: int = 4, n_cycles: int = 100
) -> dict[str, float]:
    """
    ANOVA F-test: e(det(H)) vs Tr(H) mod p.

    Returns F-statistic and p-value from scipy if available,
    otherwise returns effect size estimate.
    """
    nc = NonAbelianCRTDual(k, p)
    groups: dict[int, list[float]] = {}

    for _ in range(n_cycles):
        mats = [
            [[random.randrange(0, nc.mod_full) for _ in range(2)] for _ in range(2)]
            for _ in range(cycle_length)
        ]
        inv = nc.invariants(mats)
        if inv["det_even"]:
            continue
        e_det = inv["e_det"]
        if e_det is None:
            continue
        trace = inv["trace_modp"]
        if trace not in groups:
            groups[trace] = []
        groups[trace].append(float(e_det))

    # One-way ANOVA
    all_vals = [v for vals in groups.values() for v in vals]
    grand_mean = float(np.mean(all_vals)) if all_vals else 0.0

    ss_between = 0.0
    ss_within = 0.0
    for vals in groups.values():
        n_group = len(vals)
        group_mean = float(np.mean(vals))
        ss_between += n_group * (group_mean - grand_mean) ** 2
        ss_within += sum((v - group_mean) ** 2 for v in vals)

    df_between = len(groups) - 1
    df_within = len(all_vals) - len(groups)

    if df_within > 0 and df_between > 0:
        ms_between = ss_between / df_between
        ms_within = ss_within / df_within
        f_stat = ms_between / ms_within if ms_within > 0 else 0.0
    else:
        f_stat = 0.0

    return {
        "F_stat": f_stat,
        "df_between": df_between,
        "df_within": df_within,
        "n_samples": len(all_vals),
    }


def exponentvaluation_profile(k: int, n_samples: int = 500) -> dict[int, float]:
    """
    Profile the distribution of v₂(e_true) across random odd weights.

    This is the genuinely graded stability measure (the ghost
    ratio encodes only the α-bit).
    """
    v2_counts: dict[int, int] = {}
    for _ in range(n_samples):
        w = random.randrange(1, 1 << k)
        if w & 1 == 0:
            continue
        w_adj = w if (w & 3) == 1 else (-w) & bitmask(k)
        result = two_adic_dlog(w_adj, k)
        if result is not None:
            _, e_true = result
            v2: int
            if e_true != 0:
                v2_val = valuation(e_true)
                v2 = k - 2 if v2_val is None else v2_val
            else:
                v2 = k - 2
            v2_counts[v2] = v2_counts.get(v2, 0) + 1

    total = sum(v2_counts.values())
    return {v: count / total if total > 0 else 0.0 for v, count in sorted(v2_counts.items())}
