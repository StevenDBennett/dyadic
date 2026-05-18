"""
padic_roots.py
--------------
Multi-order p-adic root finding for x^3 ≡ a (mod p^k), p ≠ 2, 3.

Implements Newton (order 2), Halley (order 3), composed-Newton
(order 4), and triple-composed Newton (order 8) iterations.

The convergence law in Z_p is exact: v_p(x_n - x*) = m^n · v_p(x_0 - x*)
where m is the method order.  At finite precision mod p^k the valuation
saturates at k.
"""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Callable

import numpy as np

# ── p-adic helpers (self-contained for general primes) ──────────────────────


def _vp(n: int, p: int) -> int | None:
    """p-adic valuation for prime p. Returns None for zero."""
    if n == 0:
        return None
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def _modinv(a: int, pk: int) -> int | None:
    """Modular inverse modulo pk (pk need not be a prime power)."""
    # Extended Euclidean algorithm
    g, x, y = _ext_gcd(a, pk)
    if g != 1:
        return None
    return x % pk


def _ext_gcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = _ext_gcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def _pk(p: int, k: int) -> int:
    return int(p**k)


# ── Hensel lifting ──────────────────────────────────────────────────────────


def lift_root(a: int, p: int, k: int) -> int | None:
    """
    Hensel lift a cube root from mod p to mod p^k.

    Returns x such that x^3 ≡ a (mod p^k), or None if no root exists.
    Requires p ≠ 2, 3 and a not divisible by p.
    """
    # Find root mod p by brute force
    x0 = None
    for t in range(p):
        if (t * t * t) % p == a % p:
            x0 = t
            break
    if x0 is None:
        return None

    # Hensel lift
    x = x0
    mod = p
    for _ in range(1, k):
        mod *= p
        f = (x * x * x - a) % mod
        df = (3 * x * x) % mod
        inv_df = _modinv(df, mod)
        if inv_df is None:
            return None
        x = (x - f * inv_df) % mod

    return x


# ── Step functions ──────────────────────────────────────────────────────────


def newton_step(x: int, a: int, pk: int) -> int:
    """Newton (order 2): x_{n+1} = x_n - f/f'."""
    f = (x * x * x - a) % pk
    df = (3 * x * x) % pk
    inv_df = _modinv(df, pk)
    if inv_df is None:
        return x
    return (x - f * inv_df) % pk


def halley_step(x: int, a: int, pk: int) -> int:
    """Halley (order 3): x_{n+1} = x_n - 2f f' / (2f'^2 - f f'')."""
    f = (x * x * x - a) % pk
    df = (3 * x * x) % pk
    ddf = (6 * x) % pk
    num = (2 * f * df) % pk
    denom = (2 * df * df - f * ddf) % pk
    inv_denom = _modinv(denom, pk)
    if inv_denom is None:
        return x
    return (x - num * inv_denom) % pk


def newton2_step(x: int, a: int, pk: int) -> int:
    """Composed-2-Newton (effective order 4)."""
    x1 = newton_step(x, a, pk)
    return newton_step(x1, a, pk)


def newton3_step(x: int, a: int, pk: int) -> int:
    """Composed-3-Newton (effective order 8)."""
    x1 = newton_step(x, a, pk)
    x2 = newton_step(x1, a, pk)
    return newton_step(x2, a, pk)


# ── Convergence analysis ────────────────────────────────────────────────────


def convergence_profile(
    x0: int,
    a: int,
    p: int,
    k: int,
    step_fn: Callable[[int, int, int], int],
    x_true: int | None = None,
) -> list[int]:
    """
    Track v_p(x_n - x*) at each Newton iteration.

    Returns list of valuations at each step.
    """
    if x_true is None:
        x_true = lift_root(a, p, k)
        if x_true is None:
            return []

    pk = _pk(p, k)
    x = x0
    init_v = _vp(abs(x - x_true), p)
    profile: list[int] = [k] if init_v is None else [init_v]

    for _ in range(k + 1):
        x_new = step_fn(x, a, pk)
        diff = abs(x_new - x_true)
        if diff == 0:
            profile.append(k)  # saturated
            break
        v_val = _vp(diff, p)
        assert v_val is not None  # diff > 0, so valuation is finite
        profile.append(v_val)
        if v_val >= k:
            break
        x = x_new

    return profile


def compare_methods(
    p: int,
    k: int,
    n_trials: int = 20,
    seed: int | None = None,
) -> dict[str, float]:
    """
    Run all methods and report convergence rates.

    Returns dict with mean final v_p per method.
    """
    if seed is not None:
        random.seed(seed)
    methods: dict[str, Callable[[int, int, int], int]] = {
        "Newton (ord 2)": newton_step,
        "Halley (ord 3)": halley_step,
        "Comp-Newton (ord 4)": newton2_step,
        "Comp×3 (ord 8)": newton3_step,
    }
    results: dict[str, list[float]] = {name: [] for name in methods}

    pk = _pk(p, k)
    for _ in range(n_trials):
        a = random.randrange(2, min(pk, 10000))
        if a % p == 0:
            continue
        x_true = lift_root(a, p, k)
        if x_true is None:
            continue
        x0 = random.randrange(1, pk)

        for name, step_fn in methods.items():
            prof = convergence_profile(x0, a, p, k, step_fn, x_true)
            if prof:
                results[name].append(float(prof[-1]) if prof[-1] is not None else float(k))

    summary = {name: float(np.mean(vals)) if vals else 0.0 for name, vals in results.items()}
    return summary


def verify_order(
    primes: list[int] | None = None,
    k: int = 8,
    n_trials: int = 10,
    seed: int | None = None,
) -> dict[str, dict[int, float]]:
    """
    Verify that v_p(x_{n+1} - x*) / v_p(x_n - x*) = m (method order).

    Returns nested dict: method -> {trial_index: observed_order}.
    """
    if seed is not None:
        random.seed(seed)
    if primes is None:
        primes = [5, 7, 11, 13]

    methods: dict[str, Callable[[int, int, int], int]] = {
        "Newton": newton_step,
        "Halley": halley_step,
    }
    results: dict[str, dict[int, float]] = {name: {} for name in methods}

    for p in primes:
        pk = _pk(p, k)
        for trial in range(n_trials):
            a = random.randrange(2, min(pk, 5000))
            if a % p == 0:
                continue
            x_true = lift_root(a, p, k)
            if x_true is None:
                continue
            x0 = random.randrange(2, pk)
            if x0 % p == 0:
                continue

            for name, step_fn in methods.items():
                prof = convergence_profile(x0, a, p, k, step_fn, x_true)
                if len(prof) >= 3:
                    ratio = prof[-1] / prof[-2] if prof[-2] > 0 else 0
                    results[name][trial] = ratio

    return results


def newton_correction_uniformity(
    p: int,
    k: int,
    n_seeds: int = 1000,
    seed: int | None = None,
) -> dict[str, float]:
    """
    Test whether first-step Newton corrections are uniformly
    distributed modulo p (chi-square goodness-of-fit test).

    Reference: the first Newton correction term is
    delta = (x^3 - a) / (3 x^2) mod p, and the claim is this is
    uniform over F_p for random (x, a).
    """
    if seed is not None:
        random.seed(seed)
    pk = _pk(p, k)
    corrections = []

    for _ in range(n_seeds):
        a = random.randrange(2, min(pk, 10000))
        if a % p == 0:
            continue
        x = random.randrange(1, pk)
        if x % p == 0:
            continue

        f = (x * x * x - a) % pk
        df = (3 * x * x) % pk
        inv_df = _modinv(df, pk)
        if inv_df is None:
            continue
        delta = (f * inv_df) % pk
        corrections.append(delta % p)

    n = len(corrections)
    expected = n / p if n > 0 else 0.0
    counts = Counter(corrections)
    chi2 = (
        sum(((counts.get(r, 0) - expected) ** 2) / expected for r in range(p))
        if expected > 0
        else 0.0
    )
    return {
        "chi2_stat": chi2,
        "n_samples": n,
        "df": p - 1,
    }
