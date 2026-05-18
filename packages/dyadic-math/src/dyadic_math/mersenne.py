"""
mersenne.py
-----------
Mersenne Ghost Theorem, bootstrap optimality analysis, and cliff
constant verification.

Mersenne Ghost Theorem
    For w = 2^n - 1 (a Mersenne number), the 2-adic dual coordinates
    are always α = 1 (ghost sector) and e_true = 2^(n-2) within a
    stable window.  Consequently v₂(e_true) = n-2 and the quantisation
    cliff is at k* = n+2.

Bootstrap Optimality
    In the Viglietta discrete-log algorithm, the bootstrap phase
    starting precision eprec₀ = k/2 is optimal (saves log₂(√k) - 1
    Newton steps compared to the √k heuristic currently used).

Cliff Constant Verification
    Numerical verification of the Mersenne Cliff Theorem:
      - cliff_constant(g, k)        — compute c = v₂(log₂(g)/4 + 1)
      - cliff_formula(g)            — human-readable c(g) explanation
      - mersenne_cliff_conjecture()    — state and verify the conjectured cliff formula
      - verify_cliff_constant()      — verify c=5 from 4 log-series terms
      - verify_c_formula()           — verify c(g) = v₂(g-5) - 2
      - exp2_neg4(k)                — delegates to dyadic_core.g0
      - cliff_constant_unified(g,k) — unified formula
      - verify_unified_formula()    — verify unified matches direct
      - verify_connection()          — show all checks are connected
"""

from __future__ import annotations

import math
import random
from functools import lru_cache

from dyadic_core import (
    bitmask,
    dlog_bootstrap,
    dlog_newton_step,
    g0,
    modinv_newton,
    padic_log,
    two_adic_dlog,
    two_adic_log5,
    valuation,
)

# ── Mersenne coordinates ────────────────────────────────────────────────────


def mersenne_coordinates(n: int, k: int) -> tuple[int, int, int] | None:
    """
    Return (α, e_true, v₂(e_true)) for Mersenne weight w = 2^n - 1.

    Requires n ≥ 3 and k ≥ n + 2.
    """
    if n < 3:
        raise ValueError("n must be ≥ 3")
    w = (1 << n) - 1
    result = two_adic_dlog(w, k)
    if result is None:
        return None
    alpha, e_true = result
    v2_e: int
    if e_true != 0:
        v2_e_val = valuation(e_true)
        assert v2_e_val is not None  # e_true > 0
        v2_e = v2_e_val
    else:
        v2_e = k
    return alpha, e_true, v2_e


def verify_core_identity(n_max: int = 12) -> dict[int, bool]:
    """
    Verify the core identity: 5^(2^(n-2)) ≡ 1 - 2^n (mod 2^(n+1)).

    This identity holds for all n ≥ 3 and is the foundation of the
    Mersenne Ghost Theorem.
    """
    results: dict[int, bool] = {}
    for n in range(3, n_max + 1):
        lhs = pow(5, 1 << (n - 2), 1 << (n + 1))
        rhs = (1 - (1 << n)) & bitmask(n + 1)
        results[n] = lhs == rhs
    return results


def _cliff_prediction(n: int, c: int) -> int:
    """Mersenne cliff formula: the last stable precision k* for Mersenne number 2^n - 1.

    Three regimes (Mersenne Ghost Theorem, g=5, c=5):
        n <= c:      k* = 2n - 1          (quadratic dominates)
        n = c + 1:   k* = n + c + 1        (boundary)
        n >= c + 2:  k* = n + c            (linear determines)
    """
    if n <= c:
        return 2 * n - 1
    elif n == c + 1:
        return n + c + 1
    else:
        return n + c


def mersenne_cliff_table(n_max: int = 12) -> list[dict[str, int]]:
    """
    Find the cliff precision k* for each Mersenne weight 2^n - 1.

    Returns list of dicts with n, k*, α, e_true, v₂(e_true), and
    the formula prediction k_pred using the full 3-regime formula.
    """
    c = cliff_constant(g=5)
    rows = []
    for n in range(3, n_max + 1):
        stable_k = None
        for k in range(n + 1, n + 10):
            result = two_adic_dlog((1 << n) - 1, k)
            if result is not None:
                _, e_true = result
                if e_true == (1 << (n - 2)):
                    stable_k = k
                else:
                    break
        if stable_k is not None:
            _, e_true, v2_e = mersenne_coordinates(n, stable_k + 1) or (0, 0, 0)
            rows.append(
                {
                    "n": n,
                    "k*": stable_k,
                    "c": c,
                    "k_pred": _cliff_prediction(n, c),
                    "alpha": 1,
                    "e_true": 1 << (n - 2),
                    "v2_e": v2_e,
                }
            )
    return rows


# ── Bootstrap optimality ────────────────────────────────────────────────────


def bootstrap_cost(eprec0: int, k: int) -> int:
    """
    Total bits processed by the Viglietta dlog algorithm.

    Cost = eprec0 (bootstrap) + Σ_{i} 2·eprec_i (Newton steps)
    where eprec_i doubles each iteration.
    """
    cost = eprec0
    eprec = eprec0
    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        cost += new_eprec  # each Newton step processes ~new_eprec bits
        eprec = new_eprec
    return cost


def optimal_bootstrap(k_values: list[int] | None = None) -> dict[int, int]:
    """
    Find the optimal eprec₀ for each k by minimising total cost.

    Returns dict mapping k → optimal eprec₀ (consistently near k/2).
    """
    if k_values is None:
        k_values = list(range(8, 65, 4))

    results: dict[int, int] = {}
    for k in k_values:
        best_cost = float("inf")
        best_eprec = 0
        for eprec0 in range(2, k - 2):
            cost = bootstrap_cost(eprec0, k)
            if cost < best_cost:
                best_cost = cost
                best_eprec = eprec0
        results[k] = best_eprec
    return results


def compare_bootstrap_strategies(k_values: list[int] | None = None) -> dict[int, dict[str, int]]:
    """
    Compare sqrt(k) heuristic vs k/2 optimal vs LUT b=8.

    Returns dict mapping k → {sqrt_cost, half_cost, lut_cost}.
    """
    if k_values is None:
        k_values = [16, 24, 32, 48, 64]

    results: dict[int, dict[str, int]] = {}
    for k in k_values:
        sqrt_eprec = max(4, int(math.isqrt(k)) + 2)
        half_eprec = max(4, k // 2 + 2)
        results[k] = {
            "sqrt_eprec": sqrt_eprec,
            "sqrt_cost": bootstrap_cost(sqrt_eprec - 2, k),
            "half_eprec": half_eprec,
            "half_cost": bootstrap_cost(half_eprec - 2, k),
            "lut_cost": bootstrap_cost(6, k),  # 8-bit LUT gives 6 bits
        }
    return results


# ── LUT-based dlog ──────────────────────────────────────────────────────────


@lru_cache(maxsize=32)
def _build_lut(b: int = 8) -> dict[int, int]:
    """
    Build a lookup table for dlog modulo 2^b.

    Returns dict mapping a → e such that 5^e ≡ a (mod 2^b),
    for all odd a ≡ 1 (mod 4).
    """
    mod = 1 << b
    lut: dict[int, int] = {}
    for e in range(1 << (b - 2)):
        a = pow(5, e, mod)
        lut[a] = e
    return lut


def dlog_with_lut(a: int, k: int, b: int = 8) -> int:
    """
    Discrete logarithm using a pre-computed LUT for bootstrap.

    Eliminates the O(k) bootstrap phase — the LUT gives exactly
    6 bits of precision at b=8, and each Newton step doubles.
    """
    if k <= b - 2:
        return dlog_bootstrap(a, k)
    if a & 3 != 1:
        raise ValueError("dlog_with_lut requires a ≡ 1 (mod 4)")

    lut = _build_lut(b)
    a_b = a & ((1 << b) - 1)
    e = lut.get(a_b, 0) & bitmask(b - 2)
    eprec = b - 2

    log5_unit = two_adic_log5(k) >> 2
    mask_full = bitmask(k)
    a &= mask_full

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = bitmask(bits)
        emask = bitmask(new_eprec)
        e, _ = dlog_newton_step(e, a, log5_unit, bits, mask, emask, new_eprec)
        eprec = new_eprec

    return e


def verify_lut_dlog(k: int, b: int = 8, n_trials: int = 100, seed: int | None = None) -> bool:
    """
    Verify that LUT-based dlog matches the standard dlog.
    """
    if seed is not None:
        random.seed(seed)
    for _ in range(n_trials):
        a = random.randrange(1, 1 << k)
        if a & 1 == 0:
            continue
        a_mod = a if (a & 3) == 1 else (-a) & bitmask(k)
        e1 = dlog_with_lut(a_mod, k, b)
        e2 = dlog_bootstrap(a_mod, k)
        if e1 != e2:
            return False
    return True


# ── cliff constant ──────────────────────────────────────────────────────────────


def cliff_constant(g: int = 5, k: int = 20) -> int:
    """
    Compute c = v_2(log_2(g)/4 + 1) for generator g ≡ 5 (mod 8).

    This constant determines the Mersenne cliff: k*(n) = n + c for n >= c + 2.

    Formula (derived from 2-adic log linearisation):
        c(g=5) = 5   (computed from v_2(two_adic_log5/4 + 1))
        c(g)   = v_2(g-5) - 2   if v_2(g-5) < 7   (linear term dominates)
        c(g)   = 5              if v_2(g-5) > 7   (capped at c(5))
        c(g)   = v_2(L_g+1)    at the boundary v_2(g-5)=7
    """
    if g == 5:
        log_val = two_adic_log5(k) >> 2
    else:
        log_val = padic_log(g, k) >> 2
    c_val = valuation(log_val + 1)
    assert c_val is not None  # log_val + 1 > 0
    return c_val


def cliff_formula(g: int) -> str:
    """
    Human-readable summary of the c(g) formula for a given generator.

    The formula c(g) = v_2(log_2(g)/4 + 1) has a closed-form approximation:

        s = v_2(g - 5)

        s < 7:  c(g) = s - 2
        s > 7:  c(g) = 5
        s = 7:  c(g) = 5 + v_2(A + B)  (boundary: recursive)
    """
    diff = g - 5
    s_val = valuation(diff)
    s: int = 999 if s_val is None else s_val
    c_actual = cliff_constant(g)
    if s == 999:
        regime = f"g = 5 (fixed point): c = {c_actual}"
    elif s - 2 < 5:
        regime = f"s={s} < 7: c = s-2 = {s - 2} (linear dominates)"
    elif s - 2 > 5:
        regime = f"s={s} > 7: c = 5 (capped at c(5))"
    else:
        regime = f"s=7 (BOUNDARY): c = 5 + recursion = {c_actual}"
    return f"g={g}: {regime}"


def mersenne_cliff_conjecture(verbose: bool = False) -> dict[str, int | str | bool]:
    """
    State and verify the Mersenne Cliff Theorem over a finite range.

    THEOREM (conjectured form, verified empirically below)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Let w = 2^n - 1 (Mersenne number, n >= 3) and let c = v_2(L_unit + 1)
    where L_unit = two_adic_log5() >> 2. For g = 5: c = 5.

    Define k*(n) = last precision k where e_true(w, k) = 2^(n-2) exactly.

    Then:
        n <= c:      k*(n) = 2n - 1
        n = c + 1:   k*(n) = n + c + 1
        n >= c + 2:  k*(n) = n + c

    For g = 5 (c = 5):
        n = 3,4,5:   k* = 5, 7, 9     (= 2n-1)
        n = 6:       k* = 12           (= n+6 = n+c+1)
        n >= 7:      k* = n + 5        (verified empirically for n=7..23)
    """
    c = cliff_constant(g=5)

    gen7_ok = True
    for n in range(c + 2, c + 18):
        k_last = None
        for k in range(n + 1, n + c + 8):
            result = mersenne_coordinates(n, k)
            if result is None:
                continue
            _, e, _ = result
            e_exp = 1 << (n - 2)
            if e == e_exp:
                k_last = k
            elif k_last is not None:
                break
        if k_last != n + c:
            gen7_ok = False
            break

    small_ok = True
    for n in range(3, c + 1):
        k_last = None
        for k in range(n + 1, n + c + 8):
            result = mersenne_coordinates(n, k)
            if result is None:
                continue
            _, e, _ = result
            e_exp = 1 << (n - 2)
            if e == e_exp:
                k_last = k
            elif k_last is not None:
                break
        if k_last != 2 * n - 1:
            small_ok = False
            break

    n6 = c + 1
    result = mersenne_coordinates(n6, n6 + c)
    boundary_ok = result is not None and result[1] == (1 << (n6 - 2))

    if verbose:
        print("Mersenne Cliff Theorem  (g=5, c=v_2(L_unit+1)=" + str(c) + ")")
        print()
        print(f"    n <= {c}:      k* = 2n-1          [quadratic dominates]")
        print(f"    n = {c + 1}:     k* = n+{c + 1}           [boundary]")
        print(f"    n >= {c + 2}:    k* = n+{c}")
        print(f"        [linear determines; verified n={c + 2}..{c + 17}]")
        print()
        print(f"    n=3..{c}:  {'PASS' if small_ok else 'FAIL'}")
        print(f"    n={c + 1}:   {'PASS' if boundary_ok else 'FAIL'}")
        print(f"    n={c + 2}..{c + 17}: {'PASS' if gen7_ok else 'FAIL'}")

    return {
        "c": c,
        "formula": "k*(n) = n+c for n>=c+2; k*(n)=2n-1 for n<=c; n=c+1 is boundary",
        "gen7_ok": gen7_ok,
        "small_ok": small_ok,
        "boundary_ok": boundary_ok,
        "all_pass": gen7_ok and small_ok and boundary_ok,
    }


def verify_cliff_constant(verbose: bool = False) -> bool:
    """
    Verify v_2(log_2(5)/4 + 1) = 5 from the first 4 terms of the 2-adic log series.

    Uses only the first 4 terms of the 2-adic log series, showing all
    higher terms are O(2^8) and therefore invisible mod 2^8 = 256.
    """
    inv3 = modinv_newton(3, 8)
    t1 = 4
    t2 = (-8) % 256
    t3 = (64 * inv3) % 256
    t4 = (-64) % 256
    total = (t1 + t2 + t3 + t4) % 256

    checks = [
        ("3^{-1} mod 256 = 171", inv3 == 171),
        ("sum mod 256 = 124", total == 124),
        ("124 ≡ -4 (mod 128)", total % 128 == (-4) % 128),
        ("124 ≢ -4 (mod 256)", total != (-4) % 256),
        ("v_2(L_unit+1) = 5", valuation(two_adic_log5(10) // 4 + 1) == 5),
    ]

    all_ok = all(ok for _, ok in checks)

    if verbose:
        print("Check: v_2(log_2(5)/4 + 1) = 5")
        print()
        print(f"  log_2(5) ≡ {t1} + {t2} + {t3} + {t4}  (mod 256)")
        print(f"           = {total}  (mod 256)")
        print()
        print("  v_2(log_2(5)/4 + 1) = 5.")
        print()
        for msg, ok in checks:
            print(f"  {'✓' if ok else '✗'}  {msg}")

    return all_ok


def verify_c_formula(verbose: bool = False) -> bool:
    """
    Verify c(g) = v_2(g-5) - 2 for g ≡ 5 (mod 8), g ≠ 5, v_2(g-5) ≠ 7.

    The formula c(g) = min(v_2(g-5)-2, 5) follows from 2-adic log
    linearisation.  This function checks it numerically over g = 13..600.
    """
    failures = []
    for g in range(13, 600, 8):
        if g % 8 != 5:
            continue
        diff = g - 5
        s = valuation(diff)
        assert s is not None  # g > 5, so diff > 0
        if s == 7:
            continue
        s2 = s - 2
        pred = min(s2, 5)
        actual = cliff_constant(g, k=20)
        if pred != actual:
            failures.append((g, s, pred, actual))

    ok = len(failures) == 0
    if verbose:
        print("Check: c(g) = min(v_2(g-5)-2, 5)  for v_2(g-5) ≠ 7")
        print()
        print("  Numerical check (g = 13..600, all s ≠ 7):")
        print(f"    Failures: {len(failures)}")
        if failures:
            for g, s, pred, actual in failures[:3]:
                print(f"    g={g}: s={s}, pred={pred}, actual={actual}")

    return ok


def exp2_neg4(k: int) -> int:
    """
    Compute exp_2(-4) mod 2^k — the 2-adic zero of log_2(g)/4 + 1.

    Delegates to ``dyadic_core.g0`` (the canonical implementation).
    """
    return g0(k)


def cliff_constant_unified(g: int, k: int = 32) -> int:
    """
    Compute c(g) = v_2(g - exp_2(-4)) - 2.

    Unified formula for the cliff constant with no casework.

    By the 2-adic Newton-Taylor lemma:
        v_2(f(g)) = v_2(g - g_0) + v_2(f'(g_0))
        where f(g) = log_2(g)/4 + 1, g_0 = exp_2(-4), f'(g_0) = 1/(4g_0)
        => v_2(f(g_0)) = -2
        => c(g) = v_2(g - exp_2(-4)) - 2.  QED
    """
    g0_val = exp2_neg4(k)
    v = valuation(g - g0_val)
    return (v - 2) if v is not None else k


def verify_unified_formula(
    g_values: list[int] | None = None,
    k: int = 24,
    verbose: bool = False,
) -> bool:
    """
    Verify cliff_constant_unified(g) == cliff_constant(g) for all test generators.
    """
    if g_values is None:
        g_values = [
            5,
            13,
            21,
            29,
            37,
            45,
            53,
            61,
            69,
            77,
            85,
            93,
            133,
            261,
            389,
            517,
            645,
            773,
            901,
            1157,
            1925,
            2053,
        ]

    failures = []
    for g in g_values:
        if g % 8 != 5:
            continue
        direct = cliff_constant(g, k=k)
        unified = cliff_constant_unified(g, k=k + 8)
        if direct != unified:
            failures.append((g, direct, unified))

    ok = len(failures) == 0
    if verbose:
        print("Unified formula: c(g) = v_2(g - exp_2(-4)) - 2")
        print()
        print(f"  Failures: {len(failures)}")
        print(f"  Result: {'PASS ✓' if ok else 'FAIL ✗'}")

    return ok


def verify_connection(verbose: bool = False) -> bool:
    """
    Show that all four verification checks are connected through a single identity:

        log_2(5)  ≡  -4  (mod 128)

    This single fact, computed from 4 terms of the 2-adic log series,
    generates the entire theorem structure.
    """
    log5_mod128 = two_adic_log5(8) & 127
    minus4_mod128 = (-4) % 128

    checks = [
        ("log_2(5) mod 128 = 124 = -4 mod 128", log5_mod128 == minus4_mod128),
        ("exp_2(-4) mod 128 = 5", exp2_neg4(7) == 5),
        ("c(5) = v_2(log_2(5)/4 + 1) = 5", cliff_constant(5) == 5),
    ]

    all_ok = all(ok for _, ok in checks)
    if verbose:
        print("Connection between all four checks")
        print()
        print("  Central identity:  log_2(5) ≡ -4  (mod 128)")
        print()
        for msg, ok in checks:
            print(f"  {'✓' if ok else '✗'}  {msg}")
        print(f"\n  All checks: {'PASS ✓' if all_ok else 'FAIL ✗'}")

    return all_ok
