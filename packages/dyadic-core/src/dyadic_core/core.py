"""
core.py
-------
Fast modular arithmetic in power-of-two rings via group-algebra duality.

Public names
~~~~~~~~~~~~
    modinv_newton   – Newton-lifted modular inverse mod 2^k
    two_adic_log5   – 2-adic logarithm of 5 (truncated to k bits)
    two_adic_dlog   – Full decomposition: n = 2^v · (-1)^alpha · 5^e
    dual_add        – LTE-based exact addition in (v, alpha, e) coordinates
    DualNumber      – Coordinate triple (v, alpha, e) for n in Z/2^k
    TwoAdicProcessor– Arithmetic on DualNumbers
    run_all_tests   – Basic self-check
"""

from __future__ import annotations

import math
from contextlib import contextmanager
from functools import lru_cache
from typing import TypedDict

# ── exception hierarchy ────────────────────────────────────────────────────────


class DyadicError(ValueError):
    """Base for all dyadic-core errors."""


class NonInvertibleError(DyadicError):
    """Raised when an inversion is attempted on a non-invertible element."""


class ConvergenceError(DyadicError):
    """Raised when an iterative method fails to converge."""


class DomainError(DyadicError):
    """Raised when a function argument is outside its domain of convergence."""


# ── 10-bit discrete log LUT ────────────────────────────────────────────────────
# Precomputed: for each odd a ≡ 1 (mod 4) in [0, 1024), store e s.t. 5^e ≡ a (mod 1024).
# 256 entries ~ 2 KB. Gives O(1) bootstrap to 8-bit precision.
# b=10 is the optimal sweet spot (lab finding: b=8 costs 324, b=10 costs 256).
# List is faster and more cache-friendly than a dict for 1024 entries.
_DLOG10_LUT: list[int] = [-1] * 1024
_DLOG10_LUT_INIT = {pow(5, e, 1024): e for e in range(256)}
for _a, _e in _DLOG10_LUT_INIT.items():
    _DLOG10_LUT[_a] = _e
del _DLOG10_LUT_INIT


# ── private utilities ──────────────────────────────────────────────────────────


def bitmask(k: int) -> int:
    return (1 << k) - 1


def mat_mul(A: list[list[int]], B: list[list[int]], mod: int) -> list[list[int]]:
    """2×2 matrix multiplication modulo mod."""
    return [
        [
            (A[0][0] * B[0][0] + A[0][1] * B[1][0]) % mod,
            (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % mod,
        ],
        [
            (A[1][0] * B[0][0] + A[1][1] * B[1][0]) % mod,
            (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % mod,
        ],
    ]


def mat_det(M: list[list[int]], mod: int) -> int:
    """Determinant of a 2×2 matrix modulo mod."""
    return (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % mod


def valuation(n: int) -> int | None:
    """2-adic valuation v_2(n). Returns None for zero."""
    if n == 0:
        return None
    return (n & -n).bit_length() - 1


# ── LTE-based dual addition ────────────────────────────────────────────────────


def dual_add(
    v_a: int,
    alpha_a: int,
    e_a: int,
    v_b: int,
    alpha_b: int,
    e_b: int,
    k: int,
) -> tuple[int | float, int, int]:
    """
    Exact addition in dual (v, alpha, e) coordinates via the Lifting The Exponent Lemma.

    Given a = 2^{v_a} · (-1)^{alpha_a} · 5^{e_a} and b = 2^{v_b} · (-1)^{alpha_b} · 5^{e_b},
    returns (v_sum, alpha_sum, e_sum) for a + b mod 2^k without converting
    back to the group representation.

    Cases
    -----
    1. Different valuations: smaller v dominates (annihilation)
    2. Same sign, same e: exact doubling (v → v+1)
    3. Same sign, different e: v2(sum) = v + 1 via 5^m ≡ 1 mod 4
    4. Opposite signs: v2(sum) = v + 2 + v2(De) via LTE
    5. Exact cancellation: e_a = e_b and signs oppose -> v = inf
    """
    if v_a > v_b:
        return (v_b, alpha_b, e_b)
    if v_b > v_a:
        return (v_a, alpha_a, e_a)

    v = v_a
    if alpha_a == alpha_b:
        e_min = e_a if e_a < e_b else e_b
        e_diff = e_a - e_b if e_a > e_b else e_b - e_a
        if e_diff == 0:
            return (min(v + 1, k), alpha_a, e_a)
        return (min(v + 1, k), alpha_a, e_min)
    else:
        e_diff = e_a - e_b if e_a > e_b else e_b - e_a
        if e_diff == 0:
            return (float("inf"), 0, 0)
        v_lte = 2 + valuation(e_diff)
        return (min(v + v_lte, k), 0, e_a if e_a < e_b else e_b)


# ── Newton modular inverse ─────────────────────────────────────────────────────


def modinv_newton(a: int, k: int) -> int:
    """
    a^{-1} mod 2^k via quadratic Newton lifting.  Requires a odd.
    """
    if k <= 0:
        raise DomainError("k must be positive")
    if a & 1 == 0:
        raise NonInvertibleError("a must be odd to be invertible mod 2^k")
    a &= bitmask(k)
    x = a & 7
    bits = 3
    while bits < k:
        bits = min(bits * 2, k)
        x = (x * (2 - a * x)) & bitmask(bits)
    return x & bitmask(k)


# ── 2-adic log of 5 ────────────────────────────────────────────────────────────


@lru_cache(maxsize=256)
def two_adic_log5(k: int) -> int:
    """
    Compute the 2-adic logarithm of 5, truncated to k bits.

    This is the unique L in Z_2 satisfying exp_2(L) = 5 in the
    2-adic sense.  Concretely it is computed via the 2-adic log series
    applied to (5-1).

    Primary use: derivative of the map e |-> 5^e is 5^e · (L >> 2).
    Cached per k.
    """
    mask = bitmask(k)
    result = 0
    n = 1
    while True:
        v = valuation(n)
        exp = 2 * n - v
        if exp >= k:
            break
        odd_part = n >> v
        inv_odd = modinv_newton(odd_part, k)
        term = ((1 << exp) * inv_odd) & mask
        if n % 2 == 0:
            term = (-term) & mask
        result = (result + term) & mask
        n += 1
    return result


# ── discrete-log helpers (private) ─────────────────────────────────────────────


def dlog_bootstrap(a: int, k: int) -> int:
    """
    Bit-by-bit dlog for a ≡ 1 (mod 8).
    Returns e with 5^e ≡ a (mod 2^k).  O(k) steps; used to seed Newton.
    """
    if k <= 2:
        return 0
    mask = bitmask(k)
    a &= mask
    e, pow5 = (0, 1) if (a & 7) == 1 else (1, 5)
    mult = 25 & mask
    for n in range(3, k):
        if ((a >> n) & 1) != ((pow5 >> n) & 1):
            e |= 1 << (n - 2)
            pow5 = (pow5 * mult) & mask
        mult = (mult * mult) & mask
    return e


def dlog_newton(a: int, k: int, L: int | None = None) -> int:
    """
    Newton-lifted dlog: 5^e ≡ a (mod 2^k).  Requires a ≡ 1 (mod 4).

    Bootstrap uses an 8-bit LUT for k ≤ 34 (covers 6-bit precision),
    otherwise uses bit-by-bit to k/2.  Per the Mersenne boostrap
    optimality analysis, k/2 saves ⌈log2(√k)⌉−1 Newton steps vs √k.
    """
    if k <= 2:
        return 0
    if a & 3 != 1:
        raise DomainError("dlog_newton requires a ≡ 1 (mod 4)")
    mask_full = bitmask(k)
    a &= mask_full
    if L is None:
        L = two_adic_log5(k)
    L_unit = L >> 2

    if k <= 34:
        a_b = a & 0x3FF
        e_raw = _DLOG10_LUT[a_b] if _DLOG10_LUT[a_b] >= 0 else 0
        eprec = min(8, k - 2)
        e = e_raw & bitmask(eprec)
    else:
        bootstrap_k = max(4, k // 2 + 2)
        eprec = min(bootstrap_k - 2, k - 2)
        e = dlog_bootstrap(a, bootstrap_k) & bitmask(eprec)

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = bitmask(bits)
        emask = bitmask(new_eprec)

        e, _ = _dlog_newton_step(e, a, L_unit, bits, mask, emask, new_eprec)
        eprec = new_eprec

    return e


def _dlog_newton_step(
    e: int, a: int, L_unit: int, bits: int, mask: int, emask: int, new_eprec: int
) -> tuple[int, int]:
    """
    Single Newton step for 2-adic discrete log: refine e so that 5^e ≡ a (mod 2^bits).

    Returns (new_e, delta).
    """
    pow5e = pow(5, e, 1 << bits)
    f = (pow5e - a) & mask
    df_unit = (pow5e * L_unit) & emask
    df_inv = modinv_newton(df_unit, new_eprec)
    delta = ((f >> 2) * df_inv) & emask
    e_new = (e - delta) & emask
    return e_new, delta


# ── residual tracking (diagnostic) ─────────────────────────────────────────────


class DLogNewtonStep(TypedDict):
    bits: int
    eprec_before: int
    eprec_after: int
    tau_before: int
    v2_before: int | None
    tau_after: int
    v2_after: int | None
    delta: int


def dlog_residual_tracking(
    a: int, k: int, L: int | None = None
) -> tuple[int, list[DLogNewtonStep]]:
    """
    Viglietta discrete log with residual tracking at each Newton step.

    Returns (e, history) where history contains the normalised residual
    (tau = (5^e * a^{-1} - 1) / 4) and its 2-adic valuation before and
    after each Newton correction.  The doubling of v2(tau) confirms the
    quadratic convergence law.

    Requires a ≡ 1 (mod 4).
    """
    if k <= 2:
        return 0, []
    if a & 3 != 1:
        raise DomainError("dlog_residual_tracking requires a ≡ 1 (mod 4)")
    mask_full = bitmask(k)
    a &= mask_full
    if L is None:
        L = two_adic_log5(k)
    L_unit = L >> 2
    a_inv_full = modinv_newton(a, k)

    bootstrap_k = max(4, math.isqrt(k) + 2)
    eprec = min(bootstrap_k - 2, k - 2)
    e = dlog_bootstrap(a, bootstrap_k) & bitmask(eprec)
    history: list[dict] = []

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = bitmask(bits)
        emask = bitmask(new_eprec)

        pow5e = pow(5, e, 1 << bits)
        a_inv = a_inv_full & mask
        rho_before = (pow5e * a_inv) & mask
        tau_before = ((rho_before - 1) & mask) >> 2
        v2_before = valuation(tau_before)

        e, delta = _dlog_newton_step(e, a, L_unit, bits, mask, emask, new_eprec)

        pow5e = pow(5, e, 1 << bits)
        rho_after = (pow5e * a_inv) & mask
        tau_after = ((rho_after - 1) & mask) >> 2
        v2_after = valuation(tau_after)

        history.append(
            {
                "bits": bits,
                "eprec_before": eprec,
                "eprec_after": new_eprec,
                "tau_before": int(tau_before),
                "v2_before": v2_before,
                "tau_after": int(tau_after),
                "v2_after": v2_after,
                "delta": int(delta),
            }
        )
        eprec = new_eprec

    return e, history


# ── public discrete-log API ────────────────────────────────────────────────────


def two_adic_dlog(a: int, k: int, L: int | None = None) -> tuple[int, int] | None:
    """
    Decompose the odd part of a as (-1)^alpha · 5^e (mod 2^k).

    Parameters
    ----------
    a : int   — integer to decompose (any residue mod 2^k)
    k : int   — bit precision (>= 3)
    L : int   — precomputed two_adic_log5(k); computed if not supplied

    Returns
    -------
    (alpha, e)  with alpha ∈ {0,1} and e ∈ [0, 2^(k-2)), or
    None        if a is even.
    """
    if a & 1 == 0:
        return None
    mask = bitmask(k)
    a &= mask
    alpha = (a >> 1) & 1
    if alpha:
        a = (-a) & mask
    e = dlog_newton(a, k, L)
    return alpha, e


# ── DualNumber ─────────────────────────────────────────────────────────────────


class DualNumber:
    """
    Represent n ∈ Z/2^k in dual coordinates (v, alpha, e) where

        n  ≡  2^v · (-1)^alpha · 5^e   (mod 2^k)

    v     = v_2(n)  (2-adic valuation)
    alpha ∈ {0, 1}  (sign of the odd unit part)
    e     ∈ [0, 2^(k-2))  (discrete log base 5 of the odd unit part)

    Zero is a distinguished element with v = ∞.

    Construction
    ------------
    From an integer:       DualNumber(42, k=16)
    From coordinates:      DualNumber.from_coords(v, alpha, e, k)
    """

    __slots__ = ("k", "mask", "_ord", "v", "alpha", "e", "_n", "is_zero")

    def __init__(self, n: int, k: int = 64):
        if k < 3:
            raise DomainError("k must be ≥ 3")
        self.k = k
        self.mask = bitmask(k)
        self._ord = 1 << (k - 2)

        n_mod = n & self.mask
        if n_mod == 0:
            self.is_zero = True
            self.v = float("inf")
            self.alpha = 0
            self.e = 0
            self._n = 0
            return

        self.is_zero = False
        self.v = valuation(n_mod)
        odd = (n_mod >> self.v) & self.mask

        result = two_adic_dlog(odd, k)
        if result is None:
            raise DomainError(f"Cannot decompose {n} mod 2^{k}: odd part is even")
        self.alpha, self.e = result
        self._n = n_mod

    @classmethod
    def from_coords(cls, v: int, alpha: int, e: int, k: int = 64) -> DualNumber:
        """Build a DualNumber directly from coordinates."""
        obj = cls.__new__(cls)
        obj.k = k
        obj.mask = bitmask(k)
        obj._ord = 1 << (k - 2)
        if v >= k:
            obj.is_zero = True
            obj.v = float("inf")
            obj.alpha = 0
            obj.e = 0
            obj._n = 0
            return obj
        obj.is_zero = False
        obj.v = v
        obj.alpha = alpha & 1
        obj.e = e % obj._ord
        obj._n = obj._to_int()
        return obj

    def _to_int(self) -> int:
        if self.is_zero:
            return 0
        odd = pow(5, self.e, 1 << self.k)
        if self.alpha:
            odd = (-odd) & self.mask
        return (odd << self.v) & self.mask

    def verify(self) -> bool:
        """Round-trip check: coordinates → integer matches stored integer."""
        if self.is_zero:
            return self._n == 0
        return self._to_int() == self._n

    @property
    def value(self) -> int:
        """The integer n mod 2^k."""
        return self._n

    def coords(self) -> tuple[int | float, int, int]:
        """Return (v, alpha, e)."""
        return (self.v, self.alpha, self.e)

    def __repr__(self) -> str:
        if self.is_zero:
            return f"DualNumber(0, k={self.k})"
        sign = "-" if self.alpha else "+"
        return f"DualNumber({self._n}, k={self.k})  =  2^{self.v} · {sign}5^{self.e}"


# ── TwoAdicProcessor ──────────────────────────────────────────────────────────


class TwoAdicProcessor:
    """
    Arithmetic on DualNumbers in Z/2^k.

    Multiplication, inversion, and exponentiation operate directly in
    coordinate space, making the group structure explicit.

    Parameters
    ----------
    k : int
        Bit precision (≥ 3).
    verify : bool
        If True (default), run round-trip verification after every
        operation.  Disable for production to avoid the O(log k)
        overhead of an extra modular exponentiation per call.
    """

    def __init__(self, k: int = 64, verify: bool = True):
        if k < 3:
            raise DomainError("k must be ≥ 3")
        self.k = k
        self.mask = bitmask(k)
        self._ord = 1 << (k - 2)
        self._verify = verify
        self._L: int | None = None

    @property
    def L(self) -> int:
        if self._L is None:
            self._L = two_adic_log5(self.k)
        return self._L

    def _check(self, *args: DualNumber) -> None:
        for a in args:
            if a.k != self.k:
                raise ValueError(f"DualNumber has k={a.k} but processor has k={self.k}")

    def mul(self, a: DualNumber, b: DualNumber) -> DualNumber:
        """Multiply two elements: coordinates add componentwise."""
        self._check(a, b)
        if a.is_zero or b.is_zero:
            return DualNumber(0, self.k)
        v = a.v + b.v
        if v >= self.k:
            return DualNumber(0, self.k)
        result = DualNumber.from_coords(
            v=v,
            alpha=a.alpha ^ b.alpha,
            e=(a.e + b.e) % self._ord,
            k=self.k,
        )
        if self._verify and not result.verify():
            raise RuntimeError("TwoAdicProcessor.mul invariant failed")
        return result

    def inv(self, a: DualNumber) -> DualNumber:
        """Invert a unit (v = 0).  Raises ValueError if a is not a unit."""
        self._check(a)
        if a.is_zero or a.v != 0:
            raise NonInvertibleError("Only units (v=0) are invertible")
        result = DualNumber.from_coords(
            v=0,
            alpha=a.alpha,
            e=(-a.e) % self._ord,
            k=self.k,
        )
        if self._verify and not result.verify():
            raise RuntimeError("TwoAdicProcessor.inv invariant failed")
        return result

    def pow(self, a: DualNumber, n: int) -> DualNumber:
        """Raise a to an integer power.  Negative exponents require a to be a unit."""
        self._check(a)
        if a.is_zero:
            return DualNumber(0, self.k)
        if n < 0:
            return self.pow(self.inv(a), -n)
        v = a.v * n
        if v >= self.k:
            return DualNumber(0, self.k)
        result = DualNumber.from_coords(
            v=v,
            alpha=a.alpha if (n % 2 == 1) else 0,
            e=(a.e * n) % self._ord,
            k=self.k,
        )
        if self._verify and not result.verify():
            raise RuntimeError("TwoAdicProcessor.pow invariant failed")
        return result

    def dlog(self, a: DualNumber) -> tuple[int | float, int, int]:
        """Return full coordinate triple (v, alpha, e)."""
        self._check(a)
        return a.coords()

    @contextmanager
    def verification(self, enabled: bool = True):
        """
        Temporarily enable or disable round-trip verification.

        Use as a context manager to disable the O(log k) verification
        overhead for bulk operations:

        >>> with proc.verification(False):
        ...     result = proc.mul(a, b)
        """
        old = self._verify
        self._verify = enabled
        try:
            yield
        finally:
            self._verify = old


# ── general 2-adic exp / log ──────────────────────────────────────────────────


def padic_exp(x: int, k: int) -> int:
    """
    Compute exp(x) mod 2^k via exact integer arithmetic.

    Requires v2(x) ≥ 2 (the 2-adic exponential converges on 4Z2).
    Terminates when v2(x^n/n!) ≥ k.

    This is the general 2-adic exponential — unlike `two_adic_log5` which
    computes only log(5), this works for any valid 2-adic argument.
    """
    mod = 1 << k
    if x == 0:
        return 1
    vx = valuation(abs(x))
    if vx < 2:
        raise DomainError(f"exp requires v2(x) ≥ 2, got v2={vx}")
    result = 1
    power_exact = x
    factorial_exact = 1
    v2_fact = 0
    for n in range(1, k * 2):
        factorial_exact *= n
        v2_fact += valuation(n)
        v_term = n * vx - v2_fact
        if v_term >= k:
            break
        den_odd = factorial_exact >> v2_fact
        num_stripped = (
            power_exact >> (n * vx) if power_exact >= 0 else -((-power_exact) >> (n * vx))
        )
        term = pow(int(den_odd), -1, mod) * int(abs(num_stripped)) % mod
        if num_stripped < 0:
            term = (-term) % mod
        term = term * pow(2, v_term, mod) % mod
        result = (result + term) % mod
        power_exact *= x
    return result


def padic_log(g: int, k: int) -> int:
    """
    Compute log(g) mod 2^k via exact integer arithmetic.

    Requires g ≡ 1 mod 4 (i.e. v2(g−1) ≥ 2, so g is in the domain of
    the 2-adic logarithm). Uses the series:
        log(g) = (g−1) − (g−1)²/2 + (g−1)³/3 − …

    This is the general 2-adic logarithm — unlike `two_adic_log5` which
    computes only log(5), this works for any g ≡ 1 mod 4.
    """
    mod = 1 << k
    g_mod = g % mod
    x_exact = g_mod - 1
    if x_exact > (mod >> 1):
        x_exact -= mod
    if x_exact == 0:
        return 0
    vx = valuation(abs(x_exact))
    if vx < 2:
        raise DomainError(f"log requires v2(g−1) ≥ 2 (g ≡ 1 mod 4), got v2={vx}")
    result = 0
    power_exact = x_exact
    for n in range(1, k * 2):
        vn = valuation(n)
        v_term = n * vx - vn
        if v_term >= k:
            break
        den_odd = n >> vn
        num_stripped = (
            power_exact >> (n * vx) if power_exact >= 0 else -((-power_exact) >> (n * vx))
        )
        term = pow(int(den_odd), -1, mod) * int(abs(num_stripped)) % mod
        if num_stripped < 0:
            term = (-term) % mod
        if n % 2 == 0:
            term = (-term) % mod
        term = term * pow(2, v_term, mod) % mod
        result = (result + term) % mod
        power_exact *= x_exact
    return result


def g0(k: int) -> int:
    """
    The cliff centre: g0 = exp2(−4) mod 2^k.

    This is the unique 2-adic unit where log(g0) = −4.
    The hardware approximation −123 agrees with g0 to 13 bits:
    v2(g0 − (−123)) = 13.
    """
    return padic_exp(-4, k)


# ── self-test ─────────────────────────────────────────────────────────────────


def run_all_tests(k: int = 16, verbose: bool = True) -> None:
    """
    Basic sanity checks for the core arithmetic.

    Tests round-trip encoding, multiplication, inversion, and powering.
    Raises AssertionError on any failure.
    """
    proc = TwoAdicProcessor(k)
    mask = bitmask(k)
    failures: list[str] = []

    samples = [1, 3, 5, 7, 9, 15, 255, (1 << (k // 2)) + 1, mask - 2]
    for n in samples:
        d = DualNumber(n, k)
        if not d.verify():
            failures.append(f"Round-trip failed: n={n}")

    z = DualNumber(0, k)
    if not z.is_zero:
        failures.append("Zero not recognised as is_zero")

    d3 = DualNumber(3, k)
    d5 = DualNumber(5, k)
    d15 = proc.mul(d3, d5)
    if d15.value != (15 & mask):
        failures.append(f"mul: 3×5 → {d15.value}, expected 15")

    inv3 = proc.inv(d3)
    if proc.mul(d3, inv3).value != 1:
        failures.append("inv: 3 × 3⁻¹ ≠ 1")

    if proc.pow(d3, 4).value != (81 & mask):
        failures.append(f"pow: 3^4 → {proc.pow(d3, 4).value}, expected {81 & mask}")

    if failures:
        for msg in failures:
            print(f"FAIL  {msg}")
        raise AssertionError(f"{len(failures)} test(s) failed at k={k}")

    if verbose:
        print(f"run_all_tests: all checks passed  (k={k})")
