"""
_series.py
----------
2-adic exponential and logarithm via exact integer series.
"""

from __future__ import annotations

from dyadic_core._exceptions import DomainError
from dyadic_core._util import valuation

__all__ = [
    "g0",
    "padic_exp",
    "padic_log",
]


def _series_accumulate(
    x_exact: int,
    vx: int,
    k: int,
    *,
    use_factorial: bool,
    start: int,
    sign_alternate: bool,
) -> int:
    """
    Shared series accumulator for 2-adic exp and log.

    Iterates n = 1.., computing v2(term) = n*vx - v2(denominator).
    When v2(term) >= k the series is truncated.

    All intermediate values are kept modulo 2^k to avoid
    super-linear integer growth.

    Parameters
    ----------
    x_exact : int
        The exact integer value of x (may be negative).
    vx : int
        2-adic valuation of x (guaranteed >= 2).
    k : int
        Bit precision.
    use_factorial : bool
        If True, denominator is n! (padic_exp); if False, denominator is n (padic_log).
    start : int
        Initial value of result (1 for exp, 0 for log).
    sign_alternate : bool
        If True, apply (-1)^{n+1} sign factor (padic_log).
    """
    mod = 1 << k
    result = start

    # Odd part of x_exact: x_exact = 2^vx * odd_x
    odd_x = x_exact >> vx
    odd_x_abs = abs(odd_x)
    odd_x_mod = odd_x_abs % mod
    # power_odd_abs tracks |odd_x|^n mod 2^k
    power_odd_abs = odd_x_mod

    # For factorial mode: fact_odd tracks the odd part of n! mod 2^k
    fact_odd = 1
    v2_fact = 0  # v2(n!) accumulated incrementally

    for n in range(1, k * 2):
        vn = valuation(n)
        assert vn is not None  # n >= 1
        if use_factorial:
            fact_odd = (fact_odd * (n >> vn)) % mod
            v2_fact += vn
            den_odd = fact_odd
            v_term = n * vx - v2_fact
        else:
            den_odd = n >> vn
            v_term = n * vx - vn
        if v_term >= k:
            break

        # num_abs = |odd_x|^n mod 2^k
        num_abs = power_odd_abs % mod
        term = pow(den_odd, -1, mod) * num_abs % mod
        # Sign of odd_x^n: negative when odd_x < 0 and n is odd
        if odd_x < 0 and n % 2 == 1:
            term = (-term) % mod
        if sign_alternate and n % 2 == 0:
            term = (-term) % mod
        term = term * pow(2, v_term, mod) % mod
        result = (result + term) % mod

        # Advance power_odd_abs for next iteration
        power_odd_abs = (power_odd_abs * odd_x_mod) % mod
    return result


def padic_exp(x: int, k: int) -> int:
    """
    Compute exp(x) mod 2^k via exact integer arithmetic.

    Requires v2(x) >= 2 (the 2-adic exponential converges on 4Z2).
    """
    if x == 0:
        return 1
    vx = valuation(abs(x))
    if vx is None or vx < 2:
        raise DomainError(f"exp requires v2(x) >= 2, got v2={vx}")
    return _series_accumulate(x, vx, k, use_factorial=True, start=1, sign_alternate=False)


def padic_log(g: int, k: int) -> int:
    """
    Compute log(g) mod 2^k via exact integer arithmetic.

    Requires g == 1 mod 4 (i.e. v2(g-1) >= 2, so g is in the domain of
    the 2-adic logarithm).  Uses the series:
        log(g) = (g-1) - (g-1)^2/2 + (g-1)^3/3 - ...
    """
    mod = 1 << k
    g_mod = g % mod
    x_exact = g_mod - 1
    if x_exact > (mod >> 1):
        x_exact -= mod
    if x_exact == 0:
        return 0
    vx = valuation(abs(x_exact))
    if vx is None or vx < 2:
        raise DomainError(f"log requires v2(g-1) >= 2 (g == 1 mod 4), got v2={vx}")
    return _series_accumulate(x_exact, vx, k, use_factorial=False, start=0, sign_alternate=True)


def g0(k: int) -> int:
    """
    The cliff centre: g0 = exp2(-4) mod 2^k.

    This is the unique 2-adic unit where log(g0) = -4.
    The hardware approximation -123 agrees with g0 to 13 bits:
    v2(g0 - (-123)) = 13.
    """
    return padic_exp(-4, k)
