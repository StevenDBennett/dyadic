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
    power_exact = x_exact
    factorial_exact = 1
    v2_denom = 0
    for n in range(1, k * 2):
        vn = valuation(n)
        assert vn is not None  # n >= 1
        if use_factorial:
            factorial_exact *= n
            v2_denom += vn
            den_odd = factorial_exact >> v2_denom
        else:
            v2_denom = vn
            den_odd = n >> vn
        v_term = n * vx - v2_denom
        if v_term >= k:
            break
        num_stripped = (
            power_exact >> (n * vx) if power_exact >= 0 else -((-power_exact) >> (n * vx))
        )
        term = pow(int(den_odd), -1, mod) * int(abs(num_stripped)) % mod
        if num_stripped < 0:
            term = (-term) % mod
        if sign_alternate and n % 2 == 0:
            term = (-term) % mod
        term = term * pow(2, v_term, mod) % mod
        result = (result + term) % mod
        if use_factorial:
            power_exact *= x_exact
        else:
            power_exact *= x_exact
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
