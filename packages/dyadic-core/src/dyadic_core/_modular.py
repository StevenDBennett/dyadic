"""
_modular.py
-----------
Newton-lifted modular inverse and LTE-based dual-coordinate addition.
"""

from __future__ import annotations

from dyadic_core._exceptions import DomainError, NonInvertibleError
from dyadic_core._util import bitmask, valuation

__all__ = [
    "dual_add",
    "modinv_newton",
]


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


def dual_add(
    v_a: int,
    alpha_a: int,
    e_a: int,
    v_b: int,
    alpha_b: int,
    e_b: int,
    k: int,
) -> tuple[int | None, int, int]:
    """
    Exact addition in dual (v, alpha, e) coordinates via the Lifting The Exponent Lemma.

    Given a = 2^{v_a} . (-1)^{alpha_a} . 5^{e_a} and b = 2^{v_b} . (-1)^{alpha_b} . 5^{e_b},
    returns (v_sum, alpha_sum, e_sum) for a + b mod 2^k without converting
    back to the group representation.

    v_sum is None when exact cancellation occurs (a + b == 0 mod 2^k).

    Cases
    -----
    1. Different valuations: smaller v dominates (annihilation)
    2. Same sign, same e: exact doubling (v -> v+1)
    3. Same sign, different e: v2(sum) = v + 1 via 5^m == 1 mod 4
    4. Opposite signs: v2(sum) = v + 2 + v2(De) via LTE
    5. Exact cancellation: e_a = e_b and signs oppose -> v = None (infinite)
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
            return (None, 0, 0)
        v2 = valuation(e_diff)
        assert v2 is not None  # e_diff > 0, so valuation is finite
        v_lte = 2 + v2
        return (min(v + v_lte, k), 0, e_a if e_a < e_b else e_b)
