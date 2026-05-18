"""
_modular.py
-----------
Newton-lifted modular inverse and dual-coordinate addition.
"""

from __future__ import annotations

from dyadic_core._exceptions import DomainError, NonInvertibleError
from dyadic_core._util import bitmask

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
    Add two numbers given in dual (v, alpha, e) coordinates.

    Given a = 2^{v_a} . (-1)^{alpha_a} . 5^{e_a} and
          b = 2^{v_b} . (-1)^{alpha_b} . 5^{e_b},
    returns (v_sum, alpha_sum, e_sum) for a + b mod 2^k.

    Internally reconstructs integers, adds, and re-decomposes.
    LTE formulas for v_2(sum) hold in the equal-valuation cases
    but the full (alpha, e) coordinates always require re-decomposition.

    v_sum is None when exact cancellation occurs.
    """
    from dyadic_core._dual import DualNumber

    mask = bitmask(k)
    a = pow(5, e_a, 1 << k)
    if alpha_a:
        a = (-a) & mask
    a = (a << v_a) & mask

    b = pow(5, e_b, 1 << k)
    if alpha_b:
        b = (-b) & mask
    b = (b << v_b) & mask

    s = (a + b) & mask
    if s == 0:
        return (None, 0, 0)
    d = DualNumber(s, k)
    return (d.v, d.alpha, d.e)
