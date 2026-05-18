"""
iterates.py
-----------
Newton map iterates for N(x) = (2x^3+1)/(3x^2).

The d-th iterate is a rational function N^d(x) = A_d(x) / B_d(x)
with numerators and denominators defined by the recurrence:

    A_{d+1} = 2*A_d^3 + B_d^3
    B_{d+1} = 3*A_d^2 * B_d
    A_0 = x,  B_0 = 1
"""
from __future__ import annotations

from math import gcd

from .poly import poly_mul, poly_add, poly_scalar_mul, poly_pow

__all__ = ["mobius", "compute_iterates"]


def mobius(n: int) -> int:
    """Möbius function μ(n)."""
    if n == 1:
        return 1
    factors = {}
    temp = n
    d = 2
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    if any(v > 1 for v in factors.values()):
        return 0
    return -1 if len(factors) % 2 else 1


def compute_iterates(
    max_k: int,
) -> list[tuple[List[int], list[int]]]:
    """
    Compute (A_d, B_d) for d = 0, 1, ..., max_k.

    Returns a list of (A_d, B_d) pairs where each is a coefficient list.
    """
    A = [0, 1]  # x
    B = [1]     # 1
    iterates = [(A, B)]
    for _ in range(max_k):
        A_cubed = poly_add(poly_scalar_mul(poly_pow(A, 3), 2), poly_pow(B, 3))
        A_sq = poly_pow(A, 2)
        A_sq_B = poly_mul(A_sq, B)
        B_new = poly_scalar_mul(A_sq_B, 3)
        A, B = A_cubed, B_new
        iterates.append((A, B))
    return iterates
