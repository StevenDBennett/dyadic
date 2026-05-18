"""
poly.py
-------
Polynomial arithmetic over ℤ (coefficient-list representation).

A polynomial is represented as a list of coefficients [a_0, a_1, ..., a_d]
where a_i is the coefficient of x^i and a_d ≠ 0 (unless the polynomial
is zero, represented as [0]).
"""
from __future__ import annotations


__all__ = [
    "poly_mul", "poly_add", "poly_scalar_mul",
    "poly_pow", "poly_divmod",
]


def poly_mul(p: list[int], q: list[int]) -> list[int]:
    """Multiply two polynomials (convolution)."""
    result = [0] * (len(p) + len(q) - 1)
    for i in range(len(p)):
        for j in range(len(q)):
            result[i + j] += p[i] * q[j]
    return result


def poly_add(p: list[int], q: list[int]) -> list[int]:
    """Add two polynomials."""
    n = max(len(p), len(q))
    result = [0] * n
    for i in range(len(p)):
        result[i] += p[i]
    for i in range(len(q)):
        result[i] += q[i]
    return result


def poly_scalar_mul(p: list[int], c: int) -> list[int]:
    """Multiply a polynomial by a scalar integer."""
    return [a * c for a in p]


def poly_pow(p: list[int], n: int) -> list[int]:
    """Raise a polynomial to a non-negative integer power."""
    result = [1]
    for _ in range(n):
        result = poly_mul(result, p)
    return result


def poly_divmod(dividend: list[int], divisor: list[int]) -> tuple[list[int], list[int]]:
    """
    Polynomial division with integer coefficients.

    Returns (quotient, remainder) where dividend = quotient * divisor + remainder
    and deg(remainder) < deg(divisor).  Raises ValueError if divisor is zero.
    """
    while len(dividend) > 1 and dividend[-1] == 0:
        dividend = dividend[:-1]
    while len(divisor) > 1 and divisor[-1] == 0:
        divisor = divisor[:-1]
    if len(divisor) == 1 and divisor[0] == 0:
        raise ValueError("Division by zero")
    if len(dividend) < len(divisor):
        return [0], dividend
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    for i in range(len(dividend) - len(divisor), -1, -1):
        if remainder[i + len(divisor) - 1] != 0:
            coeff = remainder[i + len(divisor) - 1] // divisor[-1]
            if remainder[i + len(divisor) - 1] % divisor[-1] == 0:
                quotient[i] = coeff
                for j in range(len(divisor)):
                    remainder[i + j] -= coeff * divisor[j]
    while len(remainder) > 1 and remainder[-1] == 0:
        remainder = remainder[:-1]
    return quotient, remainder
