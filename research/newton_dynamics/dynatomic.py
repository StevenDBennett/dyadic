"""
dynatomic.py
------------
Dynatomic polynomial construction for the Newton map.

The dynatomic polynomial counts periodic points of exact period n:

    Φ_n^*(x) = ∏_{d|n} (N^d(x) - x)^{μ(n/d)}

Using the recurrence form N^d(x) = A_d(x)/B_d(x), the factor
N^d(x) - x = (A_d(x) - x * B_d(x)) / B_d(x) is cleared to a
polynomial F_d(x) = A_d(x) - x * B_d(x).

The final result is expressed in u = x^3 via the symmetry N(ζx) = ζN(x).
"""
from __future__ import annotations


from .poly import poly_mul, poly_add, poly_scalar_mul, poly_divmod
from .iterates import mobius

__all__ = ["dynatomic_polynomial"]


def dynatomic_polynomial(
    period: int,
    iterates: list[tuple[List[int], list[int]]],
) -> list[int]:
    """
    Compute the dynatomic polynomial Φ_n^*(u) in u = x^3.

    Parameters
    ----------
    period : int
        The period n.
    iterates : list of (A_d, B_d)
        Output of compute_iterates(period) — must cover at least d = period.

    Returns
    -------
    list of int
        Coefficients of Φ_n^*(u) from constant term to leading term.
    """
    # Step 1: compute Φ_n^*(x) — the dynatomic polynomial in x
    phi_x = [1]
    # Process divisors in descending order so that when we divide
    # by a small-d F_d, the accumulator already contains the
    # corresponding factor from a higher iterate.
    for d in reversed(range(1, period + 1)):
        if period % d == 0:
            mu = mobius(period // d)
            A_d, B_d = iterates[d]
            # F_d = A_d - x*B_d
            xB = [0] + B_d
            F_d = poly_add(A_d, poly_scalar_mul(xB, -1))
            if mu == 1:
                phi_x = poly_mul(phi_x, F_d)
            elif mu == -1:
                quot, rem = poly_divmod(phi_x, F_d)
                if all(r == 0 for r in rem):
                    phi_x = quot
                else:
                    raise ValueError(
                        f"Division failed for d={d}: "
                        f"remainder has non-zero entries"
                    )
    # Step 2: compress x-polynomial to u = x^3 — only keep coeffs at
    # multiples of 3 (the symmetry N(ζx) = ζN(x) guarantees all other
    # coefficients are zero).
    return phi_x[0::3]
