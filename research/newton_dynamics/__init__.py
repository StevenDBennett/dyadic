"""
newton_dynamics
===============
p-adic Newton dynamics for the rational map N(x) = (2x^3+1)/(3x^2).

Provides polynomial arithmetic, Newton iterate computation, dynatomic
polynomial construction, and clean-prime analysis.  The main object
of study is the Newton map for cube roots:

    N(x) = (2x^3 + 1)/(3x^2)

with iterates N^d(x) = A_d(x)/B_d(x) satisfying

    A_{d+1} = 2*A_d^3 + B_d^3
    B_{d+1} = 3*A_d^2 * B_d

The dynatomic polynomial Φ_n^*(u) (u = x^3) counts periodic points
of period exactly n, and its special values govern multiplier products.
"""

from .poly import poly_mul, poly_add, poly_scalar_mul, poly_pow, poly_divmod
from .iterates import mobius, compute_iterates
from .dynatomic import dynatomic_polynomial
from .clean_primes import is_cube, tonelli_shanks, check_quadratic_cube_roots
from .data import (
    COEFFS_PERIOD4,
    COEFFS_PERIOD5,
    MULTIPLIERS_PERIOD4,
    MULTIPLIERS_PERIOD5,
    load_period6_coefficients,
    PERIOD6_PREDICTED,
)

__all__ = [
    "poly_mul", "poly_add", "poly_scalar_mul", "poly_pow", "poly_divmod",
    "mobius", "compute_iterates",
    "dynatomic_polynomial",
    "is_cube", "tonelli_shanks", "check_quadratic_cube_roots",
    "COEFFS_PERIOD4", "COEFFS_PERIOD5",
    "MULTIPLIERS_PERIOD4", "MULTIPLIERS_PERIOD5",
    "load_period6_coefficients",
    "PERIOD6_PREDICTED",
]
