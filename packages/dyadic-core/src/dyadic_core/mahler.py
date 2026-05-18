"""
mahler.py
---------
Mahler basis operator calculus on the function space Z/N → Z/2^k.

Every integer-valued function f: Z → Z has a unique Mahler (binomial)
basis expansion:

    f(x) = Σ a_n · C(x, n)    where C(x, n) = x(x−1)…(x−n+1)/n!

The coefficients a_n encode the function in a form that diagonalises
the forward-difference operator.

Operators on the Mahler coefficient space:
    D (Dirac)  — the forward difference: D(e_n) = −e_{n−1}
    T (Volterra) — the right inverse: T(e_n) = −e_{n+1}

Key asymmetry (proved in tests):
    D∘T = id  on  ker(ε) = span{e_n : n ≥ 1}
    T∘D = id  on  span{e_n : n ≥ 2}   (NOT on all of ker(ε))

The asymmetry arises because D(e₁) = −e₀ exits ker(ε), so T∘D is
undefined at e₁.
"""
from __future__ import annotations


class MahlerCalculus:
    """
    Operations on the Mahler (binomial) basis.

    The Mahler polynomials C(x, n) = x(x−1)…(x−n+1)/n! form a basis
    for the ring of integer-valued functions.  Any such function f
    can be written uniquely as f(x) = Σ a_n · C(x, n).
    """

    @staticmethod
    def mahler_polynomial(n: int, x: int) -> int:
        """
        Compute the Mahler polynomial C(x, n) = binomial(x, n).

        For n < 0 returns 0; for n = 0 returns 1.
        """
        if n < 0:
            return 0
        if n == 0:
            return 1
        if x < 0:
            sign = -1 if n % 2 == 1 else 1
            return sign * MahlerCalculus.mahler_polynomial(n - x - 1, n)
        result = 1
        for i in range(n):
            result = result * (x - i) // (i + 1)
        return result

    @staticmethod
    def to_mahler(f: list[int], max_degree: int) -> list[int]:
        """
        Convert a function given by its values f[0], f[1], … to Mahler
        coefficients via the finite-difference table.

        The n-th coefficient is the n-th forward difference at 0:
            a_n = Δⁿf(0)
        """
        N = min(len(f) - 1, max_degree)
        coeffs: list[int] = []
        diff = f[:N + 1]
        for _ in range(N + 1):
            coeffs.append(diff[0])
            diff = [diff[i + 1] - diff[i] for i in range(len(diff) - 1)]
        return coeffs

    @staticmethod
    def from_mahler(coeffs: list[int], x: int) -> int:
        """
        Evaluate a function from its Mahler coefficients at point x.

        f(x) = Σ a_n · C(x, n)
        """
        result = 0
        for n, a in enumerate(coeffs):
            result += a * MahlerCalculus.mahler_polynomial(n, x)
        return result

    @staticmethod
    def dirac_operator(coeffs: list[int]) -> list[int]:
        """
        Dirac (forward-difference) operator on Mahler coefficients.

        D(a₀, a₁, a₂, …) = (−a₁, −a₂, −a₃, …)

        This lowers the degree by one: D(e_n) = −e_{n−1}.
        """
        if not coeffs:
            return []
        return [-coeffs[n] for n in range(1, len(coeffs))]

    @staticmethod
    def volterra_operator(coeffs: list[int]) -> list[int]:
        """
        Volterra (right-inverse of Dirac) operator on Mahler coefficients.

        Requires a₀ = 0 (the augmentation kernel, i.e. f ∈ ker(ε)).

        T(0, a₁, a₂, …) = (0, 0, −a₁, −a₂, …)

        This raises the degree by one: T(e_n) = −e_{n+1}.
        """
        if not coeffs:
            return []
        if coeffs[0] != 0:
            raise ValueError(
                "Volterra operator requires a₀ = 0 (augmentation kernel)"
            )
        return [0] + [-c for c in coeffs]

    @staticmethod
    def truncate(coeffs: list[int], k: int) -> list[int]:
        """Reduce coefficients modulo 2^k and truncate to length k."""
        mod = 1 << k
        return [c % mod for c in coeffs[:k]]
