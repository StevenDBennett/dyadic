"""
_dual.py
--------
DualNumber (v, alpha, e) coordinate system for Z/2^k.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager

from dyadic_core._dlog import two_adic_dlog, two_adic_log5
from dyadic_core._exceptions import DomainError, NonInvertibleError
from dyadic_core._util import bitmask, valuation

__all__ = [
    "DualNumber",
    "TwoAdicProcessor",
]


class DualNumber:
    """
    Represent n in Z/2^k in dual coordinates (v, alpha, e) where

        n  ==  2^v . (-1)^alpha . 5^e   (mod 2^k)

    v     = v_2(n)  (2-adic valuation)
    alpha in {0, 1}  (sign of the odd unit part)
    e     in [0, 2^(k-2))   (discrete log base 5 of the odd unit part)

    Zero is a distinguished element with v = None.
    """

    __slots__ = ("k", "mask", "_ord", "v", "alpha", "e", "_n", "is_zero")

    def __init__(self, n: int, k: int = 64):
        if k < 3:
            raise DomainError("k must be >= 3")
        self.k = k
        self.mask = bitmask(k)
        self._ord = 1 << (k - 2)

        n_mod = n & self.mask
        if n_mod == 0:
            self.is_zero = True
            self.v = None
            self.alpha = 0
            self.e = 0
            self._n = 0
            return

        self.is_zero = False
        self.v = valuation(n_mod)
        assert self.v is not None  # n_mod > 0, so valuation is finite
        odd = (n_mod >> self.v) & self.mask

        result = two_adic_dlog(odd, k)
        if result is None:
            raise DomainError(f"Cannot decompose {n} mod 2^{k}: odd part is even")
        self.alpha, self.e = result
        self._n = n_mod

    @classmethod
    def from_coords(cls, v: int | None, alpha: int, e: int, k: int = 64) -> DualNumber:
        """Build a DualNumber directly from coordinates.

        Pass v=None to represent zero (infinite valuation).
        """
        obj = cls.__new__(cls)
        obj.k = k
        obj.mask = bitmask(k)
        obj._ord = 1 << (k - 2)
        if v is None:
            obj.is_zero = True
            obj.v = None
            obj.alpha = 0
            obj.e = 0
            obj._n = 0
            return obj
        if v < 0:
            raise DomainError(f"v must be non-negative, got {v}")
        if v >= k:
            obj.is_zero = True
            obj.v = None
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
        assert self.v is not None  # guarded by is_zero
        odd = pow(5, self.e, 1 << self.k)
        if self.alpha:
            odd = (-odd) & self.mask
        return (odd << self.v) & self.mask

    def verify(self) -> bool:
        """Round-trip check: coordinates -> integer matches stored integer."""
        if self.is_zero:
            return self._n == 0
        return self._to_int() == self._n

    @property
    def value(self) -> int:
        """The integer n mod 2^k."""
        return self._n

    def coords(self) -> tuple[int | None, int, int]:
        """Return (v, alpha, e).  v is None for zero."""
        return (self.v, self.alpha, self.e)

    def __repr__(self) -> str:
        if self.is_zero:
            return f"DualNumber(0, k={self.k})"
        sign = "-" if self.alpha else "+"
        return f"DualNumber({self._n}, k={self.k})  =  2^{self.v} . {sign}5^{self.e}"


# -- TwoAdicProcessor ---------------------------------------------------------


class TwoAdicProcessor:
    """
    Arithmetic on DualNumbers in Z/2^k.

    Multiplication, inversion, and exponentiation operate directly in
    coordinate space, making the group structure explicit.
    """

    def __init__(self, k: int = 64, verify: bool = True):
        if k < 3:
            raise DomainError("k must be >= 3")
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
        assert a.v is not None and b.v is not None  # guarded by is_zero
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
        if a.is_zero:
            raise NonInvertibleError("Only units (v=0) are invertible")
        assert a.v is not None  # guarded by is_zero
        if a.v != 0:
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
        assert a.v is not None  # guarded by is_zero
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

    def dlog(self, a: DualNumber) -> tuple[int | None, int, int]:
        """Return full coordinate triple (v, alpha, e).  v is None for zero."""
        self._check(a)
        return a.coords()

    @contextmanager
    def verification(self, enabled: bool = True) -> Iterator[None]:
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


# -- self-test ----------------------------------------------------------------


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
        failures.append(f"mul: 3x5 -> {d15.value}, expected 15")

    inv3 = proc.inv(d3)
    if proc.mul(d3, inv3).value != 1:
        failures.append("inv: 3 x 3^{-1} != 1")

    if proc.pow(d3, 4).value != (81 & mask):
        failures.append(f"pow: 3^4 -> {proc.pow(d3, 4).value}, expected {81 & mask}")

    if failures:
        for msg in failures:
            print(f"FAIL  {msg}")
        raise AssertionError(f"{len(failures)} test(s) failed at k={k}")

    if verbose:
        print(f"run_all_tests: all checks passed  (k={k})")
