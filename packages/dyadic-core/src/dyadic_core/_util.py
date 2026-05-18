from __future__ import annotations

__all__ = [
    "WARN_K_LIMIT",
    "MAX_K_ENUMERATE",
    "bitmask",
    "mat_det",
    "mat_mul",
    "valuation",
]

# O(2^k) enumeration limits for Newton basin analysis.
# Enumerating all N = 2^(k-2) seeds becomes expensive past k=16.
WARN_K_LIMIT = 16
MAX_K_ENUMERATE = 20


def bitmask(k: int) -> int:
    """Bitmask of k bits: (1 << k) - 1."""
    return (1 << k) - 1


def valuation(n: int) -> int | None:
    """2-adic valuation v_2(n). Returns None for zero."""
    if n == 0:
        return None
    return (n & -n).bit_length() - 1


def mat_mul(A: list[list[int]], B: list[list[int]], mod: int) -> list[list[int]]:
    """2x2 matrix multiplication modulo mod."""
    return [
        [
            (A[0][0] * B[0][0] + A[0][1] * B[1][0]) % mod,
            (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % mod,
        ],
        [
            (A[1][0] * B[0][0] + A[1][1] * B[1][0]) % mod,
            (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % mod,
        ],
    ]


def mat_det(M: list[list[int]], mod: int) -> int:
    """Determinant of a 2x2 matrix modulo mod."""
    return (M[0][0] * M[1][1] - M[0][1] * M[1][0]) % mod
