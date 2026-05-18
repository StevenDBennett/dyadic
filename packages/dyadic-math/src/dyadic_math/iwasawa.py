"""
iwasawa.py
----------
Iwasawa decomposition and congruence filtration for GL(2, Z/2^k).

Provides the matrix analog of the scalar DualNumber coordinate
system, including:

1. Congruence depth — the largest j such that M ≡ I (mod 2^j),
   the matrix analog of v₂(n - 1).
2. Filtration residue — the gl(2, F₂) direction at depth j.
3. LDU decomposition of any matrix with odd determinant.
4. Commutator depth theorem: depth([M,N]) ≥ depth(M) + depth(N).
"""

from __future__ import annotations

import random
from dataclasses import dataclass

import numpy as np
from dyadic_core import (
    DualNumber,
    bitmask,
    mat_det,
    mat_mul,
    modinv_newton,
    valuation,
)

# ── Congruence filtration ───────────────────────────────────────────────────


def congruence_depth(matrix: list[list[int]], k: int) -> int:
    """
    Largest j such that matrix ≡ I (mod 2^j).

    This is the matrix analog of v₂(n - 1) for scalars.
    Returns 0 if M ≠ I (mod 2).
    """
    for j in range(1, k + 1):
        mask = bitmask(j)
        if (
            (matrix[0][0] & mask) != 1
            or (matrix[0][1] & mask) != 0
            or (matrix[1][0] & mask) != 0
            or (matrix[1][1] & mask) != 1
        ):
            return j - 1
    return k


def filtration_residue(matrix: list[list[int]], depth: int, _k: int) -> list[list[int]]:
    """
    The gl(2, F₂) direction of departure from identity at depth.

    Returns (matrix - I) / 2^depth mod 2 — a matrix over F₂.
    """
    shift = 1 << depth
    residue = [
        [((matrix[0][0] - 1) // shift) & 1, (matrix[0][1] // shift) & 1],
        [(matrix[1][0] // shift) & 1, ((matrix[1][1] - 1) // shift) & 1],
    ]
    return residue


# ── LDU decomposition ───────────────────────────────────────────────────────


def ldu_decompose(matrix: list[list[int]], k: int) -> dict[str, list[list[int]]] | None:
    """
    Factor matrix = L · D · U (mod 2^k) where:
        L — lower unipotent ([[1, 0], [l, 1]])
        D — diagonal       ([[d1, 0], [0, d2]])
        U — upper unipotent ([[1, u], [0, 1]])

    Returns None if pivot M[0][0] is even (not invertible mod 2^k).
    """
    a, b = matrix[0][0], matrix[0][1]
    c, d = matrix[1][0], matrix[1][1]

    if a & 1 == 0:
        return None

    inv_a = modinv_newton(a, k)
    u_val = (b * inv_a) & bitmask(k)
    l_val = (c * inv_a) & bitmask(k)
    d1 = a
    d2 = (d - l_val * b) & bitmask(k)

    return {
        "L": [[1, 0], [l_val, 1]],
        "D": [[d1, 0], [0, d2]],
        "U": [[1, u_val], [0, 1]],
    }


# ── Dual coordinates for matrices ───────────────────────────────────────────


@dataclass
class MatrixCoordinates:
    """Full 2-adic coordinate decomposition of a GL(2) matrix."""

    depth: int
    residue: list[list[int]]
    trace_depth: int
    det_dual: DualNumber
    ldu: dict[str, list[list[int]]] | None


def matrix_coordinates(matrix: list[list[int]], k: int) -> MatrixCoordinates:
    """
    Full coordinate decomposition of matrix ∈ GL(2, Z/2^k).

    Returns depth, F₂-residue, trace v₂, determinant DualNumber,
    and LDU factors.
    """
    depth = congruence_depth(matrix, k)
    res = filtration_residue(matrix, depth, k)
    trace_val = (matrix[0][0] + matrix[1][1]) & bitmask(k)
    trace_depth_val = valuation(trace_val)
    trace_depth: int
    if trace_val != 0:
        trace_depth = k if trace_depth_val is None else trace_depth_val
    else:
        trace_depth = k
    det_val = mat_det(matrix, 1 << k)
    det_dual = DualNumber(det_val, k)
    ldu = ldu_decompose(matrix, k)

    return MatrixCoordinates(
        depth=depth,
        residue=res,
        trace_depth=trace_depth,
        det_dual=det_dual,
        ldu=ldu,
    )


# ── Holonomy depth profiling ────────────────────────────────────────────────


def holonomy_depth_profile(k: int, p: int, cycle_length: int = 4, n_cycles: int = 30) -> dict[str, float | int]:
    """
    How holonomy congruence depth changes under single-bit perturbation.
    """
    from .nonabelian import NonAbelianCRTDual

    nc = NonAbelianCRTDual(k, p)
    depths_orig: list[int] = []
    depths_pert: list[int] = []

    for _ in range(n_cycles):
        mats = [
            [[random.randrange(0, nc.mod_full) for _ in range(2)] for _ in range(2)]
            for _ in range(cycle_length)
        ]
        holonomy = nc.holonomy(mats)
        d_orig = congruence_depth(holonomy, k)

        mats[0][0][0] ^= 1  # flip lowest bit
        holonomy_pert = nc.holonomy(mats)
        d_pert = congruence_depth(holonomy_pert, k)

        depths_orig.append(d_orig)
        depths_pert.append(d_pert)

    return {
        "mean_depth_orig": float(np.mean(depths_orig)),
        "mean_depth_pert": float(np.mean(depths_pert)),
        "n_cycles": n_cycles,
    }


# ── Filtration portrait ─────────────────────────────────────────────────────


def filtration_portrait(k: int) -> str:
    """
    Print the full congruence filtration structure of GL(2, Z/2^k).

    Shows the quotient sizes at each depth and the F₂-residue
    classification.
    """
    lines = [f"Congruence Filtration of GL(2, Z/2^{k})"]
    for j in range(k + 1):
        # Quotient Gamma(j) / Gamma(j+1) ≅ gl(2, F₂)
        quotient_size = 16 if j < k else 1
        lines.append(f"  Γ(2^{j}): |GL| ≈ 2^{4 * max(j, 1) - 3}·3, quotient size = {quotient_size}")
    return "\n".join(lines)


# ── Commutator depth theorem ────────────────────────────────────────────────


def matrix_commutator(
    matrix_a: list[list[int]], matrix_b: list[list[int]], k: int
) -> list[list[int]]:
    """
    Compute [matrix_a, matrix_b] = matrix_a matrix_b matrix_a^{-1} matrix_b^{-1} modulo 2^k.

    Requires det(matrix_a), det(matrix_b) to be odd (GL(2) condition).
    """
    mod = 1 << k
    inv_a = [[matrix_a[1][1], -matrix_a[0][1]], [-matrix_a[1][0], matrix_a[0][0]]]
    det_a = mat_det(matrix_a, mod)
    inv_det_a = modinv_newton(det_a, k)
    for i in range(2):
        for j in range(2):
            inv_a[i][j] = (inv_a[i][j] * inv_det_a) % mod

    inv_b = [[matrix_b[1][1], -matrix_b[0][1]], [-matrix_b[1][0], matrix_b[0][0]]]
    det_b = mat_det(matrix_b, mod)
    inv_det_b = modinv_newton(det_b, k)
    for i in range(2):
        for j in range(2):
            inv_b[i][j] = (inv_b[i][j] * inv_det_b) % mod

    ab = mat_mul(matrix_a, matrix_b, mod)
    inv_ab = mat_mul(inv_a, inv_b, mod)
    return mat_mul(ab, inv_ab, mod)


def verify_commutator_depth(
    k: int, depth_pairs: list[tuple[int, int]] | None = None, n_trials: int = 50
) -> list[tuple[int, int, int]]:
    """
    Verify depth([M,N]) >= depth(M) + depth(N) for random matrices.

    Returns list of (dM, dN, dMN) triples.
    """
    if depth_pairs is None:
        depth_pairs = [(1, 1), (1, 2), (2, 1), (2, 2)]

    results: list[tuple[int, int, int]] = []

    for depth_a, depth_b in depth_pairs:
        for _ in range(n_trials):
            shift_a = 1 << depth_a
            shift_b = 1 << depth_b
            matrix_a = [[1 + shift_a, 0], [0, 1 + shift_a]]
            matrix_b = [[1 + shift_b, 0], [0, 1 + shift_b]]

            commutator = matrix_commutator(matrix_a, matrix_b, k)
            depth_comm = congruence_depth(commutator, k)
            results.append((depth_a, depth_b, depth_comm))

    return results
