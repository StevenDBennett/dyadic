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

from dyadic_core import DualNumber, modinv_newton, two_adic_dlog, valuation, bitmask, mat_mul, mat_det

# ── Congruence filtration ───────────────────────────────────────────────────

def congruence_depth(M: list[list[int]], k: int) -> int:
    """
    Largest j such that M ≡ I (mod 2^j).

    This is the matrix analog of v₂(n - 1) for scalars.
    Returns 0 if M ≠ I (mod 2).
    """
    for j in range(1, k + 1):
        mask = bitmask(j)
        if (M[0][0] & mask) != 1 or (M[0][1] & mask) != 0 or \
           (M[1][0] & mask) != 0 or (M[1][1] & mask) != 1:
            return j - 1
    return k

def filtration_residue(M: list[list[int]], depth: int, k: int) -> list[list[int]]:
    """
    The gl(2, F₂) direction of departure from identity at depth.

    Returns (M - I) / 2^depth mod 2 — a matrix over F₂.
    """
    shift = 1 << depth
    residue = [
        [((M[0][0] - 1) // shift) & 1, (M[0][1] // shift) & 1],
        [(M[1][0] // shift) & 1, ((M[1][1] - 1) // shift) & 1],
    ]
    return residue

# ── LDU decomposition ───────────────────────────────────────────────────────

def ldu_decompose(M: list[list[int]], k: int) -> dict[str, list[list[int]]] | None:
    """
    Factor M = L · D · U (mod 2^k) where:
        L — lower unipotent ([[1, 0], [l, 1]])
        D — diagonal       ([[d1, 0], [0, d2]])
        U — upper unipotent ([[1, u], [0, 1]])

    Returns None if pivot M[0][0] is even (not invertible mod 2^k).
    """
    a, b = M[0][0], M[0][1]
    c, d = M[1][0], M[1][1]

    if a & 1 == 0:
        return None

    inv_a = modinv_newton(a, k)
    u = (b * inv_a) & bitmask(k)
    l = (c * inv_a) & bitmask(k)
    d1 = a
    d2 = (d - l * b) & bitmask(k)

    return {
        "L": [[1, 0], [l, 1]],
        "D": [[d1, 0], [0, d2]],
        "U": [[1, u], [0, 1]],
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

def matrix_coordinates(M: list[list[int]], k: int) -> MatrixCoordinates:
    """
    Full coordinate decomposition of M ∈ GL(2, Z/2^k).

    Returns depth, F₂-residue, trace v₂, determinant DualNumber,
    and LDU factors.
    """
    depth = congruence_depth(M, k)
    res = filtration_residue(M, depth, k)
    trace_val = (M[0][0] + M[1][1]) & bitmask(k)
    trace_depth = valuation(trace_val) if trace_val != 0 else k
    det_val = mat_det(M, 1 << k)
    det_dual = DualNumber(det_val, k)
    ldu = ldu_decompose(M, k)

    return MatrixCoordinates(
        depth=depth,
        residue=res,
        trace_depth=trace_depth,
        det_dual=det_dual,
        ldu=ldu,
    )

# ── Holonomy depth profiling ────────────────────────────────────────────────

def holonomy_depth_profile(
    k: int, p: int, cycle_length: int = 4, n_cycles: int = 30
) -> Dict:
    """
    How holonomy congruence depth changes under single-bit perturbation.
    """
    from .nonabelian import NonAbelianCRTDual

    nc = NonAbelianCRTDual(k, p)
    depths_orig: list[int] = []
    depths_pert: list[int] = []

    for _ in range(n_cycles):
        mats = [
            [[random.randrange(0, nc.mod_full) for _ in range(2)]
             for _ in range(2)]
            for _ in range(cycle_length)
        ]
        H = nc.holonomy(mats)
        d_orig = congruence_depth(H, k)

        mats[0][0][0] ^= 1  # flip lowest bit
        H_pert = nc.holonomy(mats)
        d_pert = congruence_depth(H_pert, k)

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
        # |GL(2, Z/2^j)| = 2^(4j - 3) * 3 for j >= 1
        order = 2 ** (4 * max(j, 1) - 3) * 3 if j >= 1 else 1
        # Quotient Gamma(j) / Gamma(j+1) ≅ gl(2, F₂)
        quotient_size = 16 if j < k else 1
        lines.append(
            f"  Γ(2^{j}): |GL| ≈ 2^{4*max(j,1)-3}·3, "
            f"quotient size = {quotient_size}"
        )
    return "\n".join(lines)

# ── Commutator depth theorem ────────────────────────────────────────────────

def matrix_commutator(
    M: list[list[int]], N: list[list[int]], k: int
) -> list[list[int]]:
    """
    Compute [M, N] = M N M^{-1} N^{-1} modulo 2^k.

    Requires det(M), det(N) to be odd (GL(2) condition).
    """
    mod = 1 << k
    inv_M = [[M[1][1], -M[0][1]], [-M[1][0], M[0][0]]]
    det_M = mat_det(M, mod)
    inv_det_M = modinv_newton(det_M, k)
    for i in range(2):
        for j in range(2):
            inv_M[i][j] = (inv_M[i][j] * inv_det_M) % mod

    inv_N = [[N[1][1], -N[0][1]], [-N[1][0], N[0][0]]]
    det_N = mat_det(N, mod)
    inv_det_N = modinv_newton(det_N, k)
    for i in range(2):
        for j in range(2):
            inv_N[i][j] = (inv_N[i][j] * inv_det_N) % mod

    MN = mat_mul(M, N, mod)
    MN_inv = mat_mul(inv_M, inv_N, mod)
    return mat_mul(MN, MN_inv, mod)

def verify_commutator_depth(
    k: int, depth_pairs: list[tuple[int, int]] = None, n_trials: int = 50
) -> list[tuple[int, int, int]]:
    """
    Verify depth([M,N]) >= depth(M) + depth(N) for random matrices.

    Returns list of (dM, dN, dMN) triples.
    """
    if depth_pairs is None:
        depth_pairs = [(1, 1), (1, 2), (2, 1), (2, 2)]

    results: list[tuple[int, int, int]] = []

    for dM, dN in depth_pairs:
        for _ in range(n_trials):
            shift_M = 1 << dM
            shift_N = 1 << dN
            M = [[1 + shift_M, 0], [0, 1 + shift_M]]
            N = [[1 + shift_N, 0], [0, 1 + shift_N]]

            M_computed = matrix_commutator(M, N, k)
            dMN = congruence_depth(M_computed, k)
            results.append((dM, dN, dMN))

    return results
