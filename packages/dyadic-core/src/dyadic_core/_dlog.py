"""
_dlog.py
--------
2-adic discrete logarithm: 5^e == a (mod 2^k).
"""

from __future__ import annotations

from functools import lru_cache
from typing import TypedDict, cast

from dyadic_core._exceptions import DomainError
from dyadic_core._modular import modinv_newton
from dyadic_core._series import padic_log
from dyadic_core._util import bitmask, valuation

__all__ = [
    "DLogNewtonStep",
    "dlog_bootstrap",
    "dlog_newton_step",
    "dlog_residual_tracking",
    "two_adic_dlog",
    "two_adic_log5",
]

# -- 10-bit discrete log LUT --------------------------------------------------
# Precomputed: for each odd a == 1 (mod 4) in [0, 1024), store e s.t. 5^e == a (mod 1024).
# 256 entries ~ 2 KB.  Gives O(1) bootstrap to 8-bit precision.
_DLOG10_LUT: list[int] = [-1] * 1024
_DLOG10_LUT_INIT = {pow(5, e, 1024): e for e in range(256)}
for _a, _e in _DLOG10_LUT_INIT.items():
    _DLOG10_LUT[_a] = _e
del _DLOG10_LUT_INIT


# -- 2-adic logarithm of 5 ----------------------------------------------------


@lru_cache(maxsize=256)
def two_adic_log5(k: int) -> int:
    """
    Compute the 2-adic logarithm of 5, truncated to k bits.

    Delegates to padic_log(5, k) — the general 2-adic log series.
    Cached per k for performance (used as derivative scaling factor
    in Newton dlog).
    """
    return padic_log(5, k)


# -- discrete-log helpers (private) -------------------------------------------


def dlog_bootstrap(a: int, k: int) -> int:
    """
    Bit-by-bit dlog for a == 1 (mod 8).
    Returns e with 5^e == a (mod 2^k).  O(k) steps; used to seed Newton.
    """
    if k <= 2:
        return 0
    mask = bitmask(k)
    a &= mask
    e, pow5 = (0, 1) if (a & 7) == 1 else (1, 5)
    mult = 25 & mask
    for n in range(3, k):
        if ((a >> n) & 1) != ((pow5 >> n) & 1):
            e |= 1 << (n - 2)
            pow5 = (pow5 * mult) & mask
        mult = (mult * mult) & mask
    return e


def dlog_newton(a: int, k: int, L: int | None = None) -> int:
    """
    Newton-lifted dlog: 5^e == a (mod 2^k).  Requires a == 1 (mod 4).

    Bootstrap uses an 8-bit LUT for k <= 34 (covers 6-bit precision),
    otherwise uses bit-by-bit to k/2.  Per the Mersenne bootstrap
    optimality analysis, k/2 saves ceil(log2(sqrt(k))) - 1 Newton steps
    vs sqrt(k).
    """
    if k <= 2:
        return 0
    if a & 3 != 1:
        raise DomainError("dlog_newton requires a == 1 (mod 4)")
    mask_full = bitmask(k)
    a &= mask_full
    if L is None:
        L = two_adic_log5(k)
    L_unit = L >> 2

    if k <= 34:
        a_b = a & 0x3FF
        e_raw = _DLOG10_LUT[a_b] if _DLOG10_LUT[a_b] >= 0 else 0
        eprec = min(8, k - 2)
        e = e_raw & bitmask(eprec)
    else:
        bootstrap_k = max(4, k // 2 + 2)
        eprec = min(bootstrap_k - 2, k - 2)
        e = dlog_bootstrap(a, bootstrap_k) & bitmask(eprec)

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = bitmask(bits)
        emask = bitmask(new_eprec)

        e, _ = dlog_newton_step(e, a, L_unit, bits, mask, emask, new_eprec)
        eprec = new_eprec

    return e


def dlog_newton_step(
    e: int, a: int, L_unit: int, bits: int, mask: int, emask: int, new_eprec: int
) -> tuple[int, int]:
    """
    Single Newton step for 2-adic discrete log: refine e so that 5^e == a (mod 2^bits).

    Returns (new_e, delta).
    """
    pow5e = pow(5, e, 1 << bits)
    f = (pow5e - a) & mask
    df_unit = (pow5e * L_unit) & emask
    df_inv = modinv_newton(df_unit, new_eprec)
    delta = ((f >> 2) * df_inv) & emask
    e_new = (e - delta) & emask
    return e_new, delta


# -- residual tracking (diagnostic) -------------------------------------------


class DLogNewtonStep(TypedDict):
    bits: int
    eprec_before: int
    eprec_after: int
    tau_before: int
    v2_before: int | None
    tau_after: int
    v2_after: int | None
    delta: int


def dlog_residual_tracking(
    a: int, k: int, L: int | None = None
) -> tuple[int, list[DLogNewtonStep]]:
    """
    Viglietta discrete log with residual tracking at each Newton step.

    Returns (e, history) where history contains the normalised residual
    (tau = (5^e * a^{-1} - 1) / 4) and its 2-adic valuation before and
    after each Newton correction.  The doubling of v2(tau) confirms the
    quadratic convergence law.

    Requires a == 1 (mod 4).
    """
    if k <= 2:
        return 0, []
    if a & 3 != 1:
        raise DomainError("dlog_residual_tracking requires a == 1 (mod 4)")
    mask_full = bitmask(k)
    a &= mask_full
    if L is None:
        L = two_adic_log5(k)
    L_unit = L >> 2
    a_inv_full = modinv_newton(a, k)

    bootstrap_k = max(4, k // 2 + 2)
    eprec = min(bootstrap_k - 2, k - 2)
    e = dlog_bootstrap(a, bootstrap_k) & bitmask(eprec)
    history: list[DLogNewtonStep] = []

    while eprec < k - 2:
        new_eprec = min(2 * eprec, k - 2)
        bits = new_eprec + 2
        mask = bitmask(bits)
        emask = bitmask(new_eprec)

        pow5e = pow(5, e, 1 << bits)
        a_inv = a_inv_full & mask
        rho_before = (pow5e * a_inv) & mask
        tau_before = ((rho_before - 1) & mask) >> 2
        v2_before = valuation(tau_before)

        e, delta = dlog_newton_step(e, a, L_unit, bits, mask, emask, new_eprec)

        pow5e = pow(5, e, 1 << bits)
        rho_after = (pow5e * a_inv) & mask
        tau_after = ((rho_after - 1) & mask) >> 2
        v2_after = valuation(tau_after)

        history.append(
            cast(
                DLogNewtonStep,
                {
                    "bits": bits,
                    "eprec_before": eprec,
                    "eprec_after": new_eprec,
                    "tau_before": int(tau_before),
                    "v2_before": v2_before,
                    "tau_after": int(tau_after),
                    "v2_after": v2_after,
                    "delta": int(delta),
                },
            )
        )
        eprec = new_eprec

    return e, history


# -- public discrete-log API --------------------------------------------------


def two_adic_dlog(a: int, k: int, L: int | None = None) -> tuple[int, int] | None:
    """
    Decompose the odd part of a as (-1)^alpha . 5^e (mod 2^k).

    Parameters
    ----------
    a : int   -- integer to decompose (any residue mod 2^k)
    k : int   -- bit precision (>= 3)
    L : int   -- precomputed two_adic_log5(k); computed if not supplied

    Returns
    -------
    (alpha, e)  with alpha in {0,1} and e in [0, 2^(k-2)), or
    None        if a is even.
    """
    if a & 1 == 0:
        return None
    mask = bitmask(k)
    a &= mask
    alpha = (a >> 1) & 1
    if alpha:
        a = (-a) & mask
    e = dlog_newton(a, k, L)
    return alpha, e
