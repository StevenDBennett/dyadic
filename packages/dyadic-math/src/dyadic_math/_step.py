"""
_step.py
--------
Shared Newton step for the 2-adic discrete-log map.

Extracted from basin.py to eliminate the hub dependency: separation.py
and fourier.py import newton_step_core without pulling in the full
BasinExplorer class, its portrait cache, and numpy.
"""

from __future__ import annotations

from dyadic_core import modinv_newton

__all__ = [
    "newton_step_core",
]


def newton_step_core(e: int, a: int, k: int, log5_unit: int, mask: int, exp_mask: int) -> int:
    """
    Single Newton iteration for the map e ↦ e - (5^e - a) / (5^e · log5_unit).

    Generator is hardcoded to 5 — the discrete-log infrastructure
    (two_adic_dlog, _DLOG10_LUT, two_adic_log5) only supports base 5.
    This is the core computation used by BasinExplorer, separation.py,
    and fourier.py.  Extracted to eliminate three-way code duplication.
    """
    g_val = pow(5, e, 1 << k)
    f_val = (g_val - a) & mask
    df_val = (g_val * log5_unit) & exp_mask
    df_inv = modinv_newton(df_val, k - 2)
    delta = ((f_val >> 2) * df_inv) & exp_mask
    return (e - delta) & exp_mask
