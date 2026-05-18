"""
fourier.py
----------
Discrete Fourier analysis of the Newton step-count function.

Let N = 2^(k-2) and let h : Z/N → {0, 1, 2, 3} be the step-count
function giving the number of Newton iterations needed for a seed
e to converge to the true root e_true.

Because h(e) depends only on v₂(e - e_true), its DFT has power
only at dyadic frequencies j = 0, N/2, N/4, N/8, ...  This module
computes the analytic closed-form coefficients and verifies them
against the numeric FFT.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from dyadic_core import bitmask, two_adic_log5

from dyadic_math._step import newton_step_core


def step_count_fn(k: int, e_true: int, max_steps: int = 10) -> np.ndarray:
    """
    Numerically compute the step-count function h(e).

    h(e) is the number of Newton iterations required for seed e
    to converge to e_true.  Returns an integer array of length N.

    The function only takes values in {0, 1, 2, 3} for valid inputs;
    h = max_steps + 1 indicates failure to converge within budget.
    """
    order = 1 << (k - 2)
    a_val = pow(5, e_true, 1 << k)
    log5_unit = two_adic_log5(k) >> 2
    mask = bitmask(k)
    exp_mask = order - 1
    h = np.zeros(order, dtype=np.int32)

    for e_seed in range(order):
        e = e_seed
        step_count = 0
        while e != e_true and step_count < max_steps:
            e = newton_step_core(5, e, a_val, k, log5_unit, mask, exp_mask)
            step_count += 1
        if e == e_true:
            h[e_seed] = step_count
        else:
            h[e_seed] = max_steps + 1  # did not converge
    return h


def analytic_step_count(k: int, e_true: int) -> np.ndarray:
    """
    Closed-form O(N) construction of h(e) using ultrametric balls.

    For seed e:
        h = 0  if e == e_true
        h = 1  if v₂(e - e_true) >= 2  (i.e. in ball of radius ≥ 4)
        h = 2  if v₂(e - e_true) == 1  (distance exactly 2)
        h = 3  if v₂(e - e_true) == 0  (distance odd)
    """
    order = 1 << (k - 2)
    h = np.full(order, 3, dtype=np.int32)
    offset = np.arange(order, dtype=np.intp)
    diff = (offset - e_true) & (order - 1)

    mask0 = diff == 0
    h[mask0] = 0

    if k >= 2:
        mask1 = (diff & 3) == 0  # v2 >= 2
        h[mask1 & ~mask0] = 1

    if k >= 3:
        mask2 = (diff & 1) == 0  # v2 >= 1 but < 2
        h[mask2 & ~mask0] = 2

    return h


def dft(f: np.ndarray) -> Any:
    """Discrete Fourier transform (numpy wrapper)."""
    return np.fft.fft(f)


def power_spectrum(f: np.ndarray) -> Any:
    """Power spectrum |DFT(f)|^2."""
    return np.abs(dft(f)) ** 2


def dyadic_coefficients(f: np.ndarray) -> dict[str, complex]:
    """
    Extract DFT coefficients at dyadic frequencies.

    Only dyadic frequencies j = 0, N/2, N/4, ... have non-zero
    power for the step-count function.
    """
    length = len(f)
    spectrum = dft(f)
    coeffs: dict[str, complex] = {"DC": spectrum[0]}

    m = 1
    while length // (2 * m) > 0:
        j = length // (2 * m)
        if j < len(spectrum):
            coeffs[f"N/{2 * m}"] = spectrum[j]
        m *= 2
    return coeffs


def analytic_coefficients(k: int) -> dict[str, complex]:
    """
    Closed-form Fourier coefficients for the step-count function.

    Returns dict with DC, N/2, N/4 components derived analytically.
    """
    order = 1 << (k - 2)

    n_dc = order - order // 4
    dc_coeff = complex(n_dc, 0)

    # N/2: contrast between even and odd exponents
    n_half = order // 2
    half_coeff = complex(n_half, 0)

    n_quarter = order // 4
    quarter_coeff = complex(-n_quarter, 0)

    return {
        "DC": dc_coeff,
        "N/2": half_coeff,
        "N/4": quarter_coeff,
    }


def fourier_summary(k: int, e_true: int) -> str:
    """
    Complete Fourier analysis of the step-count function.

    Compares analytic vs numeric coefficients.
    """
    h_num = step_count_fn(k, e_true)
    h_an = analytic_step_count(k, e_true)
    match = np.all(h_num == h_an)

    spec = power_spectrum(h_num)
    dyadic = dyadic_coefficients(h_num)
    analytic = analytic_coefficients(k)

    lines = [
        f"Fourier Analysis of Step-Count Function  (k={k}, N={1 << (k - 2)})",
        f"  Numeric == Analytic: {match}",
        f"  Non-zero frequencies: {len([v for v in spec if v > 0])}",
        "",
        "  Dyadic Coefficients:",
    ]
    for name, coeff in dyadic.items():
        a_coeff = analytic.get(name, 0)
        lines.append(f"    {name:8s}: numeric={coeff.real:8.1f}  analytic={a_coeff.real:8.1f}")

    return "\n".join(lines)


def ultrametric_uncertainty(k: int) -> str:
    """
    Statement of the 2-adic Fourier uncertainty principle.

    A function on Z/N with ultrametric support (ball of radius
    2^(-m)) has Fourier support only at dyadic frequencies.
    The product |support|·|freq_support| = N.
    """
    order = 1 << (k - 2)
    return (
        f"2-adic Fourier Uncertainty (k={k}, N={order}):\n"
        f"A function supported on an ultrametric ball of radius 2^(-m)\n"
        f"has Fourier support only at frequencies j = N/2^t for t >= m.\n"
        f"Product |supp| × |freq_supp| = N (saturated)."
    )
