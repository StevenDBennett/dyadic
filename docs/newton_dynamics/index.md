> This document is part of the dyadic (3-package monorepo) project, ported from dual-view.

# p-adic Newton Dynamics

## Overview

Analysis of the p-adic Newton dynamics of the rational map

$$N(x) = \frac{2x^3 + 1}{3x^2}$$

which computes cube roots via Newton's method. The map has:

- **Fixed points**: $x^3 = 1$ (cube roots of unity), all **superattracting** ($N'(\zeta) = 0$)
- **Symmetry**: $N(\zeta x) = \zeta N(x)$ for $\zeta^3 = 1$ — all dynatomic polynomials are polynomials in $u = x^3$
- **Critical points**: $x = 0$ (pole), $x^3 = 1$ (fixed points)

## Module Reference

| Module | Description |
|--------|-------------|
| `poly.py` | Polynomial arithmetic (convolution, power, division) |
| `iterates.py` | Newton iterate recurrence $N^d(x) = A_d(x)/B_d(x)$ |
| `dynatomic.py` | Dynatomic polynomial $\Phi_n^*(u)$ via Möbius inversion |
| `clean_primes.py` | Cube-root checking, Tonelli–Shanks, clean prime detection |
| `data.py` | Precomputed coefficients and multiplier data for periods 2–6 |

## Quick Start

```python
from dyadic_math.newton_dynamics import compute_iterates, dynatomic_polynomial

# Compute Newton iterates
iterates = compute_iterates(4)

# Get period-2 dynatomic polynomial
phi2 = dynatomic_polynomial(2, iterates)
# → [2, 5, 20]  = 20u² + 5u + 2

# Get period-4 dynatomic polynomial (degree 24 in u)
phi4 = dynatomic_polynomial(4, iterates)
# → 25 coefficients, constant = 4096 = 2¹²
```

```python
from dyadic_math.newton_dynamics import is_cube, tonelli_shanks

# Cube checking modulo p
is_cube(8, 7)  # True: 2³ ≡ 1 (mod 7) → 8 ≡ 1, and 1 is a cube

# Modular square root
sqrt = tonelli_shanks(5, 11)  # 4² ≡ 5 (mod 11)
```

## Key Results

| Period | $\mu_3(n)$ | deg(x) | deg(u) | $\Phi_n^*(0)$ | $\Phi_n^*(1)$ | $\prod\mu$ |
|--------|------------|--------|--------|---------------|---------------|------------|
| 2 | 6 | 6 | 2 | $2^1$ | $3^3$ | $6^1$ |
| 3 | 24 | 24 | 8 | — | — | — |
| 4 | 72 | 72 | 24 | $2^{12}$ | $3^{36}$ | $6^{12}$ |
| 5 | 240 | 240 | 80 | $2^{40}$ | $3^{120}$ | $6^{40}$ |
| 6 | 696 | 696 | 232 | $2^{116}$ | $3^{348}$ | $6^{116}$ |

## Clean Primes

**{7, 103, 181}** — verified below 100,000. Conjectured to be complete.

## Documents

| Document | Description |
|----------|-------------|
| `index.md` | This file — overview, quick-start, module reference |
| `theorems.md` | All 4 proven theorems with full proofs |
| `computation_data.md` | Period-by-period coefficients, multipliers, Galois data |
| `clean_primes.md` | Clean prime analysis and finiteness argument |
| `status.md` | Completed items, open problems, research log |

*Research completed 2026-05-11. Ported from dual-view.*
