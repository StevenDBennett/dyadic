# dyadic — 2-Adic Arithmetic & Newton Dynamics

A two-package Python library for exact 2-adic arithmetic and Newton dynamics
analysis of the discrete-log map.

```
n = 2^v · (-1)^α · 5^e  (mod 2^k)

v = v₂(n)        — 2-adic valuation
α ∈ {0, 1}      — sign sector
e ∈ [0, 2^(k-2)) — discrete logarithm base 5
```

**dyadic-core** provides the arithmetic: DualNumber decomposition, modular
inverse via Newton lifting, discrete logarithms, and p-adic exp/log — all
stdlib.

**dyadic-math** provides Newton dynamics tools built on top: basin
exploration, ghost attractor detection, the Mersenne Ghost Theorem,
trajectory separation analysis, and weight stability diagnostics.

## Packages

| Package | Description | Dependencies |
|---------|-------------|-------------|
| `dyadic-core` | 2-adic arithmetic: DualNumber, modinv, dlog, padic exp/log | stdlib only |
| `dyadic-math` | Newton dynamics: basin analysis, ghost detection, Mersenne theorems, stability diagnostics | `dyadic-core`, `numpy` |

## Install

```bash
pip install -e packages/dyadic-core
pip install -e packages/dyadic-math
```

## Quick Start

```python
from dyadic_core import DualNumber, modinv_newton

n = DualNumber(42, k=16)
print(f"v={n.v}, alpha={n.alpha}, e={n.e}")  # v=1, alpha=0, e=...

inv = modinv_newton(3, k=16)
print(f"3^{-1} mod 2^16 = {inv}")
```

```python
from dyadic_math import BasinExplorer

explorer = BasinExplorer(k=8, g=5, target_a=3)
portrait = explorer.portrait()
print(f"Converged: {len(portrait['converged'])}, Ghost: {len(portrait['cycle'])}")
```

```python
import numpy as np
from dyadic_math import WeightStabilityDiagnostics

...

thermo = WeightStabilityDiagnostics(k=8)
result = thermo.analyse(W)
print(f"Alpha fraction: {result['alpha_fraction']:.3f}")
```

## Key Results

- **Mersenne Ghost Theorem**: For Mersenne numbers `w = 2^n - 1`, the
  quantization cliff occurs at precision `k* = n + 2`. Verified n=3..11.
- **Trajectory Separation**: Two Newton trajectories starting from seeds
  separated by `2^s` are identical for exactly `ceil(log₂(s)) - 1` steps.
- **Basin Dichotomy**: The α=0 sector has no ghost attractors.

See `docs/` for full theorem proofs and API reference.

## Run tests

```bash
pytest packages/dyadic-core/tests -v
pytest packages/dyadic-math/tests -v
```
