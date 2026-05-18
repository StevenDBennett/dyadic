# dyadic-core

> This is a dyadic (2-package library) document, ported from the dual-view project.

**dyadic-core** is the 2-adic arithmetic foundation. It provides the dual-view decomposition of integers modulo powers of 2, fast modular arithmetic via Newton lifting, and the Mahler calculus for integer-valued functions.

## Modules

| Module | Description |
|--------|-------------|
| `dyadic_core.core` | DualNumber, TwoAdicProcessor, modinv_newton, padic_exp/log, cliff centre g₀ |
| `dyadic_core.exponent` | ExponentSpace — additive coordinate chart on Z/2^(k-2) |
| `dyadic_core.mahler` | MahlerCalculus — Mahler basis, Dirac/Volterra operators |

## Quick-Start Examples

### DualNumber — the core data type

```python
from dyadic_core import DualNumber

# Decompose 42 modulo 2^8
n = DualNumber(42, k=8)
print(f"v={n.v}, alpha={n.alpha}, e={n.e}")
print(f"Reconstruction: {n.verify()}")  # True
```

### Modular Inverse via Newton Lifting

```python
from dyadic_core import modinv_newton

inv = modinv_newton(3, k=16)   # 3⁻¹ mod 2¹⁶
print(f"3 * inv mod 2^16 = {(3 * inv) & ((1 << 16) - 1)}")  # 1
```

### 2-Adic Exponential and Logarithm

```python
from dyadic_core import padic_exp, padic_log, g0

# The cliff centre: unique unit with log(g₀) = -4
g0_val = g0(k=16)
print(f"log(g0) = {padic_log(g0_val, 16)}")  # -4 mod 2^16

# Exponential
x = 4  # must have v2(x) >= 2
exp_x = padic_exp(x, 16)
print(f"exp(4) mod 2^16 = {exp_x}")
```

### TwoAdicProcessor — coordinate-space arithmetic

```python
from dyadic_core import TwoAdicProcessor, DualNumber

proc = TwoAdicProcessor(k=8)
a = DualNumber(3, k=8)
b = DualNumber(5, k=8)
c = proc.mul(a, b)
print(f"3 * 5 = {c.value}")  # 15
```

### ExponentSpace — additive chart on Z/2^(k-2)

```python
from dyadic_core import ExponentSpace

es = ExponentSpace(g=5, k=8)
for e in range(4):
    val = es.lift(e)
    print(f"5^{e} mod 256 = {val}")
```

### MahlerCalculus — binomial basis

```python
from dyadic_core import MahlerCalculus

# Convert function values to Mahler coefficients
f = [0, 1, 4, 9, 16]  # f(x) = x^2 for x = 0..4
coeffs = MahlerCalculus.to_mahler(f, max_degree=4)
print(f"Mahler coeffs: {coeffs}")

# Dirac operator lowers degree
d_coeffs = MahlerCalculus.dirac_operator(coeffs)
print(f"D(coeffs): {d_coeffs}")

# Evaluate at a point
val = MahlerCalculus.from_mahler(coeffs, 5)
print(f"f(5) = {val}")  # 25
```

## Dependencies

- stdlib only

## Further Reading

- `docs/api.md` — Complete API reference
- `docs/mathematics.md` — Mathematical background
