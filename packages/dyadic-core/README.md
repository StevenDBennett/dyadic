# dyadic-core

2-adic arithmetic foundation for the **dyadic** framework.

## Overview

Implements the dual-view decomposition `n = 2^v * (-1)^alpha * 5^e (mod 2^k)`,
separating magnitude (valuation `v`), sign sector (`alpha`), and discrete-logarithm
structure (`e`). All arithmetic operates directly in coordinate space.

## Features

- **DualNumber** — coordinate triple `(v, alpha, e)` for any integer mod `2^k`
- **TwoAdicProcessor** — multiplication, inversion, and powering in coordinate space
- **modinv_newton** — quadratic Newton-lifted modular inverse mod `2^k`
- **two_adic_dlog** — discrete logarithm base 5 with LUT bootstrap + Newton refinement
- **padic_exp / padic_log** — general 2-adic exponential and logarithm
- **MahlerCalculus** — binomial basis operator calculus (Dirac/Volterra operators)
- **ExponentSpace** — additive coordinate chart on the cyclic group `Z/2^(k-2)`

Zero dependencies — stdlib only.

## Install

```bash
pip install dyadic-core
```

## Quick Start

```python
from dyadic_core import DualNumber, modinv_newton

n = DualNumber(42, k=16)
print(f"v={n.v}, alpha={n.alpha}, e={n.e}")

inv = modinv_newton(3, k=16)
print(f"3^-1 mod 2^16 = {inv}")
```

## License

MIT
