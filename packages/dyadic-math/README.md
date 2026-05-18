# dyadic-math

p-adic Newton dynamics, CRT extensions, Fourier analysis, Mersenne theory.

## Overview

Builds on **dyadic-core** to provide Newton basin analysis, CRT extensions,
p-adic root finding, and the Mersenne Ghost Theorem.

## Features

- **BasinExplorer** — classify Newton seeds by convergence fate (converged, ghost cycle, diverged)
- **GhostHunt** — precision sweep to find quantization cliff thresholds
- **WeightStabilityDiagnostics** — graded 2-adic weight stability diagnostics
- **padic_roots** — Newton, Halley, and higher-order p-adic root finding
- **CRTDualNumber** — Chinese Remainder Theorem extension to `Z/(2^k * p)Z`
- **NonAbelianCRTDual** — GL(2) holonomy with phase alignment experiments

Requires `dyadic-core` and `numpy`.

## Install

```bash
pip install dyadic-math
```

## Quick Start

```python
from dyadic_math import BasinExplorer

explorer = BasinExplorer(k=8, g=5, target_a=3)
portrait = explorer.portrait()
print(f"Converged: {len(portrait['converged'])}, Ghost: {len(portrait['cycle'])}")
```

## License

MIT
