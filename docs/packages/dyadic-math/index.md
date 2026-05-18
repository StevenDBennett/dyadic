# dyadic-math

> This is a dyadic (2-package library) document, ported from the dual-view project.

**dyadic-math** provides p-adic Newton dynamics, CRT extensions, Fourier analysis, Mersenne theory, and the Mersenne Ghost Theorem.

## Modules

| Module | Description |
|--------|-------------|
| `dyadic_math.basin` | BasinExplorer, GhostHunt — Newton basin landscape analysis |
| `dyadic_math.crt` | CRT extension to Z/(2^k·p)Z |
| `dyadic_math.nonabelian` | GL(2) gauge theory with phase alignment |
| `dyadic_math.separation` | Trajectory Separation Theorem — exact zero-variance |
| `dyadic_math.fourier` | DFT of Newton step-count function, 2-adic Fourier uncertainty |
| `dyadic_math.padic_roots` | Multi-order p-adic root finding (Newton, Halley, composed) |
| `dyadic_math.iwasawa` | GL(2) congruence filtration, LDU decomposition |
| `dyadic_math.mersenne` | Mersenne Ghost Theorem, bootstrap optimality, cliff constants |
| `dyadic_math.isometry` | Exponential map isometry |
| `dyadic_math.thermodynamics` | SeedThermodynamics — graded 2-adic weight stability diagnostics |

## Quick-Start Examples

### BasinExplorer — ghost detection

```python
from dyadic_math import BasinExplorer

# Explore Newton basin for target 3 mod 2^8 with generator 5
k, g, target = 8, 5, 3
explorer = BasinExplorer(k, g, target)

# Classify a single seed
fate, value, path = explorer.classify(e0=0)
print(f"Seed 0: fate={fate}, converged to e={value}")

# Full portrait
portrait = explorer.portrait()
print(f"Converged: {len(portrait['converged'])} seeds")
print(f"Ghost cycle: {len(portrait['cycle'])} seeds")
print(f"Diverged: {len(portrait['diverged'])} seeds")
```

### Mersenne Ghost Theorem

```python
from dyadic_math import mersenne_cliff_table, mersenne_cliff_theorem

# Cliff thresholds for Mersenne numbers
table = mersenne_cliff_table(n_max=8)
for row in table:
    print(f"n={row[0]}: cliff at k*={row[3]}")

# Verify the full theorem
mersenne_cliff_theorem(verbose=True)
```

### p-adic Root Finding

```python
from dyadic_math import lift_root, newton_step, halley_step

# Hensel lift cube root of 2 from mod 7 to mod 7^4
a, p, k = 2, 7, 4
x = lift_root(a, p, k)
print(f"∛2 mod {p}^{k} = {x}")
print(f"Check: {pow(x, 3, p**k)} == {a % (p**k)}")
```

### Fourier Analysis

```python
from dyadic_math import fourier_summary, ultrametric_uncertainty

summary = fourier_summary(k=8, e_true=4)
print(f"Dyadic coefficients: {summary['dyadic_coeffs']}")

uncertainty = ultrametric_uncertainty(k=8)
print(f"Uncertainty bound: {uncertainty}")
```

## Dependencies

- `dyadic-core` (stdlib only)
- `numpy`

## Further Reading

- `docs/mersenne_ghost_theorem.md` — Complete Mersenne Ghost Theorem proof
- `docs/mathematics.md` — Mathematical background
- `docs/padic_newton_dynamics.md` — Higher-order methods, convergence law
- `docs/newton_dynamics/` — p-adic Newton dynamics research
