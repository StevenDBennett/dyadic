# Mathematical Background

> This is a dyadic (2-package library) document, ported from the dual-view project.

## 2-Adic Numbers

The 2-adic integers Z₂ are the completion of Z with respect to the 2-adic absolute value `|n|₂ = 2^(-v₂(n))`, where `v₂(n)` is the largest exponent such that `2^v₂(n)` divides `n`. Arithmetic in Z₂ is arithmetic modulo powers of 2, lifted to the limit.

### Key fact: Group Structure

For `k ≥ 3`, the multiplicative group of units modulo `2^k` decomposes as:

```
(Z/2^k Z)^×  ≅  Z/2 × Z/2^(k-2)
```

with canonical generator `g = 5`. This means every odd integer modulo `2^k` can be written uniquely as:

```
n ≡ (-1)^α · 5^e  (mod 2^k)
n = 2^v · (-1)^α · 5^e  (mod 2^k)   [full decomposition]
```

where:
- `v = v₂(n)` — 2-adic valuation
- `α ∈ {0, 1}` — the sign sector (0 if `n ≡ 1 mod 4`, 1 if `n ≡ 3 mod 4`)
- `e ∈ [0, 2^(k-2))` — discrete logarithm base 5

This is the **dual-view coordinate system**.

## Core Theorems

### T1: Quadratic Convergence of the Newton dlog Map

The Newton iteration for solving `5^e ≡ a (mod 2^k)` is:

```
e_{n+1} = e_n - (5^{e_n} - a) / (5^{e_n} · log₂(5))
```

where the derivative `log₂(5) = L` is the 2-adic logarithm of 5. This iteration converges with **gain law**:

```
v₂(e_{n+1} - e_true) = 2 · v₂(e_n - e_true) + 1
```

The `+1` bonus comes from dividing the residual by 4 in the correction step. This is strictly better than pure quadratic convergence (which would give `2j`). The proof uses the Lifting the Exponent Lemma together with the Taylor expansion of `(e^x - 1) - x·e^x`, showing that the leading term `-x²/2` has valuation `2j+3`, and dividing by the denominator `4·5^e·L` (valuation 2) yields `2j+1`.

#### Divisor Optimality

The Newton correction divides by `2^d · 5^e · L`. The choice `d = 2` (i.e. `>> 2`) is the unique optimal divisor:

| d | Gain | Behaviour |
|---|------|-----------|
| d < 2 | `v_new = j` | Linear, no gain |
| **d = 2** | **`v_new = 2j+1`** | **Optimal: quadratic + LTE bonus** |
| d > 2 | `v_new = j - (d-2)` | Sublinear; diverges for `d ≥ j+2` |

This is machine-verified for all `d ∈ [0, 5]`, `j ∈ [2, 6]` at `k = 24`. Two half-steps with `d=1` cannot recover the gain of a single step with `d=2`.

### T1b: General 2-Adic Exponential and Logarithm

The 2-adic exponential and logarithm converge on the domains `v₂(x) ≥ 2` and `v₂(g-1) ≥ 2` respectively:

```
exp(x) = Σ xⁿ/n!    (converges for v₂(x) ≥ 2)
log(g) = Σ (-1)^{n+1} (g-1)ⁿ/n    (converges for v₂(g-1) ≥ 2)
```

These are implemented as `padic_exp` and `padic_log` using exact integer arithmetic (tracking `v₂(xⁿ/n!)` analytically to determine termination). The cliff centre `g₀ = exp₂(-4)` is the unique 2-adic unit with `log(g₀) = -4`, and satisfies `v₂(g₀ - (-123)) = 13` — the hardware approximation `-123` agrees with `g₀` to 13 bits.

### T1c: Cliff Density Theory

For a uniformly random log-target `L mod 2^k`, the cliff constant `c(g) = v₂(g - g₀) - 2` follows a geometric distribution:

```
Pr[c = 0]  =  7/8          (no cliff, 87.5% of targets)
Pr[c = j]  =  2^{-(j+3)}   for j ≥ 1
E[c]       =  1/4          (expected cliff cost: 0.25 extra bits)
Pr[c ≥ 4]  ≈  1.56%        (high cliffs are extremely rare)
```

This makes cliff detection worthwhile for correctness but negligible for average-case performance — the expected overhead is only 0.25 bits.

### T1d: Real vs 2-Adic Reconciliation

The 2-adic gain laws correctly predict real Newton convergence rates. For `k`-bit fixed-point arithmetic, `j` bits of 2-adic precision equals `2^{-j}` absolute error in the real sense. The LTE `+1` bonus saves ~17% of iterations at `k=64` (5 steps vs 6 for multiplicative inverse). The 2-adic Newton approach is 10–12× faster than the squaring-loop alternative for `k=64`.

### T2: Trajectory Separation Theorem

Let `a` and `a'` be two targets separated by `v₂(a - a') = s`. Starting from the same seed, their Newton trajectories are identical for exactly:

```
n*(s) = ⌈log₂(s)⌉ - 1
```

steps. This is exact — zero variance — because the Newton map is a 2-adic contraction with factor 1/2. The proof follows from the ultrametric property of the 2-adic absolute value and the Lipschitz constant of the Newton correction map.

### T3: Basin Dichotomy

For targets in the α=0 sector (`a ≡ 1 mod 4`), every seed in the exponent domain `Z/2^(k-2)` converges to the true root under Newton iteration — there are no ghost attractors.

For targets in the α=1 sector (`a ≡ 3 mod 4`), the naive Newton iteration converges to a ghost fixed point `e* = dlog(a+2, k)`, not the true root. The fix is to solve `5^e ≡ -a` instead of `5^e ≡ a` (the α=1 fix).

### T4: Ghost Formula

For a target in the α=1 sector, the ghost fixed point is:

```
e* = dlog(a + 2, k)
```

This gives a false root: `5^{e*} ≡ a+2 ≠ a (mod 2^k)`. The Newton iteration is attracted to this fixed point because the true root lies in a different branch of the 2-adic logarithm.

### T5: Mersenne Cliff Theorem

For Mersenne numbers `w = 2^n - 1`:

- The 2-adic coordinates are always `α = 1` (ghost sector), `e_true = 2^(n-2)`
- `v₂(e_true) = n - 2`
- The quantization cliff occurs at `k* = n + 2`
- At `k = n + 2`, ghost seeds first appear and the weight becomes unstable under Newton iteration

The proof relies on the core identity:

```
5^(2^(n-2)) ≡ 1 - 2^n (mod 2^(n+1))
```

which is verified for `n = 3, ..., 11` by the code.

**Secondary correction**: At `n = 2^m` (powers of 2), there is a secondary correction `ε(n) = v₂(n) - 1` that modifies the cliff formula to `k* = n + 2 + ε(n)`.

### T6a: Exponential Map Isometry

The exponential map `e ↦ 5^e` scales 2-adic distances by exactly factor 4:

```
v₂(5^{e1} - 5^{e2}) = v₂(e1 - e2) + 2   (for e1 ≠ e2)
```

This follows from the 2-adic Taylor expansion `5^e = 1 + 4e + higher terms` and implies the Newton map is a 2-adic contraction.

### T6b: Operator Algebra Identities

On the function space `Z/N → Z/2^k` (where `N = 2^(k-2)`), let:

- `I` = identity operator
- `S` = shift operator: `(Sf)(e) = f(e+1 mod N)`
- `D = I - S` = forward difference
- `avg(f) = Σ_e f(e)` = summation operator

Then:

```
avg² = N · avg
D · avg = 0
avg · D = 0
```

These identities are verified empirically and follow from the telescoping property of finite differences.

### T6d: Mahler Basis and Dirac/Volterra Operators

Every integer-valued function `f: Z → Z` has a unique Mahler (binomial) basis expansion:

```
f(x) = Σ a_n · C(x, n)    where C(x, n) = x(x−1)…(x−n+1)/n!
```

The **Dirac operator** `D` (forward difference) and **Volterra operator** `T` (right inverse) act on Mahler coefficients as:

```
D(a₀, a₁, a₂, …) = (−a₁, −a₂, −a₃, …)     [lowers degree]
T(0, a₁, a₂, …)  = (0, 0, −a₁, −a₂, …)    [raises degree]
```

There is an exact asymmetry at the boundary:

```
D∘T = id  on  ker(ε) = span{e_n : n ≥ 1}
T∘D = id  on  span{e_n : n ≥ 2}   (NOT on all of ker(ε))
```

The asymmetry arises because `D(e₁) = −e₀` exits `ker(ε)`, so `T∘D` is undefined at `e₁`. But `T(e₁) = −e₂` stays in `ker(ε)`, so `D∘T` is well-defined at `e₁`.

### T6c: Trace-mod-p Independence

In GL(2) holonomy over `Z/(2^k · p)Z`, the 2-adic α-sector of the determinant and the trace modulo `p` are statistically independent. This is confirmed by chi-square test at `α = 0.05` significance level.

## The 2-Adic Fourier Uncertainty Principle

The Newton step-count function `h(e)` (number of iterations to converge from seed `e`) depends only on `v₂(e - e_true)`. This means `h` is constant on ultrametric balls, and its discrete Fourier transform has power only at dyadic frequencies:

```
j ∈ {0, N/2, N/4, N/8, ...}
```

The product of the support size and the frequency support size equals `N` — saturating the 2-adic Fourier uncertainty bound.

## Bootstrap Optimality

The Viglietta discrete-logarithm algorithm has two phases:
1. **Bootstrap**: compute `e` to `eprec₀` bits using bit-by-bit iteration
2. **Newton refinement**: double precision each step

Total bits processed:

```
cost(eprec₀, k) = eprec₀ + Σ 2·eprecᵢ  (where eprecᵢ doubles)
```

The optimal bootstrap precision is `eprec₀ ≈ k/2`, not the `√k` heuristic currently used. At `k = 64`, using `k/2` instead of `√k` saves 35.9% of total bits processed.

An 8-bit lookup table (64 entries, 512 bytes) provides O(1) bootstrap to 6-bit precision, eliminating the bootstrap phase entirely for `k ≤ 34`.

## Mersenne Ghost Theorem (Detailed)

See `mersenne_ghost_theorem.md` for the complete proof.

## GL(2) Congruence Filtration

GL(2, Z/2^k) has a natural congruence filtration:

```
Γ(2^j) = {M ∈ GL(2, Z/2^k) : M ≡ I (mod 2^j)}
```

with quotients `Γ(2^j)/Γ(2^(j+1)) ≅ gl(2, F₂)`, a 4-dimensional vector space over F₂. The depth of a matrix `M` is the largest `j` such that `M ≡ I (mod 2^j)` — the matrix analog of `v₂(n-1)`.

The commutator depth theorem states:

```
depth([M, N]) ≥ depth(M) + depth(N)
```

for any `M, N ∈ GL(2, Z/2^k)` with `depth(M) + depth(N) < k`.

## p-adic Root Finding

For solving `x³ ≡ a (mod p^k)` with `p ≠ 2, 3`, multiple methods converge at different rates:

| Method | Order | Step |
|--------|-------|------|
| Newton | 2 | `x ← x - f/f'` |
| Halley | 3 | `x ← x - 2ff' / (2f'² - ff'')` |
| Composed Newton | 4 | Two Newton steps |
| Triple Newton | 8 | Three Newton steps |

The convergence law is exact: `v_p(x_n - x*) = mⁿ · v_p(x₀ - x*)` with zero variance.

## p-adic Newton Dynamics for Cube Roots

For the rational map `N(x) = (2x³+1)/(3x²)` (Newton's method for cube roots), the iterates `N^d(x) = A_d(x)/B_d(x)` satisfy:

```
A_{d+1} = 2·A_d³ + B_d³
B_{d+1} = 3·A_d²·B_d
```

**Dynatomic polynomials** count points of exact period n:

```
Φ_n^*(x) = ∏_{d|n} (N^d(x) - x)^{μ(n/d)}
```

Via the symmetry `N(ζx) = ζN(x)` for ζ³ = 1, all dynatomic polynomials are polynomials in `u = x³`. Degree: `μ₃(n) = Σ_{d|n} μ(n/d)·3^d` (see `docs/newton_dynamics/computation_data.md` for the table).

**Special value formulas** (proven):
- `Φ_n^*(0) = 2^{μ₃(n)/6}` (T7)
- `Φ_n^*(1) = 3^{μ₃(n)/2}` (T8)
- `∏μ = 6^{μ₃(n)/6}` (T9)

See `docs/newton_dynamics/theorems.md` for full proofs.
