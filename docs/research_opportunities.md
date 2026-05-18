# Research Opportunities and Open Problems

> This is a dyadic (3-package monorepo) document, ported from the dual-view project.

## Status Dashboard

| Result | Status | Location |
|--------|--------|----------|
| Exponential map isometry T6a | Proved, zero variance | `dyadic_math.isometry` |
| Operator algebra identities T6b | Proved | `dyadic_math.isometry` |
| Trajectory separation T2 | Proved, zero variance | `dyadic_math.separation` |
| Global convergence for α=0 | Proved | `dyadic_math.basin` |
| Ghost fixed point formula T4 | Proved | `dyadic_math.basin` |
| Mersenne Ghost Theorem T5 | Verified n=3..11 | `dyadic_math.mersenne` |
| Bootstrap k/2 optimality | Verified exhaustive | `dyadic_math.mersenne` |
| Ergodic Newton corrections | Empirical (chi-square) | `dyadic_math.padic_roots` |
| Trace-mod-p independence T6c | Empirical (chi-square) | `dyadic_math.isometry` |
| Phase alignment in GL(2) | Empirical (68%, n=30) | `dyadic_math.nonabelian` |
| Commutator depth theorem | Verified | `dyadic_math.iwasawa` |
| p-adic convergence law | Proved, zero variance | `dyadic_math.padic_roots` |
| Popcount compression | Unverified correlation | `dyadic_math.padic_roots` |
| CRT stability correlation | Empirical, sign unexplained | `dyadic_math.crt` |
| Secondary correction ε(n) | Observed, not proved | `dyadic_math.mersenne` |

## Priority 1: Experiments Ready to Run

### 1.1 Confirm Mersenne Cliff at n=16

The secondary correction `ε(n) = v₂(n) - 1` predicts that at `n = 16` (a power of 2), the cliff should be `k* = 16 + 2 + (4 - 1) = 21`. Current verification goes up to `n=12`. Running at `n=16` requires `k=25`, which is feasible with the current LUT-accelerated dlog.

**Implementation**: Extend `mersenne_cliff_table(n_max=16)`. The scan range needs `k` from `n+1` to `n+10`.

### 1.2 QAT / STE Training Experiment

Train a quantized MLP on MNIST with:
- **Control**: Standard STE quantization
- **Experimental**: Ghost-regularized training using `SeedThermodynamics.analyse()` to bias weight updates toward residues with higher `v₂(e_true)`.

The hypothesis is that weights with higher `v₂(e_true)` have higher quantization cliffs and should be more stable under training dynamics.

**Implementation**: Use `dyadic_ml.training` and `dyadic_ml.thermodynamics`.

### 1.3 Trace-mod-p Bridge

The empirical independence of `α(det H)` and `Tr(H) mod p` (T6c) is statistically significant but lacks a theoretical explanation. Two approaches:

1. **Algebraic**: Show that `Tr(H) mod 2` and `det(H) mod 4` (which determines `α`) are structurally separated
2. **Geometric**: Relate to the structure of the `GL(2)` flag variety over F₂

## Priority 2: Theorems to Formalize

### 2.1 Secondary Correction Formula

Prove that for `n = 2^m`:

```
k* = n + 2 + v₂(n) - 1 = n + 2 + m - 1
```

The core identity `5^(2^(n-2)) ≡ 1 - 2^n (mod 2^(n+1))` has an extra term at `n = 2^m`:

```
5^(2^(n-2)) ≡ 1 - 2^n + 2^(2n-4) (mod 2^(n+2))
```

The correction arises because `2n-4 < n+2` for `n < 6`, and the extra term `2^(2n-4)` affects the stability threshold. The general formula should account for all `n = 2^m`.

### 2.2 Bootstrap Optimality Theorem

Prove that `eprec₀ = k/2` is the exact minimizer of the Viglietta dlog bit-cost:

```
cost(eprec₀, k) = eprec₀ + Σ_{i=1}^{⌈log₂((k-2)/eprec₀)⌉} 2ⁱ·eprec₀
```

The derivative `dcost/d(eprec₀) = 0` gives `eprec₀ = (k-2)/2`. The closed-form optimal value is:

```
eprec₀* = (k-2) / 2
```

for the continuous relaxation, and `⌊(k-2)/2⌋` or `⌈(k-2)/2⌉` for the integer case.

### 2.3 Order-m Universality

For any `m`-order root-finding method (Newton m=2, Halley m=3, Newton² m=4, Newton³ m=8), the convergence law is:

```
v_p(x_{n+1} - x*) = m · v_p(x_n - x*)
```

This holds exactly for all `n` where `v_p(x_n - x*) < k`. Prove this by induction using the Taylor expansion of the iteration function.

## Priority 3: Corrected Predictions

### 3.1 Bootstrap Residual Distribution

**Old prediction**: The bootstrap residual `e_boot - e_true` is approximately normal.

**Corrected finding**: The bootstrap residual is **uniformly distributed** over `Z/(2^(k-2))` for random `a`. This follows from the uniformity of the first Newton correction (chi-square verified).

### 3.2 Ghost Density vs. Training

**Old prediction**: Training with ghost regularization should reduce ghost density (fraction of weights in α=1 sector).

**Corrected finding**: Ghost density is always ~50% (random chance of α=0 vs α=1) and is not affected by training. The **graded** quantity is `v₂(e_true)`, which does change under training dynamics.

## Priority 4: Long-Range Questions

### 4.1 2-Adic Julia Set

The Newton map for the dlog function is linearisable (the Julia set is the whole exponent domain). Are there non-linearisable perturbations of the Newton map that produce non-trivial Julia set structure? This connects to the theory of p-adic dynamics and Rivera-Letelier's classification.

### 4.2 Berkovich Line and Tidal Dynamics

The tidal scalar `H = v + α·e` is a point in the Berkovich 2-adic line. The dynamics of weight updates under SGD can be studied as a flow on this line, with the Mersenne cliff being a "boundary" where the flow becomes unstable. This suggests a Berkovich-analytic approach to understanding quantization-aware training.

### 4.3 Halley/Composed-Newton Theorems

For the 2-adic dlog problem, higher-order methods (Halley, composed Newton) have not been studied. The analogous theorems to T1-T6 for order-m methods would include:
- Order-m convergence of the dlog map (generalising T1)
- Order-m trajectory separation: `n*(s) = ceil(log_m(s)) - 1` (generalising T2)
- Order-m ghost formula (generalising T4)

### 4.4 Newton Dynamics Open Problems (from `newton_dynamics/`)

Open problems from the p-adic Newton dynamics research:

1. **Exact period-4 Galois group**: identify in `S₆ wr V₄` using GAP/Magma
2. **Prove clean primes = {7, 103, 181}**: obstruction accumulation argument
3. **Compute period-6 multipliers**: numerical computation from the 233 coefficients
4. **Closed form for individual multipliers**: not just the product formula
5. **Relate K₂ to ray class fields** of `ℚ(ζ₃)`
6. **Generalize to degree d maps**: beyond the cubic case

### 4.5 Frobenius Coin Problem and Cliff Prediction

The Mersenne cliff `k* = n+2` is a diophantine solution to the equation `5^e ≡ 2^n - 1 (mod 2^k)`. Generalising to arbitrary weights gives a Frobenius coin problem: for which `(w, k)` does the Newton iteration have ghost attractors? The answer depends on the 2-adic expansion of `log₅(w)`.

## Methodological Principles

### Mathematical Rigour
- Every theorem must have a proof sketch in the code docstring
- Statistical results (chi-square, ANOVA) must be labelled as "empirical" not "proved"
- Zero-variance results (separation, convergence law) are gold standard

### Empirical Verification
- All theorems are accompanied by verification code that runs in tests
- Randomised tests use seeded RNG for reproducibility
- `n_trials` must be large enough for statistical significance

### Bug Prevention
- Every bug fix is documented with reproducer, impact, and fix rationale
- The `α=1` fix is documented at every level (module docstring, inline comment, bug_history.md)
- Deprecated functions retain the old API with DeprecationWarning for one version
