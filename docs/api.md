# API Reference

> This is a dyadic (2-package library) document, ported from the dual-view project.

The dyadic library is split into two packages:

- **dyadic-core** (`dyadic_core`): 2-adic arithmetic foundation — DualNumber, TwoAdicProcessor, modinv_newton, padic_exp/log, MahlerCalculus, ExponentSpace. Stdlib only.
- **dyadic-math** (`dyadic_math`): p-adic Newton dynamics, CRT extensions, Fourier analysis, Mersenne theory, thermodynamics. Depends on `dyadic-core` and `numpy`.

## dyadic-core

### `dyadic_core.core` Functions

#### `_mask(k: int) -> int`
Bit mask `(1 << k) - 1`.

#### `_valuation(n: int) -> int`
2-adic valuation `v₂(n)`. Returns `float('inf')` for `n = 0`.

#### `modinv_newton(a: int, k: int) -> int`
`a⁻¹ mod 2^k` via quadratic Newton lifting. Requires `a` odd.

#### `two_adic_log5(k: int) -> int`
2-adic logarithm of 5 truncated to `k` bits. The unique `L ∈ Z₂` satisfying `exp₂(L) = 5`. Cached per `k` (LRU, 256 entries).

Used as derivative scaling factor in Newton dlog: the derivative of `e ↦ 5^e` is `5^e · (L >> 2)`.

#### `two_adic_dlog(a: int, k: int, L: Optional[int] = None) -> Optional[Tuple[int, int]]`
Full 2-adic decomposition of odd part: `a = (-1)^α · 5^e (mod 2^k)`.
Returns `(alpha, e)` or `None` if `a` is even. Uses 10-bit LUT bootstrap for `k ≤ 34` (256 entries, optimal sweet spot), bit-by-bit to `k/2` otherwise.

#### `dual_add(v_a: int, alpha_a: int, e_a: int, v_b: int, alpha_b: int, e_b: int, k: int) -> Tuple[int | float, int, int]`
Exact addition in dual `(v, α, e)` coordinates via LTE, without converting back to group representation. Given `a = 2^{v_a}·(-1)^{α_a}·5^{e_a}` and `b = 2^{v_b}·(-1)^{α_b}·5^{e_b}`, returns `(v_sum, α_sum, e_sum)` for `a + b mod 2^k`. Handles 5 cases: different valuations (annihilation), same sign (doubling or v+1), opposite signs (v+2+v₂(Δe) via LTE), and exact cancellation.

#### `run_all_tests(k: int = 16, verbose: bool = True) -> None`
Self-check for core arithmetic: round-trip, multiplication, inversion, powering.

#### `padic_exp(x: int, k: int) -> int`
General 2-adic exponential `exp(x) mod 2^k`. Requires `v₂(x) ≥ 2`. Uses exact integer arithmetic with analytic valuation tracking for termination.

#### `padic_log(g: int, k: int) -> int`
General 2-adic logarithm `log(g) mod 2^k`. Requires `g ≡ 1 mod 4` (i.e. `v₂(g−1) ≥ 2`). Uses exact integer arithmetic.

#### `g0(k: int) -> int`
The cliff centre `g₀ = exp₂(−4) mod 2^k`. The unique 2-adic unit with `log(g₀) = −4`. The hardware approximation `−123` agrees with `g₀` to 13 bits.

#### `dlog_residual_tracking(a: int, k: int, L: Optional[int] = None) -> Tuple[int, List[Dict]]`
Viglietta discrete log with normalised residual tracking at each Newton step. Returns `(e, history)` where each entry in history contains `tau_before`/`tau_after` (the residual `(5^e · a⁻¹ − 1) / 4`) and its 2-adic valuation, confirming the quadratic convergence gain law. Requires `a ≡ 1 (mod 4)`.

### `dyadic_core.core` Classes

#### `DualNumber(n: int, k: int = 64)`
2-adic dual-view decomposition of `n` modulo `2^k`.

**Attributes**: `v` (valuation, float('inf') for zero), `alpha` (0 or 1), `e` (exponent), `is_zero`, `value` (the integer `n mod 2^k`).

**Methods**:
- `verify()` — round-trip check: coordinates → integer matches stored integer
- `coords()` — returns `(v, alpha, e)` tuple
- `from_coords(v, alpha, e, k)` — classmethod, build from coordinates

#### `TwoAdicProcessor(k: int = 64)`
Arithmetic on DualNumbers in coordinate space.

**Methods**:
- `mul(a, b)` — multiply, coordinates add componentwise
- `inv(a)` — invert unit (v=0)
- `pow(a, n)` — integer power (negative → inverse)
- `dlog(a)` — returns `(v, alpha, e)`

### `dyadic_core.exponent` Module

#### `ExponentSpace(g: int, k: int)`
Model the multiplicative group `⟨g⟩ ≅ Z/N` where `N = 2^(k-2)`.

**Methods**:
- `lift(e)` — `g^e mod 2^k`
- `difference(f, e)` — forward difference `(Df)(e) = f(e+1) - f(e)`
- `integrate(f)` — discrete integral `Σ f(e) mod 2^k`
- `is_eigenfunction(f, e)` — check `D(g^e) = (g-1)·g^e`

### `dyadic_core.mahler` Module

#### `MahlerCalculus`
Operations on the Mahler (binomial) basis for integer-valued functions.

**Static methods**:
- `mahler_polynomial(n, x)` — binomial coefficient `C(x, n)`
- `to_mahler(f, max_degree)` — convert function values to Mahler coefficients via finite-difference table
- `from_mahler(coeffs, x)` — evaluate function from Mahler coefficients at point `x`
- `dirac_operator(coeffs)` — forward difference on coefficients: `D(a₀, a₁, …) = (−a₁, −a₂, …)`
- `volterra_operator(coeffs)` — right inverse of Dirac: `T(0, a₁, …) = (0, 0, −a₁, …)`; requires `a₀ = 0`
- `truncate(coeffs, k)` — reduce coefficients modulo `2^k`

**Boundary behaviour**:
```
D∘T = id  on  ker(ε) = span{e_n : n ≥ 1}
T∘D = id  on  span{e_n : n ≥ 2}   (NOT on all of ker(ε))
```

## dyadic-math


### `dyadic_math.basin` Module

#### `BasinExplorer(k: int, g: int, target_a: int)`
Newton basin landscape analysis.

**Methods**:
- `newton_step(e)` — single Newton iteration
- `classify(e0)` — returns `(fate, value, path)`, fate ∈ {'converged', 'cycle', 'diverged'}
- `portrait()` — full basin portrait dict over all seeds
- `fate_vector()` — compact encoding (0=converged, 1=cycle, 2=diverged)

#### `LayerGhostDiagnosticV2(k: int, g: int = 5, max_iter: int = 100)`
Per-layer diagnostic for weight matrices.

**Methods**:
- `diagnostic_matrix(W)` — returns `(fate, conv_ratio, ghost_ratio, mean_e, v2_e)`

#### `GhostHunt(g: int = 5, max_iter: int = 100)`
Precision threshold hunting.

**Methods**:
- `precision_threshold_sweep(k_min, k_max, target_e)` — print sweep table
- `quantization_cliff(W, k_min, k_max)` — ghost density vs k for a matrix

#### `precision_sweep(k_min, k_max, g=5, target_e=2) -> List[Tuple[int, float]]`
Standalone precision sweep.

### `dyadic_math.crt` Module

#### `CRTDualNumber(n: int, k: int, p: int, g_p: int)`
Element of `Z/(2^k·p)Z` with dual coordinates.

#### `CRTDualProcessor(k: int, p: int, g_p: Optional[int] = None)`
CRT arithmetic.

**Methods**:
- `crt_reconstruct(r2, rp)` — CRT combine residues
- `mul(A, B)` — multiply two CRTDualNumbers
- `product(weights)` — product of raw integers
- `cycle_product(numbers)` — product of CRTDualNumber list
- `convergence_ratio_2adic(P)` — 2-adic convergence ratio

#### `combined_stability(k, p, num_cycles=50, cycle_length=4) -> Dict`
Randomised correlation test: correlates `v₂(prod)` with the change in `v₂` under 2-adic multiplicative perturbation `w → w·(1+2^t)`.

### `dyadic_math.nonabelian` Module

#### `NonAbelianCRTDual(k: int, p: int)`
GL(2) gauge theory on a cycle.

**Methods**:
- `holonomy(mats)` — product of all matrices
- `invariants(mats)` — dict with det_mod2k, alpha_det, trace_modp, crt views
- `convergence_ratio_full(mats)` — ghost ratio of determinant

#### `phase_alignment_experiment(k, p, N_cycle=4, n_cycles=30) -> Dict`
Test α-sector flip under perturbation.

### `dyadic_math.separation` Module

- `newton_trajectory(a, k, e_seed, steps=10)` — per-step Newton history
- `separation_step(a, a_prime, k, e_seed)` — first divergence step
- `predicted_separation(s, method_order=2)` — theoretical `n*(s)`
- `verify_separation(k, s_values, n_trials=50)` — zero-variance verification
- `ultrametric_ball_tree(k, e_true, depth=3)` — ASCII tree
- `step_count_profile(k, e_true)` — v₂ level counts

### `dyadic_math.fourier` Module

- `step_count_fn(k, e_true)` — numeric step-count function
- `analytic_step_count(k, e_true)` — closed-form O(N) construction
- `dft(f)` — numpy FFT wrapper
- `power_spectrum(f)` — |DFT|²
- `dyadic_coefficients(f)` — extract at N/2, N/4, ...
- `analytic_coefficients(k)` — closed-form Fourier coefficients
- `fourier_summary(k, e_true)` — complete analysis
- `ultrametric_uncertainty(k)` — uncertainty principle statement

### `dyadic_math.padic_roots` Module

- `lift_root(a, p, k)` — Hensel lift cube root from mod p to mod p^k
- `newton_step(x, a, pk)` — order 2
- `halley_step(x, a, pk)` — order 3
- `newton2_step(x, a, pk)` — order 4 (composed Newton)
- `newton3_step(x, a, pk)` — order 8 (triple Newton)
- `convergence_profile(x0, a, p, k, step_fn, x_true)` — track v_p
- `compare_methods(p, k, n_trials)` — rate comparison
- `verify_order(primes, k, n_trials)` — verify convergence ratio
- `newton_correction_uniformity(p, k, n_seeds)` — chi-square test
- `popcount_compression(k, n_trials)` — popcount correlation

### `dyadic_math.iwasawa` Module

- `congruence_depth(M, k)` — depth in GL(2) filtration
- `filtration_residue(M, depth, k)` — gl(2, F₂) direction
- `ldu_decompose(M, k)` — LDU factorisation
- `matrix_coordinates(M, k)` — full coordinate decomposition
- `holonomy_depth_profile(k, p, ...)` — depth under perturbation
- `filtration_portrait(k)` — GL(2) quotient sizes
- `matrix_commutator(M, N, k)` — `[M,N]` computation
- `verify_commutator_depth(k, ...)` — depth theorem verification
- `MatrixCoordinates` — dataclass

### `dyadic_math.mersenne` Module

- `mersenne_coordinates(n, k)` — `(α, e_true, v₂(e_true))`
- `verify_core_identity(n_max)` — 5^(2^(n-2)) ≡ 1 - 2^n
- `mersenne_cliff_table(n_max)` — cliff thresholds
- `bootstrap_cost(eprec0, k)` — Viglietta bit-cost
- `optimal_bootstrap(k_values)` — minimiser search
- `compare_bootstrap_strategies(k_values)` — sqrt vs k/2 vs LUT
- `dlog_with_lut(a, k, b=8)` — LUT-based dlog
- `verify_lut_dlog(k, b=8, n_trials)` — correctness check
- `cliff_constant(g, k)` — compute `c = v₂(log₂(g)/4 + 1)`
- `cliff_formula(g)` — human-readable c(g) formula
- `mersenne_cliff_theorem(verbose)` — state and verify the full theorem
- `prove_cliff_constant(verbose)` — prove `c=5` from 4 log-series terms
- `prove_c_formula(verbose)` — prove `c(g) = v₂(g-5) - 2`
- `exp2_neg4(k)` — compute `exp₂(-4) mod 2^k`, the zero of `log₂(g)/4+1`
- `cliff_constant_unified(g, k)` — unified formula via Newton-Taylor lemma
- `verify_unified_formula(g_values, k)` — verify unified matches direct
- `proof_connection(verbose)` — show all proofs connected via `log₂(5) ≡ -4 (mod 128)`

### `dyadic_math.isometry` Module

- `verify_isometry(k, n_trials)` — v₂(5^e-1) = v₂(e)+2
- `isometry_pair_test(k, n_trials)` — pair form
- `isometry_summary(k)` — full conditioning picture
- `trace_alpha_independence(k, p, ...)` — chi-square test
- `trace_exponent_independence(k, p, ...)` — ANOVA F-test
- `exponent_valuation_profile(k, n_samples)` — v₂(e_true) distribution

### `dyadic_math.thermodynamics` Module

#### `SeedThermodynamics(k: int = 16, g: int = 5)`
Graded 2-adic weight stability diagnostics.

**Coordinate analysis**:
- `weight_coordinates(w)` — return `(v, α, e)`
- `analyse(W)` — return dict with alpha_fraction, mean_v2_e, std_v2_e, cliff_risk
- `mersenne_cliff_score(w)` — returns `k* = v₂(e_true) + 2`
- `compare_to_random(W, n_samples)` — z-score comparison

**Precision-sweep analysis**:
- `__call__(weights, k_range)` — configure for sweep
- `compute()` — lazy computation of ghost profiles
- `profiles` — property: per-weight ghost density profiles
- `cliffs` — property: per-weight cliff thresholds
- `stable_weights(max_k=None)` — indices of stable weights
- `ghost_weights(max_k=None)` — indices of ghost-affected weights
- `cliff_histogram()` — cliff distribution
- `summary()` — aggregate statistics
- `report()` — formatted ASCII report with bar chart
