# dyadic Documentation Index

## Project Overview

**dyadic** is a 2-adic arithmetic library with Newton dynamics analysis tools.
The core insight is that every odd integer modulo `2^k` decomposes uniquely as:

```
n = 2^v · (-1)^α · 5^e   (mod 2^k)
```

where:
- `v = v₂(n)` — the 2-adic valuation
- `α ∈ {0, 1}` — sign sector (0 for `n ≡ 1 mod 4`, 1 for `n ≡ 3 mod 4`)
- `e ∈ [0, 2^(k-2))` — discrete logarithm base 5

This coordinate system reveals quantization cliffs: specific bit-precisions
where Newton iteration becomes unstable ("ghost attractors").

## Packages

| Package | Description | Key Modules |
|---------|-------------|-------------|
| **dyadic-core** | 2-adic arithmetic (stdlib only) | `core`, `exponent`, `mahler` |
| **dyadic-math** | Newton dynamics for the discrete-log map | `basin`, `mersenne`, `separation`, `fourier`, `thermodynamics`, `padic_roots`, `isometry`, `crt`, `nonabelian`, `iwasawa` |

## Module Reference

| Module | Lines | Description |
|--------|-------|-------------|
| `dyadic_core.core` | 684 | DualNumber, modular inverse, 2-adic exp/log, cliff centre g₀ |
| `dyadic_core.exponent` | 81 | ExponentSpace — additive coordinate chart on Z/2^(k-2) |
| `dyadic_core.mahler` | 123 | Mahler basis, Dirac/Volterra operators, boundary asymmetry |
| `dyadic_math.basin` | 361 | BasinExplorer, ghost detection, Newton basin portraits |
| `dyadic_math.mersenne` | 585 | Mersenne Ghost Theorem, bootstrap optimality, cliff constant proofs |
| `dyadic_math.separation` | 185 | Trajectory Separation Theorem |
| `dyadic_math.fourier` | 185 | DFT of Newton step-count function |
| `dyadic_math.padic_roots` | 311 | Multi-order p-adic root finding |
| `dyadic_math.isometry` | 257 | Exponential isometry, operator algebra theorems |
| `dyadic_math.weight_stability` | 289 | WeightStabilityDiagnostics — graded 2-adic stability diagnostics |
| `dyadic_math.crt` | 205 | CRT extension to Z/(2^k·p)Z |
| `dyadic_math.nonabelian` | 214 | GL(2) matrix utilities with phase alignment |
| `dyadic_math.iwasawa` | 235 | GL(2) congruence filtration, LDU decomposition |

## Theorem Reference

| ID | Theorem | Module | Status |
|----|---------|--------|--------|
| T1 | Gain law: `v(e_new - e_true) = 2j+1` (quadratic + LTE bonus) | `core` | Proven |
| T1b | General 2-adic exp/log via exact integer arithmetic | `core` | Verified |
| T1c | Cliff density: `Pr[c=0]=7/8`, `E[c]=1/4` | `core` | Proven |
| T2 | Trajectory separation: `n*(s) = ceil(log₂(s)) - 1` | `separation` | Proven, zero variance |
| T3 | Basin dichotomy: α=0 globally stable | `basin` | Proven |
| T4 | Ghost formula: `e* = dlog(a+2, k)` for α=1 | `basin` | Proven |
| T5 | Mersenne cliff: `k* = n+2` for `w = 2^n - 1` | `mersenne` | Verified n=3..11 |
| T6a | Exponential map isometry: `v₂(5^e-1) = v₂(e)+2` | `isometry` | Proven |
| T6b | Operator algebra: `avg² = N·avg`, `D·avg = avg·D = 0` | `isometry` | Proven |
| T6d | Mahler basis: boundary asymmetry of D∘T / T∘D | `mahler` | Proven |
| T6c | Trace-mod-p independence (GL(2) holonomy) | `isometry` | Statistical |
| — | Commutator depth: `depth([M,N]) >= depth(M)+depth(N)` | `iwasawa` | Verified |
| — | Mersenne cliff constant: `c(g) = v₂(g - exp₂(-4)) - 2` | `mersenne` | Proven |
| — | p-adic convergence law: `v_p(x_n-x*) = m^n·v_p(x₀-x*)` | `padic_roots` | Proven, zero variance |

## Documents

| Document | Description |
|----------|-------------|
| `tutorial.md` | Getting-started walkthrough |
| `api.md` | Complete API reference |
| `mathematics.md` | Full mathematical background |
| `mersenne_ghost_theorem.md` | Mersenne Ghost Theorem proof |
| `research_opportunities.md` | Open problems, experiments, future directions |
| `CHANGELOG.md` | Version history |
