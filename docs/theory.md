# Unified 2-Adic Operator Theory

This document connects the four pillars of the butterfly ecosystem — **operator
calculus**, **Witt ghosts**, **Mersenne ghosts**, and **thermodynamic
classification** — under a single framework: operations on power series and
integer sequences in the ring ℤ₂[[t]] of 2-adic formal power series.

---

## 1. The Core Object: ℤ₂[[t]]

All four projects operate in the ring of formal power series over the 2-adic
integers:

```
f(t) = Σ_{i=0}^{∞} a_i · t^i,    a_i ∈ ℤ₂
```

A computer represents `f(t)` by truncating both the series (at degree N) and
the coefficients (at word size W). Every precision window, every ghost
identity, and every classification in this ecosystem follows from this
double truncation.

---

## 2. Operator Calculus (dyadic)

Three operators on ℤ₂[[t]] form the foundation.

### 2.1 Derivative D

```
D(f)(t) = df/dt = Σ i · a_i · t^{i-1}
```

Linear, satisfies Leibniz, and is **nilpotent** on degree-N polynomials:
`D^N = 0`.

### 2.2 Forward Difference Δ

```
Δ(f)(t) = f(t+1) - f(t)
```

Discrete analogue of the derivative. On monomials:

```
Δt^n = (t+1)^n - t^n = Σ_{k=0}^{n-1} C(n,k) · t^k
```

### 2.3 The Central Identity: Δ = e^D - I

The forward difference is the exponential of the derivative:

```
Δ = e^D - I = Σ_{k=1}^{∞} D^k / k!
```

This is an exact identity in ℤ₂[[t]] when both sides act on polynomials:
the series terminates at degree N because `D^N = 0`.

**Proof** (closed-form binomial): For a monomial `t^n`:

```
(e^D - I)(t^n) = Σ_{k=1}^{n} C(n,k) · t^{n-k}
               = (t+1)^n - t^n
               = Δ(t^n)
```

The equivalent logarithm identity:

```
D = ln(I + Δ) = Σ_{k=1}^{∞} (-1)^{k+1} Δ^k / k
```

### 2.4 Carry Chain C

The carry chain `C` is the **same operator** as synthetic division by `t-1`:

```
C(a)_i = a_i + C(a)_{i-1}    (carry recurrence)
       ≡ a(t) mod (t-1)^i    (synthetic division interpretation)
```

where `a(t) = Σ a_i · t^i`. The operator triple `{C, D, Δ}` forms a closed
algebraic system:

| Identity | Domain |
|----------|--------|
| `Δ = e^D - I` | All polynomials |
| `D = ln(I + Δ)` | All polynomials |
| `D ∘ Δ = Δ ∘ D` | All polynomials |
| `D^N = 0` | Degree < N |
| `C(a) · C(b) = C(C(a) · b)` | Carry chain |
| `C(a) ≡ a mod (t-1)^∞` | Ghost identity |

### 2.5 Artin-Schreier Operator

The Artin-Schreier operator `℘(x) = x² - x` on ℤ₂ satisfies `℘(0) = ℘(1) = 0`
and `℘(x) = ℘(1-x)` for all x. It classifies the additive group modulo the
Frobenius: `ℤ₂/℘(ℤ₂) ≅ ℤ/2ℤ`.

**Reference**: `include/dyadic/core.h` (artin_schreier),
`include/dyadic/verify.h` (compile-time Artin-Schreier proofs).

---

## 3. Witt Ghosts (dyadic)

### 3.1 Ghost Map

The ghost map sends a Witt vector `a = (a₀, a₁, …, a_{N-1})` to its ghost
components:

```
ghost_j(a) = Σ_{i=0}^{j} 2^i · a_i^{2^{j-i}}
```

The ghost map is a **ring isomorphism** from the Witt ring `W_{N}(ℤ/2^Wℤ)` to
`(ℤ/2^Wℤ)^N`:

```
ghost(a + b)_j = ghost(a)_j + ghost(b)_j
ghost(a · b)_j = ghost(a)_j · ghost(b)_j
```

### 3.2 Newton Recovery

Recovery of Witt components from ghost values uses Newton iteration:

```
r_0 = G_0  (mod 2^W)
r_j = (G_j - Σ_{i<j} 2^i · r_i^{2^{j-i}}) / 2^j  (mod 2^W)
```

This succeeds exactly when `r_j < 2^{W-j}` (the precision window).

### 3.3 Witt Exponential and Logarithm

The Witt exponential applies ℤ₂'s exponential componentwise to ghost values:

```
ghost(exp_W(a))_j = exp(ghost_j(a))
```

where `exp(x) = Σ x^n/n!` converges in ℤ₂ for `v₂(x) ≥ 2`. The Witt logarithm
is the inverse:

```
ghost(log_W(a))_j = log(ghost_j(a))
```

where `log(1+y) = Σ (-1)^{n+1} y^n/n` converges for `v₂(y) ≥ 1`.

The **ghost homomorphism** property of the Witt exponential:

```
exp_W(a + b) = exp_W(a) × exp_W(b)
```

follows directly from `exp(x+y) = exp(x)·exp(y)` on the ghost side.

### 3.4 Connection to Operator Calculus

The ghost-map formula `ghost_j(a) = Σ 2^i · a_i^{2^{j-i}}` is the **same
structure** as the carry chain recurrence: each ghost component accumulates
lower-index contributions multiplied by powers of 2, exactly as synthetic
division accumulates carries.

In both cases, the operation is:

```
result_j = f(j, inputs_{i≤j})
```

with an **invertible triangular system**: given the outputs, the inputs are
recovered by iterative subtraction and division by a power of 2.

**Reference**: `include/dyadic/witt.h` (Witt operations, exp/log),
`include/dyadic/verify.h` (compile-time Witt proofs).

---

## 4. Mersenne Ghosts (dual-view)

### 4.1 Core Identity

For the Mersenne numbers `M_n = 5^{2^{n-2}}`:

```
M_n ≡ 1 - 2^n  (mod 2^{n+1})
```

The exponent `2^{n-2}` is the **ghost sector** boundary: for any
`k < 2^{n-2}`, `5^k` lives in the ghost sector `α = 1` where it behaves
like 1; for `k ≥ 2^{n-2}`, it submerges into the bulk sector `α = 0`.

### 4.2 Four-Part Theorem

1. **Sector**: `α = 1` for all `k ≥ n+1` (ghost sector)
2. **Exponent**: `e_true = 2^{n-2}` for `n+1 ≤ k ≤ n+2` (stable window)
3. **Valuation**: `v₂(e_true) = n - 2`
4. **Cliff**: quantization cliff at `k* = n + 2 + max(0, v₂(n) - 1)`

### 4.3 Connection to Ghost Recovery

The Mersenne ghost phenomenon is the **same precision window** as Witt ghost
recovery. In both cases:

- A "ghost value" (`ghost_j(a)` or `5^k mod 2^W`) is computed at full precision.
- Recovery extracts a component by dividing by a power of 2.
- The division loses bits when the intermediate value exceeds the available window.

The bootstrap formula for optimal Mersenne computation:

```
eprec_i = min(2^i · eprec_0, k - 2),    eprec_0 = k/2
```

mirrors the Witt Newton recovery: each bootstrap step doubles precision until
the target is reached, exactly as each `r_j` recovery step doubles the ghost
index.

### 4.4 Cliff Density

```
Pr[c = 0] = 7/8         87.5% of targets have no cliff
Pr[c = j] = 2^{-(j+3)}  for j ≥ 1
E[c]      = 1/4         expected cliff cost: 0.25 extra bits
```

The cliff is the Mersenne analogue of the precision window overflow in Witt
recovery: most computations stay within the window, and when they don't, the
extra cost is small.

**Reference**: `dual-view/docs/mersenne_ghost_theorem.md` (full proof),
`dual-view/docs/mathematics.md` (background).

---

## 5. Thermodynamic Classification (butterfly-cpp, butterfly-compiler)

### 5.1 Operator Spectrum

Given a linear operator `T` on ℤ₂[[t]] with eigenvalues `λ_i`, classify by
spectral radius `ρ = max |λ_i|`:

| Class | Criterion | Behavior |
|-------|-----------|----------|
| Unitary | `ρ = 1`, `|λ_i| = 1` | Conserves 2-adic norm |
| Contractive | `ρ < 1` | Shrinks perturbations |
| Expansive | `ρ > 1` | Amplifies perturbations |
| Nilpotent | `T^n = 0` | Finite extinction |
| Conservative | `ρ ≤ 1` | Bounded (contractive + unitary) |

### 5.2 Thermodynamic Formalism

For a transfer operator `Ł` acting on a space of 2-adic functions:

```
Lyapunov exponent:   λ(Ł) = lim_{n→∞} (1/n) · log ||Ł^n||
Escape rate:         γ(Ł) = -lim_{n→∞} (1/n) · log |Ł^n(1)|
Topological entropy: h_top = log(ρ(Ł))
```

These quantities are related to the operator calculus by:

```
Ł = e^D      (forward evolution)
Δ = Ł - I    (finite-time difference)
D = ln(Ł)   (infinitesimal generator)
```

The **thermodynamic classification** of `Ł` is equivalent to the spectral
classification of `D` via the exponential map.

### 5.3 Butterfly Seed

The butterfly seed matrix `S` (size `d × d`, eigenvalues `λ_i`) generates the
Kronecker-power operator `S^{⊗n}`. Its thermodynamic classification governs:

- **Spectral butterfly**: `f(S^{⊗n})·v` via diagonalization
- **Quantum circuits**: unitary seeds give quantum gates
- **p-adic Fourier**: contractive seeds give 2-adic convolution kernels

**Reference**: `butterfly-cpp/include/butterfly/thermodynamics.hpp`
(classification + properties), `butterfly-cpp/include/butterfly/spectral.hpp`
(product spectrum + butterfly polynomial).

---

## 6. Unification

### 6.1 The Four-Fold Correspondence

| Domain | Structure | Generator | Recovery | Window |
|--------|-----------|-----------|----------|--------|
| Operator calculus | `{C, D, Δ}` on ℤ₂[[t]] | `D = d/dt` | `C = synthetic_divide` | Degree N |
| Witt ghosts | `W(ℤ/2^Wℤ)` | `ghost_j` sum | Newton iteration | `r_j < 2^{W-j}` |
| Mersenne ghosts | `5^k mod 2^W` | `2^{n-2}` exponent | Bootstrap | `k* = n + 2 + ...` |
| Thermodynamic | `Ł = e^D` | Lyapunov exponent | Spectral radius | Classification |

### 6.2 Unified Principle

Every operation in this ecosystem follows the same pattern:

1. **Embed**: lift components to doubled precision (`gw_t`, `dword_t`).
2. **Compute**: apply the operation in the lifted ring.
3. **Extract**: divide by powers of 2 to recover original precision.

The precision windows (ghost recovery `r_j < 2^{W-j}`, factorial constraint
`j!·FF_j < 2^W`, Mersenne cliff `k*`, bootstrap optimality `eprec_0 = k/2`,
exp truncation `v₂(x) ≥ 2`) are all the same phenomenon: **step 3 fails when
step 1 didn't provide enough headroom**.

### 6.3 The Spectral Sequence

The operator `D` (derivative) acts on ℤ₂[[t]] as a nilpotent operator with
index N. Its spectral decomposition:

```
D = Σ_{i=0}^{N-1} E_i
```

where `E_i` projects onto the `i`-th degree component, corresponds to the
ghost sector decomposition in Mersenne theory. The exponential map
`e^D = I + Δ` sends the additive structure of derivatives to the multiplicative
structure of differences, exactly as `exp_W` sends additive Witt vectors to
multiplicative Witt vectors.

### 6.4 The Six Axioms

The six axioms implemented in dyadic span all four domains:

| # | Axiom | Domain | Implementation |
|---|-------|--------|----------------|
| 1 | Carry chain = synthetic division | Operator | `synthetic_divide` in `core.h` |
| 2 | Formal derivative D | Operator | `formal_derivative` in `core.h` |
| 3 | Forward difference Δ | Operator | `forward_difference` in `core.h` |
| 4 | Exact D | Operator | `exact_formal_derivative` in `core.h` |
| 5 | Exact Δ | Operator | `exact_forward_diff` in `core.h` |
| 6 | Basis change | Unification | `change_basis` in `basis.h` |

---

## 7. Open Problems

1. ~~**Witt exp convergence**~~ **Resolved**: Dynamic term counting (2× bit-width
   budget) achieves `v₂(a₀) ≥ 2` convergence (the theoretical optimum).
   The old fixed 16-term budget required `v₂(ghost_j(a)) ≥ 9`.

2. **Mersenne-Butterfly duality**: Is the Mersenne ghost sector `α = 1`
   isomorphic to the nilpotent sector of the thermodynamic classification?

3. **Operator calculus in W(ℤ/2^Wℤ)**: The carry chain operator `C` acts on
   integer sequences. Does it lift to an operator on Witt vectors that
   commutes with the ghost map?

4. **Thermodynamic Witt exp**: Does `exp_W` have a natural interpretation as
   a transfer operator on a 2-adic function space, with Lyapunov exponents
   determined by the valuations of the ghost inputs?

**Reference**: `dual-view/docs/research_opportunities.md` (full list of
open problems across all projects).

---

## 8. Further Reading

| Project | Document |
|---------|----------|
| dyadic | `README.md`, `docs/precision.md`, `docs/development.md` |
| butterfly-compiler | `README.md` (18 modules, 12 experiments) |
| butterfly-cpp | `thermodynamics.hpp`, `spectral.hpp` |
| dual-view | `docs/mersenne_ghost_theorem.md`, `docs/mathematics.md` |
| BigInteger/carryline | `DRAFT.md`, `validation_paper.md`, `REVIEW.md` |
