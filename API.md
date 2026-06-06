# dyadic — 2-Adic Operator Calculus API Reference

## Overview

dyadic is a **header-only C++20** library implementing the six-axiom 2-adic operator
calculus. All operations are `constexpr`, support compile-time verification, and
work on any C++20 compiler with zero external dependencies.

All symbols live in `namespace dyadic`. The umbrella header is `dyadic.h`.

---

## 1. `core.h` — Core Types & Primitives

### Word Types

| Type | Description |
|---|---|
| `Ring<W, Tag>` | Typed word with a ring-theoretic tag (`Exact`, `Modular<W>`, `Unnormalized`) |
| `ExactRing_t<W>` | `Ring<W, Exact>` |
| `NormalizedRing_t<W>` | `Ring<W, Modular<W>>` |
| `UnnormalizedRing_t<W>` | `Ring<W, Unnormalized>` |
| `dword_t<W>` | Double-width type (8→16, 16→32, 32→64, 64→`detail::uint128_t`) |
| `quad_width<W>` | Quad-width (`widen_t<dword_t<W>>`), used for unsaturated accumulation |

### `detail::uint128_t`

Software 128-bit unsigned pair (`{lo, hi}`). Full constexpr arithmetic:
`+`, `-`, `*`, `<<`, `>>`, `|`, `&`, `^`, `~`, `/`, `%`, `+=`, `-=`, `*=`, `<<=`, `>>=`, `%=`,
comparisons, `explicit operator unsigned __int128()` (where available).

Division and modulo use bit-by-bit long division (slow but constexpr-safe).

### 2-Adic Primitives

| Function | Signature | Description |
|---|---|---|
| `v2(x)` | `W -> int` | 2-adic valuation: `max{k: 2^k \| x}`, returns bit width for 0 |
| `modinv_odd(a)` | `W -> W` | Modular inverse of odd `a` modulo 2^W via Newton iteration |
| `exact_divide(x, d)` | `W -> W` | `x / d` when `d` is odd (returns `x * modinv_odd(d)`) |
| `div_2k_adic(val, k)` | `U -> U` | Signed 2-adic right-shift by `k`: sign-extends when high bit set |
| `artin_schreier(x)` | `W -> W` | `x² - x`. Artin-Schreier operator: kernel = {0, 1}. Symmetric: `℘(x) = ℘(1-x)` |

### Basis Tags

- `MonomialBasis` — standard power basis: `{a₀, a₁, a₂, …}`
- `FallingFactorialBasis` — falling factorial basis: `{(t)₀, (t)₁, (t)₂, …}`
- `TaylorBasis` — divided-power basis: `{T₀, T₁, T₂, …}` where `T_k = k!·FF_k`

### Combinatorial Caches (`namespace detail`)

- `StirlingCache<N, W>` — precomputes Stirling numbers S₂(n,k), s₁(n,k) for n,k < N
- `PascalCache<N, W>` — precomputes binomial coefficients C(n,k) for n,k < N

Both are `inline constexpr` globals (`STIRLING_CACHE`, `PASCAL_CACHE`).

### `Polynomial<N, W, Basis>`

Core polynomial type. Inherits from `std::array<W, N>`.

```cpp
template<int N, std::unsigned_integral W, typename Basis = MonomialBasis>
struct Polynomial : std::array<W, N> { … };
```

**Members:**
- `word_type` — aliases `W`
- `basis_type` — aliases `Basis`
- `static constexpr int max_degree = N - 1`
- `constexpr int actual_degree()` — highest index with non-zero coefficient, or -1
- `eval(x)` — Horner evaluation (**monomial only**; FF/Taylor via `dyadic/eval.h` or `change_basis`)
- `operator+(…)`, `operator-(…)` — coefficient-wise addition/subtraction
- `operator*(…)` — **carry-chain multiplication** (big-integer semantic); see `poly_mul_cw` for standard ring

**Arithmetic semantics:**
- `+`, `-` are always coefficient-wise
- `*` is the **carry-chain** (big-integer) ring — use `poly_mul_cw` for standard convolution

### `synthetic_divide` (Axiom 1)

```cpp
Ring<W, Modular<W>> synthetic_divide(Ring<W, Unnormalized>);
Polynomial<N, W, Basis> synthetic_divide(Polynomial<N, W, Basis>);
```

Implements carry propagation as `(I − N)⁻¹`. On polynomials, splits each word at
its midpoint and carries the top half into the next limb. Also called
"carry-normalize the middle."

### `carry_normalize`

```cpp
Ring<W, Modular<W>> carry_normalize(Ring<W, Unnormalized>);
Ring<W, Modular<W>> carry_normalize(Ring<W, Modular<W>>);
Ring<W, Modular<W>> carry_normalize(Ring<W, Exact>);
```

Converts any tagged ring to its normalized form.

### `carry_chain`

```cpp
template<std::unsigned_integral W, typename Accum = dword_t<W>>
constexpr Accum carry_chain(W* r, const Accum* v, int N);
```

Full-width carry propagation: `r[i] = low(v[i] + carry)`, `carry = high(v[i] + carry)`.
Returns the final carry. Accumulator type defaults to `dword_t<W>`.

```cpp
constexpr W carry_chain_word(W hi, W lo);
```
Returns `low(hi·2^W + lo)` — the low word of a double-width pair.

### `poly_mul`

```cpp
void poly_mul(W* r, const W* a, int na, const W* b, int nb);
```

Carry-chain polynomial multiplication. Uses `quad_width<W>` accumulators and
processes in 256-bit chunks to keep stack small. For short products (≤ 4 quad words)
uses a single-buffer path.

### Standard Ring Multiplication (coefficient-wise)

```cpp
Polynomial<N+M-1, W, Basis> poly_mul_cw(const Polynomial<N,W,Basis>& a,
                                         const Polynomial<M,W,Basis>& b);
Polynomial<N+M-1, W, Basis> operator*(const Polynomial<N,W,Basis>& a,
                                       const Polynomial<M,W,Basis>& b);
```

- `poly_mul_cw` — standard polynomial convolution (the usual coefficient ring)
- `operator*` — **carry-chain** (big-integer) multiplication (different from `poly_mul_cw`!)

### Division Verification

```cpp
bool verify_divmod_cw(const Polynomial<N,W,MonomialBasis>& A,
                      const Polynomial<N,W,MonomialBasis>& Q,
                      const Polynomial<M,W,MonomialBasis>& B,
                      const Polynomial<M,W,MonomialBasis>& R);
```

Verifies `A = Q·B + R` in the coefficient-wise ring.

### Formal Derivative (Axiom 2)

```cpp
Polynomial<N-1, W, MonomialBasis> formal_derivative(const Polynomial<N,W,MonomialBasis>& p);
```

Monomial basis: `D(t^n) = n·t^{n-1}`. FF/Taylor overloads require `dyadic/basis.h`.

### Forward Difference (Axiom 3)

```cpp
Polynomial<N-1, W, MonomialBasis> forward_difference(const Polynomial<N,W,MonomialBasis>& p);
```

Monomial basis: `(Δp)(t) = p(t+1) - p(t)`. FF/Taylor overloads require `dyadic/basis.h`.

### Exact Variants (Axioms 4 & 5)

```cpp
Polynomial<N-1, W, MonomialBasis> exact_forward_difference(…);  // same as Δ
Polynomial<N-1, W, MonomialBasis> exact_formal_derivative(…);   // same as D
```

In the monomial basis these are identical to the non-exact versions. Distinction
matters in FF and Taylor bases.

### Internal Helpers (`namespace detail`)

| Function | Description |
|---|---|
| `adc(a, b, carry_in)` | Add with carry: returns `(sum, carry_out)` |
| `add_overflow(a, b, result)` | Overflow-checked addition (wraps `__builtin_add_overflow`) |
| `mul_overflow(a, b, result)` | Overflow-checked multiplication (wraps `__builtin_mul_overflow`) |
| `poly_mul_unsaturated(r, a, na, b, nb)` | Unsaturated multiply into accumulator array |

---

## 2. `basis.h` — Basis Conversion (Axiom 6)

### `change_basis`

```cpp
template<typename ToBasis, int N, std::unsigned_integral W, typename FromBasis>
constexpr Polynomial<N, W, ToBasis> change_basis(const Polynomial<N, W, FromBasis>& p);
```

Converts between bases:
- `MonomialBasis` ↔ `FallingFactorialBasis` (via Stirling numbers)
- `MonomialBasis` ↔ `TaylorBasis` (via Stirling numbers + factorials)
- FF ↔ Taylor: goes through monomial

### Low-level conversions

| Function | From | To |
|---|---|---|
| `monomial_to_falling(p)` | Monomial | FallingFactorial |
| `falling_to_monomial(p)` | FallingFactorial | Monomial |
| `monomial_to_taylor(p)` | Monomial | Taylor |
| `taylor_to_monomial(p)` | Taylor | Monomial |

### Eval in FF/Taylor bases

```cpp
W eval(const Polynomial<N, W, FallingFactorialBasis>& p, W x);
W eval(const Polynomial<N, W, TaylorBasis>& p, W x);
```

Evaluates polynomial in falling factorial or Taylor basis (uses 2-adic division
for factorial denominators).

### Formal Derivative — FF/Taylor bases

```cpp
Polynomial<N-1, W, FallingFactorialBasis> formal_derivative(…);
Polynomial<N-1, W, TaylorBasis> formal_derivative(…);
```

Converts to monomial, differentiates, converts back.

### Forward Difference — FF/Taylor bases

```cpp
Polynomial<N-1, W, FallingFactorialBasis> forward_difference(…);
Polynomial<N-1, W, TaylorBasis> forward_difference(…);
```

In FF basis: `Δᵢ = i·p[i]` (simple shift). In Taylor/monomial: converts.

---

## 3. `calculus.h` — Taylor Shift, Indefinite Sum, Exp/Log

### `taylor_shift`

```cpp
Polynomial<N, W, Basis> taylor_shift(const Polynomial<N, W, Basis>& p, W delta);
```

Computes `P(t + delta)` using closed-form binomial transform. Works in all three
bases:
- **Monomial**: uses Pascal's triangle and `δⁿ⁻ᵏ`
- **FF**: uses falling factorial of delta
- **Taylor**: converts to monomial, shifts, converts back

### `indefinite_sum`

```cpp
Polynomial<N+1, W, Basis> indefinite_sum(const Polynomial<N, W, Basis>& p);
```

Computes Σ (inverse of Δ). In FF basis: `Σ(t_{n−1}) = (t)_n / n`.
Uses `modinv_odd` or `div_2k_adic` + `modinv_odd` for the division by `n`.

### Power Series Exponential — `poly_exp`

```cpp
template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis>
poly_exp(const Polynomial<N, W, MonomialBasis>& P);
```

Computes `exp(P) mod x^N` using recurrence:
```
E₀ = 1
Eₙ = (1/n) · Σ_{k=1}^{n} k·Pₖ·Eₙ₋ₖ
```

**Precondition:** `P[0] = 0`, `P[1]` even (valuation ≥ 1 in ℤ₂).
Uses `dword_t` intermediates with 2-adic division by `n`.

### Power Series Logarithm — `poly_log`

```cpp
template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis>
poly_log(const Polynomial<N, W, MonomialBasis>& P);
```

Computes `log(1+P) mod x^N` using recurrence:
```
L₀ = 0
Lₙ = Pₙ − (1/n) · Σ_{k=1}^{n−1} k·Lₖ·Pₙ₋ₖ
```

**Precondition:** `P[0] = 0`, `P[1]` even (valuation ≥ 1 in ℤ₂).
**Roundtrip:** `poly_log` and `poly_exp` compose: `exp(log(1+P)) = 1+P`.

---

## 4. `arith.h` — Unsigned Big-Integer Arithmetic on Polynomial

This header treats `Polynomial<N, W, MonomialBasis>` as an N-limb unsigned big
integer in little-endian order (limb 0 = least significant). All operations are
in `namespace dyadic` (public) or `namespace dyadic::detail` (internal).

### Internal Helpers (`namespace detail`)

| Function | Description |
|---|---|
| `zero_extend<TargetN>(p)` | Extend polynomial to `TargetN ≥ N` limbs, zero-filled |
| `unsigned_lt(a, b)` | Lexicographic less-than (big-integer compare) |
| `shift_left(p, shift)` | Left shift by `shift` bits (with cross-limb carry) |
| `clz_normalize(x)` | Count leading zeros of a single word |

### `mul_unsigned`

```cpp
template<int N, int M, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N + M, W, Basis>
mul_unsigned(const Polynomial<N, W, Basis>& a, const Polynomial<M, W, Basis>& b);
```

Full-width unsigned multiplication returning `N+M` limbs. Uses `quad_width<W>`
unsaturated convolution, then a full carry-chain that can produce a carry into
the highest limb.

### `divmod_single`

```cpp
template<int NL, std::unsigned_integral W, typename Basis>
constexpr std::pair<Polynomial<NL, W, Basis>, W>
divmod_single(const Polynomial<NL, W, Basis>& u, W v);
```

Single-limb divisor: divides a multi-limb dividend `u` by a single word `v` using
schoolbook long division (top-to-bottom). Returns `(quotient, remainder)`.

### `div_unsigned_knuth`

```cpp
template<int NL, int NR, std::unsigned_integral W>
constexpr std::pair<…> div_unsigned_knuth(const Polynomial<NL,W,MonomialBasis>& u,
                                          const Polynomial<NR,W,MonomialBasis>& v);
```

Knuth Algorithm D for unsigned division. Precondition: `v[NR-1] ≠ 0`, `NR ≥ 2`.
Normalizes divisor so its top bit is set, then performs trial division with
correction (the standard multi-limb long-division algorithm).

### `CarrySavePair`

```cpp
template<int N, std::unsigned_integral W, typename Basis>
struct CarrySavePair {
    Polynomial<N, W, Basis> sum;
    Polynomial<N + 1, W, Basis> carry;
};
```

Result of carry-save addition: `value = sum + carry`. The `carry` vector is
offset by one limb (`carry[0] = 0`) so each `carry[i]` is the carry-out from
limb `i−1`.

### Carry-Save Operations

| Function | Description |
|---|---|
| `carry_save_add(a, b)` | Decompose `a+b` into `(sum, carry)`; `carry[i] ∈ {0,1}` |
| `csa3(a, b, c)` | 3-2 carry-save adder: `a+b+c = sum + carry`; `carry[i] ∈ {0,1,2}` |
| `combine_carry_save(sum, carry)` | Propagate carries to get a normal integer (final overflow discarded) |
| `batched_add_carry_save(ops, count)` | Sum multiple operands with single final carry propagation; `count=0` returns zero, `count=1` returns a copy |

### `div_unsigned`

```cpp
template<int NL, int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<NL, W, MonomialBasis>, Polynomial<NR, W, MonomialBasis>>
div_unsigned(const Polynomial<NL, W, MonomialBasis>& u,
             const Polynomial<NR, W, MonomialBasis>& v);
```

Polymorphic unsigned division wrapper:
- Strips leading zero limbs from `v` (determines `eff_nr`)
- `eff_nr = 0` → returns `{}` (zero-divisor sentinel)
- `eff_nr = 1` → delegates to `divmod_single`
- `eff_nr < NR` → recurses with shortened divisor
- Otherwise → calls `div_unsigned_knuth`

### Newton Division

#### `reciprocal_newton`

```cpp
template<int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<2*NR, W, MonomialBasis>, int>
reciprocal_newton(const Polynomial<NR, W, MonomialBasis>& d);
```

Computes `r ≈ B^{2·NR} / d` using Newton iteration:
`r_{k+1} = r_k · (2·B^{2n} − d·r_k) / B^{2n}`.

Returns `(reciprocal, norm_shift)` where the reciprocal has `2·NR` limbs and
`norm_shift` records the normalization shift applied to `d`.

**Algorithm details:**
1. Normalizes `d` so `d[NR-1]` has its top bit set
2. Initial estimate: `r[NR] = floor(2^W / d_hi)` (single-word)
3. Iterates `ceil(log₂((NR+1)·W))` times
4. Each iteration: `t = d·r`, `error = 2·B^{2n} − t`, `r_adj = r·error` (top 2n limbs become new r)

#### `div_newton`

```cpp
template<int NL, int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<NL, W, MonomialBasis>, Polynomial<NR, W, MonomialBasis>>
div_newton(const Polynomial<NL, W, MonomialBasis>& u,
           const Polynomial<NR, W, MonomialBasis>& v);
```

Unsigned division using Newton reciprocal. O(N log N) complexity.

**Algorithm:**
1. Strips leading zero limbs from `v`
2. If `eff_nr < NR` — delegates to `div_unsigned`
3. If `NR == 1` — uses `divmod_single` directly
4. If `u < v` — returns `(0, u)`
5. Computes `recip = reciprocal_newton(v)`
6. Normalizes `u` by the same shift, multiplies by `recip`
7. Extracts quotient (top NL limbs)
8. Computes remainder via `rem = u − q·v`
9. If `carry = 0` (q too large by 1): decrements q, adds v to rem
10. If `rem ≥ v`: re-divides to correct, adds correction to q

---

## 5. `witt.h` — Witt Vectors

### `WittVector<N, W>`

```cpp
template<int N, std::unsigned_integral W>
struct WittVector {
    std::array<W, N> a;
    // …
};
```

Represents a length-`N` Witt vector over ℤ/2^W.

**Members:**
- `ghost(j)` — the `j`-th ghost component: `a₀^{2ʲ} + 2·a₁^{2ʲ⁻¹} + … + 2ʲ·aⱼ`
- `ghost_vector()` — all ghost components as an array
- `frobenius()` — `F(a) = (a₀², a₁², …, a_{N−1}²)`
- `verschiebung()` — `V(a) = (0, a₀, a₁, …, a_{N−2})`
- `FV()` — Frobenius then Verschiebung
- `VF()` — Verschiebung then Frobenius
- `operator+` — Witt addition (aliases `witt_add`)
- `operator*` — Witt multiplication (aliases `witt_mul`)
- `zero()` — zero vector
- `one()` — multiplicative identity `(1, 0, …, 0)`

### Ghost Map

```cpp
WittVector<N,W>::ghost(j);
```

Computes `Σ_{i=0}^{j} 2ⁱ·aᵢ^{2ʲ⁻ⁱ}` using precomputed powers (O(N²) total).

### `check_ghost_ring`

```cpp
bool check_ghost_ring(const WittVector<N,W>& a, const WittVector<N,W>& b);
```

Verifies the ghost-ring homomorphism: `ghost(witt_add(a,b)) = ghost(a) + ghost(b)`.

### Witt Addition / Multiplication

```cpp
WittVector<N,W> witt_add(const WittVector<N,W>& a, const WittVector<N,W>& b);
WittVector<N,W> witt_mul(const WittVector<N,W>& a, const WittVector<N,W>& b);
```

Uses ghost-map + Newton recovery (O(N²)). `ghost_op` computes ghost components
at `quad_width<W>` precision, applies the operation (`+` or `×`), then recovers
Witt components via `ghost_recover`.

### Adams Operation & Teichmüller Lift

```cpp
WittVector<N,W> adams_operation(const WittVector<N,W>& a, int n);
WittVector<N,W> teichmueller_lift(W x);
```

- `adams_operation`: `ghostⱼ(ψⁿ(a)) = ghostⱼ(a)ⁿ`. Ring endomorphism; `ψⁿ ∘ ψᵐ = ψⁿᵐ`.
- `teichmueller_lift`: `τ(x) = (x, 0, …, 0)`. Embeds residue field; `τ(xy) = τ(x)·τ(y)`.

### `check_witt_recovery_precision`

```cpp
bool check_witt_recovery_precision(const WittVector<N,W>& w);
```

Checks that ghost → Witt recovery is lossless for the given vector.

### Witt Logarithm / Exponential / Inverse

```cpp
WittVector<N,W> witt_log(const WittVector<N,W>& a);      // requires a[0] odd
WittVector<N,W> witt_exp(const WittVector<N,W>& a);      // requires a[0] ≡ 0 (mod 4)
WittVector<N,W> witt_inverse(const WittVector<N,W>& a);  // requires a[0] odd
```

- `witt_log`: applies `log(1+y) = Σ (−1)^{n+1} yⁿ/n` to each ghost value, then recovers
- `witt_exp`: applies `exp(x) = Σ xⁿ/n!` to each ghost value, then recovers
- `witt_inverse`: computes `modinv_odd(ghostⱼ(a))` for each ghost, then recovers

## Verification Helpers

```cpp
bool check_witt_inverse(const WittVector<N,W>& a);           // a·a⁻¹ = 1
bool check_witt_exp_log_roundtrip(const WittVector<N,W>& a); // exp(log(a)) = a
```

---

## 6. `compose.h` — Composition & Reversion

### `compose`

```cpp
template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<(N-1)*(M-1)+1, W, MonomialBasis>
compose(const Polynomial<N, W, MonomialBasis>& P,
        const Polynomial<M, W, MonomialBasis>& Q);
```

Power series composition: `P(Q(t))` truncated to degree `(N-1)·(M-1)`.
Complexity: O(N²·M²). The result degree must be ≤ 4095 (stack guard).

### `reversion`

```cpp
template<int N, std::unsigned_integral W> requires (N >= 2)
constexpr Polynomial<N, W, MonomialBasis>
reversion(const Polynomial<N, W, MonomialBasis>& P);
```

Lagrange reversion: finds `R(t)` such that `P(R(t)) = t`.
**Precondition:** `P[1]` odd. O(N³).

---

## 7. `gcd.h` — Polynomial GCD & Friends

All operations are in the **coefficient-wise** (standard) ring, not the carry-chain ring.

### `pseudo_remainder_cw`

```cpp
Polynomial<N, W, MonomialBasis>
pseudo_remainder_cw(const Polynomial<N, W, MonomialBasis>& A,
                    const Polynomial<M, W, MonomialBasis>& B);
```

Pseudo-remainder for polynomials over ℤ₂ (works with any leading coefficient).
Returns `A·lc(B)^(deg(A)−deg(B)+1) mod B` at coefficient precision.

### `poly_divmod_cw`

```cpp
std::pair<Polynomial<N, W, MonomialBasis>, Polynomial<M, W, MonomialBasis>>
poly_divmod_cw(const Polynomial<N, W, MonomialBasis>& A,
               const Polynomial<M, W, MonomialBasis>& B);
```

Polynomial division with remainder. **Precondition:** leading coefficient of `B` must be odd.

### `poly_gcd_cw`

```cpp
Polynomial<K, W, MonomialBasis>
poly_gcd_cw(const Polynomial<N, W, MonomialBasis>& A,
            const Polynomial<M, W, MonomialBasis>& B);
```

Polynomial GCD (greatest common divisor) using Euclidean algorithm with pseudo-remainder.
Normalizes result to have odd leading coefficient.

### `det_laplace` (`namespace detail`)

```cpp
W det_laplace(const std::array<std::array<W, 6>, 6>& M, int n);
```

Laplace expansion determinant for dimensions ≤ 6.

### `polynomial_resultant_cw`

```cpp
W polynomial_resultant_cw(const Polynomial<N, W, MonomialBasis>& A,
                          const Polynomial<M, W, MonomialBasis>& B);
```

Sylvester matrix resultant. Returns 0 for shared roots. Maximum dimension 6.

### `poly_discriminant_cw`

```cpp
W poly_discriminant_cw(const Polynomial<N, W, MonomialBasis>& P);
```

The discriminant `Δ(P) = res(P, D(P)) / lc(P)`. Returns 0 for non-odd leading coefficient.

### `poly_is_square_free_cw`

```cpp
bool poly_is_square_free_cw(const Polynomial<N, W, MonomialBasis>& P);
```

Checks if `gcd(P, D(P))` is constant.

---

## 8. `combinatorial.h` — Binomial & Stirling Numbers

### `binom`

```cpp
template<std::unsigned_integral W, int MaxN = 67>
constexpr W binom(int n, int k);
```

Binomial coefficient `C(n, k)` using multiplicative GCD-based reduction.
MaxN=67 (fits in uint64) or 100 (with `__int128`). Out-of-range returns 0.

### `stirling_2`

```cpp
template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_2(int n, int k);
```

Stirling numbers of the second kind `S₂(n,k)` using DP recurrence.

### `stirling_1`

```cpp
template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_1(int n, int k);
```

Signed Stirling numbers of the first kind `s(n,k)`. Returned as unsigned (2-adic
wrapping: `s(2,1) = −1` wraps to `UINT64_MAX`).

### `stirling_1_unsigned`

```cpp
template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_1_unsigned(int n, int k);
```

Unsigned Stirling numbers of the first kind `|s(n,k)|`.

---

## 9. `dynamic_polynomial.h` — Runtime-Degree Polynomial

### `DynamicPolynomial<W, Basis>`

Runtime-polynomial backed by `std::vector<W>`. Supports all bases, conversion
to/from compile-time `Polynomial<N,W,Basis>`, carry-chain multiplication, eval
(all bases), basis conversion, formal derivative, and forward difference.

**Key functions:**
- `to_dynamic(p)` — convert compile-time to runtime
- `to_static<N>(p)` — convert runtime to compile-time (truncates/pads)
- `change_basis<ToBasis>(p)` — basis conversion over runtime
- `formal_derivative(p)`, `forward_difference(p)` — derivative/difference over runtime

### `DynamicStirlingCache<W>`

Runtime analogue of `StirlingCache`, built with `std::vector`.

---

## 10. `matrix.h` — Matrix / Linear Algebra over ℤ₂

### `Matrix<M, N, W>`

```cpp
template<int M, int N, std::unsigned_integral W>
struct Matrix { … };
```

**Members:**
- `identity()` / `zero()` — named constructors
- `transpose()` — matrix transpose
- `trace()` — sum of diagonal (requires square)
- `determinant()` — Bareiss-like 2-adic determinant (avoids division by even pivots)
- `inverse()` — augmented-matrix Gauss-Jordan
- `rref()` — reduced row echelon form; returns `(RREF, rank, det_sign)`
- `rank()` — rank over ℤ₂ (parity-only pivoting)
- `solve(b)` — linear system solver
- `is_singular()` — determinant is even
- `operator+`, `operator-`, `operator*` — standard matrix arithmetic
- `operator==`, `operator!=` — element-wise comparison

---

## 11. `pade.h` — Padé Approximants

### `pade_approximant`

```cpp
template<int M, int N, std::unsigned_integral W>
constexpr std::pair<Polynomial<M+1, W, MonomialBasis>,
                    Polynomial<N+1, W, MonomialBasis>>
pade_approximant(const Polynomial<M+N+1, W, MonomialBasis>& series);
```

Computes `[M/N]` Padé approximant of a power series over ℤ₂[[t]].
Uses Gaussian elimination on the Toeplitz system. Returns `(P, Q)` where
`P(t) / Q(t) = series(t) + O(t^{M+N+1})`. Returns empty pairs if singular.

---

## 12. `continued_fractions.h` — Continued Fractions

### `cf_expand`

```cpp
std::vector<W> cf_expand(const DynamicPolynomial<W, MonomialBasis>& series, int max_terms);
```

Extract continued fraction coefficients from a power series (extracts `s[0]`,
shifts, repeats).

### `cf_convergent`

```cpp
std::pair<DynamicPolynomial<W, MonomialBasis>, DynamicPolynomial<W, MonomialBasis>>
cf_convergent(const std::vector<W>& c, int n);
```

Compute the `n`-th convergent numerator/denominator from CF coefficients.

### `cf_eval`

```cpp
DynamicPolynomial<W, MonomialBasis> cf_eval(const std::vector<W>& c, int max_terms, W x);
```

Evaluate the `max_terms`-th convergent at `x` (returns empty if denominator even at `x`).

---

## 13. `verify.h` — Compile-Time Verification

### `check_taylor_roundtrip_precision`

```cpp
template<int N, std::unsigned_integral W>
constexpr bool check_taylor_roundtrip_precision(const Polynomial<N,W,MonomialBasis>& p);
```

Returns true if `Mono → Taylor → Mono` is lossless. Fails when factorial
multiplication wraps: `FF_k ≥ 2^W / k!`.

### `check_witt_recovery_precision`

See §5 — checks if ghost→Witt recovery is lossless.

### Proof Engine

```cpp
struct prove_all<Pred, Ns<NN…>, Ws<WW…>>;
```

Template metaprogram for auto-generating proof instances across combinations
of polynomial degrees `NN…` and word widths `WW…`. Each proof value is a
`static_assert`. Controlled by `DYADIC_HEAVY_PROOFS` to enable/disable
exhaustive checks.

### Compile-Time Proofs (24 named proofs, ~55 assertions)

Covers:
- `PROOF_GHOST_HOM_*` — ghost-ring addition homomorphism (various N, W)
- `PROOF_TAYLOR_ROUNDTRIP_*` — Taylor basis roundtrip for small coefficients
- `PROOF_BINOM_*` — binomial identity `C(n,k) = C(n-1,k) + C(n-1,k-1)`
- `PROOF_STIRLING_2_*`, `PROOF_STIRLING_1_*` — Stirling number recurrence identities
- `PROOF_WITT_INVERSE_*` — Witt vector inverse
- `PROOF_WITT_LOG_EXP_*` — Witt exp(log(·)) roundtrip
- `PROOF_ADAMS_TEICH_*` — Adams operation on Teichmüller lifts
- `PROOF_CARRY_*` — carry-chain idempotency
- `PROOF_SYNTHETIC_DIVIDE_*` — synthetic divide properties
- `PROOF_ARTIN_SCHREIER_*` — Artin-Schreier symmetry verification
