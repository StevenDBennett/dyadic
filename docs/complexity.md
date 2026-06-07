# dyadic — Algorithmic Complexity Reference

Comprehensive breakdown of time and space complexity for every function in the
library. `N`, `M`, `P`, `K` denote compile-time template parameters (polynomial degree,
matrix dimensions). `da`, `db` denote actual degrees at runtime. `B` = 8 × sizeof(W)
(bit width of the word type). All operations are `constexpr`.

---

## 1. `core.h` — Core Types & Primitives

### 2-Adic Primitives

| Function | Time | Space | Notes |
|---|---|---|---|
| `v2(W x)` | O(1) | O(1) | `std::countr_zero` → single CPU instruction (`tzcnt`/`bsf`) |
| `modinv_odd(W a)` | O(log B) ≈ O(1) | O(1) | Newton/Hensel iteration: `⌈log₂(B)⌉` iterations, 2 multiplies per iteration. For uint64: 6 iterations. |
| `exact_divide(W x, W d)` | O(log B) ≈ O(1) | O(1) | 1× `modinv_odd` + 1× multiply |
| `div_2k_adic(U val, int k)` | O(1) | O(1) | 3 branches + shift + optional sign-extension mask |
| `artin_schreier(W x)` | O(1) | O(1) | `x*x − x` |

### `detail::uint128_t` (Software 128-bit)

| Function | Time | Space | Notes |
|---|---|---|---|
| `+`, `-`, `*`, `<<`, `>>`, `&`, `\|`, `^`, `~`, comparisons | O(1) | O(1) | 1–2 word ops |
| `*` (multiply) | O(1) | O(1) | 12 half-width multiplies via decomposed 64×64→128 schoolbook |
| `/`, `%` (division) | O(B) = O(128) | O(1) | Schoolbook shift-and-subtract; 128 iterations |

### Compile-Time Caches

| Type | Constructor | Time | Space | Notes |
|---|---|---|---|---|
| `StirlingCache<N,W>` | ctor | O(N²) compile-time | O(N²) | Two N×N tables (S₂, s₁); DP recurrence; computed once per (N,W) |
| `PascalCache<N,W>` | ctor | O(N²) compile-time | O(N²) | One N×N table (binomial); computed once per (N,W) |

Both are `inline constexpr` globals — zero runtime cost after instantiation.

### `Polynomial<N, W, Basis>`

| Function | Time | Space | Notes |
|---|---|---|---|
| `Polynomial()` | O(1) | O(N) (embedded) | Value-initialized array |
| `Polynomial(std::array)` | O(N) | O(N) | Copy from array |
| `actual_degree()` | O(N) worst, O(1) best | O(1) | Linear scan from N−1 down to 0 |
| `eval(W x)` | O(N) | O(1) | Horner's method; monomial-only (`static_assert` guard) |
| `operator+`, `operator−` | O(N) | O(N) | Element-wise |
| `operator*` | O(N·M) | O(N+M) | **Carry-chain** multiply → `poly_mul`; NOT the coefficient-wise ring |

### Carry Chain

| Function | Time | Space | Notes |
|---|---|---|---|
| `synthetic_divide` (Ring) | O(1) | O(1) | Tag dispatch, no computation |
| `synthetic_divide` (Polynomial) | O(N) | O(N) in-place | Half-word carry propagation: (I−N)⁻¹; one pass over N coefficients |
| `carry_normalize` (all) | O(1) | O(1) | Tag conversion |
| `carry_chain(W* r, Accum* v, int N)` | O(N) | O(N) | Single pass: `r[i] = low(v[i]+carry)`, `carry = high(v[i]+carry)` |
| `carry_chain_word(W hi, W lo)` | O(1) | O(1) | Single-word carry |
| `adc` / `add_overflow` / `mul_overflow` | O(1) | O(1) | Compiler builtins or manual compare |

### Axiom Operations

| Function | Time | Space | Notes |
|---|---|---|---|
| `formal_derivative` (Monomial) | O(N) | O(N−1) | D(xⁿ) = n·xⁿ⁻¹; single loop, coefficient scaling |
| `forward_difference` (Monomial) | O(N²) | O(N−1) | Double nested loop: Σ_{i,j} p[j]·C(j,i) using Pascal cache |
| `exact_formal_derivative` | O(N) | O(N−1) | ≡ `formal_derivative` in monomial basis |
| `exact_forward_difference` | O(N²) | O(N−1) | ≡ `forward_difference` in monomial basis |

### Polynomial Multiplication

| Function | Time | Space | Notes |
|---|---|---|---|
| `poly_mul_unsaturated` (detail) | O(na·nb) | O(na+nb) | `quad_width` accumulators (hw `unsigned __int128` for uint64_t when available), `#pragma GCC ivdep` for auto-vectorization |
| `poly_mul(r, a, na, b, nb)` | O(na·nb) | O(1) chunked | `CHUNK_COUNT = 256/sizeof(accum_t)`; for uint64_t with `__SIZEOF_INT128__` uses hw `unsigned __int128` (~2.7× speedup deg 63 vs software); single-pass + one carry_chain per chunk |
| `poly_mul_cw` | O(N·M) | O(N+M) | Coefficient-wise convolution (standard ring); zero-skip for sparsity |
| `verify_divmod_cw` | O(N·M) | O(N+M) | One `poly_mul_cw` + verification loop |

---

## 2. `basis.h` — Basis Conversion

### Low-Level Conversions

All conversions are triangular matrix-vector products using Stirling number caches.

| Function | Time | Space | Pattern |
|---|---|---|---|
| `monomial_to_falling` | O(N²) | O(N) | Upper-triangular mat-vec multiply by S₂(n,k); `dword_t` accum |
| `falling_to_monomial` | O(N²) | O(N) | Upper-triangular mat-vec multiply by s₁(n,k) (signed, stored unsigned) |
| `monomial_to_taylor` | O(N²) | O(N) | S₂(n,k) × k!; factorial accumulated incrementally |
| `taylor_to_monomial` | O(N²) | O(N) | s₁(j,k) × p[j]/j!; parity branch: `modinv_odd` vs `div_2k_adic`+`modinv_odd`; `quad_width` accum |

### High-Level Dispatch

| Function | Time | Space | Notes |
|---|---|---|---|
| `change_basis<ToBasis>` | O(N²) | O(N) | Template dispatch; FF↔Taylor routes through Monomial (2 conversions) |

### Basis-Specific Operations

| Function | Time | Space | Notes |
|---|---|---|---|
| `eval` (FallingFactorial) | O(N) | O(1) | Incremental falling product: Π(x−j) |
| `eval` (Taylor) | O(N) | O(1) | Incremental binomial coefficient: (x choose i) |
| `formal_derivative` (FF) | O(N²) | O(N) | Convert → Monomial D → convert back (no direct FF derivative formula) |
| `formal_derivative` (Taylor) | O(N²) | O(N) | Convert → Monomial D → convert back |
| `forward_difference` (FF) | O(N) | O(N) | **Direct**: Δ(FFₙ) = n·FFₙ₋₁ — FF diagonalises Δ |
| `forward_difference` (Taylor) | O(N) | O(N) | **Direct**: Δ(Tₙ) = Tₙ₋₁ — Taylor nearly diagonalises Δ |

---

## 3. `calculus.h` — Taylor Shift, Indefinite Sum, Exp/Log

### Taylor Shift

| Overload | Time | Space | Notes |
|---|---|---|---|
| `taylor_shift` (Monomial) | O(N²) | O(N) | Double nested loop using Pascal cache and δⁿ⁻ᵏ powers |
| `taylor_shift` (FallingFactorial) | O(N²) | O(N) | Uses falling factorial of delta: (δ)₍ₖ₋ᵢ₎ |
| `taylor_shift` (Taylor) | O(N²) | O(N) | Converts to monomial, shifts, converts back (3× constant factor) |

### Indefinite Sum (Σ = Δ⁻¹)

| Overload | Time | Space | Notes |
|---|---|---|---|
| `indefinite_sum` (FallingFactorial) | O(N) | O(N) | Σ(FFₙ₋₁) = FFₙ/n; single loop with 2-adic division per term |
| `indefinite_sum` (Monomial) | O(N²) | O(N) | Convert → FF sum → convert back |
| `indefinite_sum` (Taylor) | O(N²) | O(N) | Convert → Monomial → FF sum → Monomial → Taylor (4× constant) |

### Power Series Exp/Log

| Function | Time | Space | Notes |
|---|---|---|---|
| `poly_exp` | O(N²) | O(N) | Recurrence: Eₙ = (1/n)·Σ k·Pₖ·Eₙ₋ₖ; `dword_t` intermediates; precondition P[1] even |
| `poly_log` | O(N²) | O(N) | Recurrence: Lₙ = Pₙ − (1/n)·Σ k·Lₖ·Pₙ₋ₖ; `dword_t` intermediates; precondition P[1] even |

---

## 4. `witt.h` — Witt Vectors

### Ghost Map Infrastructure

| Function | Time | Space | Notes |
|---|---|---|---|
| `GhostPowers<N,W>::init` | O(N²) | O(N²) | Precomputes a[i]^(2^k) for all i, k; stores N×N W-precision table |
| `GhostPowersFull<N,W>::init` | O(N²) | O(N²) | Same but at `quad_width<W>` precision |
| `GhostPowersFull::ghost(j)` | O(j) ⊆ O(N) | O(1) | Ghost component: Σ 2ⁱ·aᵢ^(2ʲ⁻ⁱ) |
| `ghost_j_widen` | O(j) ⊆ O(N) | O(1) | Widened ghost: `quad_width` accumulation |
| `ghost_j_full` | O(j²) ⊆ O(N²) | O(1) | Full-precision ghost (recomputes powers, no precomputed cache) |

### WittVector Operations

| Function | Time | Space | Notes |
|---|---|---|---|
| `WittVector::ghost_vector()` | O(N²) | O(N) | N ghost components, each O(N) → total O(N²) |
| `WittVector::ghost(j)` | O(N²) | O(N) | Calls `ghost_vector()[j]` |
| `WittVector::frobenius()` | O(N) | O(N) | Componentwise square |
| `WittVector::verschiebung()` | O(N) | O(N) | Shift right (insert 0 at index 0) |
| `WittVector::FV()` / `VF()` | O(N) | O(N) | Two O(N) operations |

### Witt Ring Operations

| Function | Time | Space | Notes |
|---|---|---|---|
| `ghost_recover` | O(N²) | O(N²) | Newton recovery: rⱼ = (Gⱼ − Sⱼ)/2ʲ; nested loops; stores N×N `quad_width` powers |
| `ghost_op` | O(N²) | O(N²) | Two `GhostPowersFull::init` + N ghost comps + `ghost_recover` |
| `witt_add` | O(N²) | O(N²) | `ghost_op` with `std::plus` |
| `witt_mul` | O(N²) | O(N²) | `ghost_op` with `std::multiplies` |
| `adams_operation` | O(N² + N·n) | O(N²) | Ghost component ^ n powering; early return for n ≤ 1 |
| `teichmueller_lift` | O(N) | O(N) | Zero-fills after index 0 |

### Witt Log/Exp/Inverse

| Function | Time | Space | Notes |
|---|---|---|---|
| `p_adic_log_1plus_impl` | O(B) | O(1) | T = 2·B ≈ 128–256 terms; Σ (−1)ⁿ⁺¹yⁿ/n; early exit on vanishing. For W=uint64_t with `__SIZEOF_INT128__` uses hardware `unsigned __int128` multiply-accumulate and `modinv_odd_128` (full 128-bit inverse) for ~4× speedup over software `detail::uint128_t` |
| `p_adic_exp_impl` | O(B) | O(1) | T ≤ 2·B terms; Σ xⁿ/n!; valuation-aware budget; early exit. Same `unsigned __int128` optimization as `p_adic_log` |
| `witt_log` | O(N² + N·B) | O(N²) | `GhostPowersFull::init` + N × `p_adic_log` + `ghost_recover`; requires a[0] odd |
| `witt_exp` | O(N² + N·B) | O(N²) | `GhostPowersFull::init` + N × `p_adic_exp` + `ghost_recover`; requires a[0] ≡ 0 mod 4 |
| `witt_inverse` | O(N²) | O(N²) | `GhostPowersFull::init` + N × `modinv_odd` on ghost + `ghost_recover` |

### Verification Helpers

| Function | Time | Space | Notes |
|---|---|---|---|
| `check_ghost_ring` | O(N²) | O(N²) | 3× `ghost_vector` + O(N) loop |
| `check_witt_recovery_precision` | O(N²) | O(N²) | 2× `ghost_vector` + `ghost_recover` |
| `check_witt_inverse` | O(N²) | O(N²) | `witt_inverse` + `operator*` + O(N) verify |
| `check_witt_exp_log_roundtrip` | O(N² + N·B) | O(N²) | `witt_log` + `witt_exp` |

---

## 5. `compose.h` — Composition & Reversion

| Function | Time | Space | Notes |
|---|---|---|---|
| `compose` | O(N²·M²) | O(N·M) | Naive power accumulation: maintain Qᵏ, accumulate Pₖ·Qᵏ; result degree (N−1)(M−1); guarded ≤ 4095 |
| `reversion` | O(N³) | O(N²) | Classical Lagrange inversion (NOT Newton). Phase 1: O(N³) 2D power table. Phase 2: O(N³) incremental update. Requires P[1] odd. |

---

## 6. `gcd.h` — Polynomial GCD & Friends

| Function | Time | Space | Notes |
|---|---|---|---|
| `pseudo_remainder_cw` | O(da·(da+db)) ⊆ O(N·(N+M)) | O(N) | Subresultant PRS; multiplies by (lcB)ᵏ; no division by lc; works with any lc |
| `poly_divmod_cw` | O(da·(da+db)) ⊆ O(N·(N+M)) | O(max(N,M)) | Exact long division; requires odd lc (`modinv_odd` on leading coefficient) |
| `poly_gcd_cw` | O(K³) worst case, K=max(N,M) | O(K) | Euclidean PRS; up to K `pseudo_remainder_cw` calls of O(d²) each → Σd² ≈ K³/3; normalizes result to monic |
| `det_laplace` (detail) | O(n!) ≤ O(720) | O(n²) stack | Recursive Laplace expansion; bounded by all callers to n ≤ 6 |
| `polynomial_resultant_cw` | O(1) bounded | O(1) | Sylvester matrix (dim ≤ 6) + `det_laplace`; returns 0 for dim > 6 |
| `poly_discriminant_cw` | O(N) + O(1) | O(N) | `formal_derivative` + `polynomial_resultant_cw`; returns 0 for even lc or deg > 6 |
| `poly_is_square_free_cw` | O(N³) | O(N) | `poly_gcd_cw(P, P′)`, checks if degree = 0 |

---

## 7. `combinatorial.h` — Binomial & Stirling Numbers

| Function | Time | Space | Notes |
|---|---|---|---|
| `detail::binom_gcd` | O(log min(a,b)) | O(1) | Euclidean GCD, iterative |
| `binom<W,MaxN>(n,k)` | O(k·log n) ⊆ O(1) bounded | O(1) | Product-of-fractions with GCD per iteration to avoid overflow; MaxN=67 (uint64) or 100 (with __int128) |
| `stirling_2<W,MaxN>(n,k)` | O(n²) | O(MaxN) | 1D DP table; recurrence: S(n,k) = S(n−1,k−1) + k·S(n−1,k) |
| `stirling_1<W,MaxN>(n,k)` | O(n²) | O(MaxN) | Same pattern; signed recurrence: s(n,k) = s(n−1,k−1) − (n−1)·s(n−1,k) |
| `stirling_1_unsigned<W,MaxN>(n,k)` | O(n²) | O(MaxN) | Same as `stirling_2` but unsigned first kind |

---

## 8. `arith.h` — Big-Integer Arithmetic on Polynomial

### Internal Helpers

| Function | Time | Space | Notes |
|---|---|---|---|
| `zero_extend<M>` | O(N) | O(M) | Copy N limbs, zero-fill to M |
| `unsigned_lt` | O(N) | O(1) | Lexicographic MSB comparison |
| `shift_left` | O(N) | O(N) | Word-level shift + intra-word bit-level shift with carry |
| `clz_normalize` | O(1) | O(1) | `std::countl_zero` |

### Unsigned Multiplication

| Function | Time | Space | Notes |
|---|---|---|---|
| `mul_unsigned` | O(N·M) | O(1) chunked | `quad_width` or `unsigned __int128` unsaturated convolution + per-chunk carry chain; produces N+M limbs. For uint64_t with `__SIZEOF_INT128__` uses hardware `unsigned __int128` (same chunked pattern as `poly_mul`); else falls back to software `quad_width<W>` full-buffer accumulation |

### Unsigned Division

| Function | Time | Space | Notes |
|---|---|---|---|
| `divmod_single` | O(NL) | O(NL) | Schoolbook long division by single-word divisor; `dword_t` accum |
| `div_unsigned_knuth` | O(NL·NR) | O(NL+NR) | Knuth Algorithm D (TAOCP §4.3.1); normalization, per-digit quotient estimation with correction, denormalization. For W=uint64_t with `__SIZEOF_INT128__` uses hardware `unsigned __int128` for `dw_t` (trial division, multiply-subtract inner loop) giving ~3-4× speedup over software `detail::uint128_t` |
| `div_unsigned` | O(NL·NR) worst / O(NL) (NR=1) | O(NL+NR) | Dispatch wrapper; strips leading zeros; delegates to `divmod_single` or `div_unsigned_knuth` |

### Newton Division

| Function | Time | Space | Notes |
|---|---|---|---|
| `reciprocal_newton` | O(NR²·log NR) | O(NR) | xₖ₊₁ = xₖ·(2·B²ⁿ − d·xₖ); ⌈log₂((NR+1)·B)⌉ iterations, each with 2 full convolutions |
| `div_newton` | O(NR²·log NR + NL·NR) | O(NL+NR) | Normalize both operands; q = u_shifted × recip >> 2·NR; correction step for exact remainder |

### Carry-Save Arithmetic

| Function | Time | Space | Notes |
|---|---|---|---|
| `carry_save_add` | O(N) | O(N) | 2→2 CSA; carry[i] ∈ {0,1} |
| `csa3` | O(N) | O(N) | 3→2 CSA; carry[i] ∈ {0,1,2} |
| `combine_carry_save` | O(N) | O(N) | Final carry propagation; overflow beyond N discarded |
| `batched_add_carry_save` | O(count·N) | O(N) | Iterative CSA of `count` operands; single final carry propagation |

---

## 9. `dynamic_polynomial.h` — Runtime-Degree Polynomial

| Function | Time | Space | Notes |
|---|---|---|---|
| `DynamicStirlingCache::DynamicStirlingCache(N)` | O(N²) | O(N²) | Builds 2D `std::vector` tables for S₂ and s₁ |
| `DynamicPolynomial` constructors | O(size) | O(size) | Copy/initialise from various sources |
| `operator[]` | O(1) | O(1) | Direct vector subscript |
| `degree()` | O(size) | O(1) | Linear scan from back to strip trailing zeros |
| `operator+=`, `operator-=` | O(max(n,m)) | O(max(n,m)) amortized | Element-wise ring op; reallocates if needed |
| `eval` (Monomial) | O(size) | O(1) | Horner's method |
| `eval` (FallingFactorial) | O(size) | O(1) | Iterative falling factorials |
| `eval` (Taylor) | O(size) | O(1) | Iterative binomial coefficients + 2-adic division |
| `operator*` | O(size_a·size_b) | O(size_a+size_b) | Delegates to `poly_mul`; carry-chain sematics |
| `change_basis` | O(N²) | O(N²) | Builds `DynamicStirlingCache` + triangular mat-vec |
| `formal_derivative` (Monomial) | O(size) | O(size) | Scale-and-shift |
| `formal_derivative` (FF/Taylor) | O(N²) | O(N²) | Convert → differentiate → convert back |
| `forward_difference` (Monomial) | O(N²) | O(N) | Binomial transform from cache |
| `forward_difference` (FF) | O(size) | O(size) | Direct: Δ(FFₙ) = n·FFₙ₋₁ |
| `to_dynamic` / `to_static<N>` | O(N) | O(N) / O(1) | Copy between compile-time and runtime representations |

---

## 10. `matrix.h` — Matrix / Linear Algebra over ℤ₂

All complexity in terms of dimensions M, N, P. All matrices are compile-time sized.

| Function | Time | Space | Notes |
|---|---|---|---|
| `identity()` (M==N) | O(M) | O(M·N) | Diagonal 1s |
| `zero()` | O(1) | O(M·N) | Value-initialized |
| `operator+`, `operator−` | O(M·N) | O(M·N) | Element-wise |
| `operator==`, `operator!=` | O(M·N) | O(1) | Element-wise; early exit on mismatch |
| `operator*` (M×N × N×P) | O(M·N·P) | O(M·P) | i-k-j loop with zero-skip on A[i][k] |
| `transpose()` | O(M·N) | O(N·M) | Element-wise copy |
| `trace()` (M==N) | O(M) | O(1) | Sum of diagonal |
| `is_singular()` | O(M³) | O(M·N) | Calls `determinant()` |
| `rref()` | O(M·N·min(M,N)) | O(M·N) | Gauss–Jordan; pivot search with `modinv_odd`; returns (RREF, rank, det_sign) |
| `rank()` | O(M·N·min(M,N)) | O(M·N) | ℤ₂-only: parity checks, XOR-style row ops; per-column early-out |
| `determinant()` (M==N) | O(M³) | O(M·N) | **Bareiss** fraction-free elimination: 2×2 determinant recurrence, divides by previous pivot via `dword_t`; parity sign tracking |
| `inverse()` (M==N) | O(M³) | O(M²) | Gauss–Jordan on [A\|I] with 2×M row width; `modinv_odd` row scaling |
| `solve(b)` (M==N) | O(M³) | O(M·(N+1)) | Gauss–Jordan on [A\|b] augmented system |

---

## 11. `pade.h` — Padé Approximants

| Function | Time | Space | Notes |
|---|---|---|---|
| `pade_approximant<M,N>` | O(N³ + M·N) | O(N² + M+N) | Builds N×N Toeplitz-like system from series; Gaussian elimination with odd-pivot `modinv_odd`; back-substitute for denominator q; convolution for numerator p. N=0 special case is O(M). |

---

## 12. `continued_fractions.h` — Continued Fractions

| Function | Time | Space | Notes |
|---|---|---|---|
| `cf_expand(series, max_terms)` | O(max_terms·S) | O(S) | Extract s[0], shift, repeat; S = `series.size()` |
| `cf_convergent(c, n)` | O(n²) | O(n) | Classical recurrence: Pₖ = cₖ·Pₖ₋₁ + Pₖ₋₂; degree grows linearly; total work Σk = O(n²) |
| `cf_eval(c, max_terms, x)` | O(n²) | O(n) | `cf_convergent` + evaluate P(x), Q(x); returns P(x)/Q(x) if Q(x) odd, else empty |

---

## 13. `prng.h` — XorShift64 PRNG

| Function | Time | Space | Notes |
|---|---|---|---|
| `XorShift64(s)` | O(1) | O(1) | Seed single uint64_t |
| `next()` | O(1) | O(1) | Triple xorshift: `^= <<13`, `^= >>7`, `^= <<17`; full cycle 2⁶⁴−1 |
| `is_valid_seed(s)` | O(1) | O(1) | `static constexpr`, checks s ≠ 0 |

---

## Summary by Complexity Class

| Class | Functions |
|---|---|
| **O(1)** | `v2`, `modinv_odd`, `exact_divide`, `div_2k_adic`, `artin_schreier`, `carry_chain_word`, `adc`, `uint128_t` arithmetic, `XorShift64::next`, ring tag dispatches |
| **O(N)** | `eval` (all bases), `formal_derivative` (Monomial), `forward_difference` (FF/Taylor), `indefinite_sum` (FF), `synthetic_divide`, `carry_chain`, `actual_degree`, `WittVector::frobenius/verschiebung`, `teichmueller_lift`, `divmod_single`, `carry_save_add`, `csa3`, `combine_carry_save`, `to_dynamic`/`to_static` |
| **O(N·M)** | `poly_mul` / `operator*`, `poly_mul_cw`, `mul_unsigned`, `div_unsigned_knuth`, `div_unsigned`, `Matrix::operator*`, `verify_divmod_cw` |
| **O(N²)** | `forward_difference` (Monomial), `taylor_shift` (all bases), `indefinite_sum` (Monomial/Taylor), `poly_exp`, `poly_log`, `change_basis`, `formal_derivative` (FF/Taylor), `GhostPowers::init`, `ghost_vector`, `witt_add`/`witt_mul`/`witt_inverse`, `witt_log`/`witt_exp` (plus N·B), `adams_operation`, `binom` runtime, `stirling_2`/`stirling_1`, `cf_convergent`, `StirlingCache`/`PascalCache` ctor (compile-time) |
| **O(N²·M²)** | `compose` |
| **O(N³)** | `reversion`, `Matrix::determinant`/`inverse`/`rref`/`solve`/`rank`/`is_singular`, `pade_approximant` (dominates at O(N³)), `poly_gcd_cw` (worst case) |

---

## Key Design Decisions Affecting Complexity

1. **No FFT/NTT over ℤ₂**: Carry-chain multiplication is inherently O(N·M). The standard NTT requires a principal root of unity, which does not exist in ℤ/2^Wℤ for W > 1. Chunked convolution with auto-vectorization is the practical limit.

2. **Newton vs Lagrange for reversion**: The library uses classical Lagrange inversion (O(N³)) rather than Newton/Rothe–Henry (O(N log² N)). This is a deliberate choice for `constexpr` simplicity — Newton requires fast multiplication, which would need a sub-quadratic algorithm.

3. **Ghost-op + Newton recovery for Witt ops**: All Witt operations go through the ghost map (diagonalise the operation), then Newton-recover the Witt components. This is O(N²) but avoids per-component Hensel lifting.

4. **Bareiss elimination for determinant**: The fraction-free Bareiss algorithm avoids division by even pivots (which would fail in ℤ/2^Wℤ). The division by the previous pivot is done via `dword_t` → 2-adic division → `modinv_odd`.

5. **Stirling/Pascal caches at compile time**: Both are computed at compile time via `inline constexpr` globals computed once per (N,W). This avoids redundant O(N²) recomputation across multiple basis conversions.

6. **Chunked stack buffer in `poly_mul` and `mul_unsigned`**: Constant stack space O(1) regardless of polynomial size, using a `CHUNK_COUNT`-sized buffer (≤ 32 elements for uint64). Eliminates full-size `quad_width` array allocation for large products.

7. **`quad_width<W>` accumulator precision**: Used for ghost-map accumulation and unsaturated polynomial products. Provides 4× word width headroom, preventing overflow in intermediate sums. On uint64_t hot paths (`poly_mul`, `mul_unsigned`, `div_unsigned_knuth`, p-adic series), `unsigned __int128` replaces `detail::uint128_t` where available for 3-4× hardware-accelerated arithmetic.

8. **All `constexpr`**: No runtime-only code paths. The entire library, including `poly_mul`, carry chains, Witt arithmetic, and division, supports compile-time evaluation.

---

## Memory Usage Notes

- **Polynomial storage**: `Polynomial<N,W>` is `std::array<W,N>` — exactly `N × sizeof(W)` bytes embedded in the object. No heap allocation.
- **Compile-time caches**: `StirlingCache<N,W>` = `2 × N × N × sizeof(W)` bytes in `.rodata`. `PascalCache<N,W>` = `N × N × sizeof(W)` bytes in `.rodata`. Both are `inline constexpr` shared across all TUs.
- **`poly_mul` / `mul_unsigned` chunked buffer**: `CHUNK_COUNT × sizeof(accum_t)` bytes on stack. For uint64_t with `unsigned __int128`: 32 × 16 = 512 bytes. For smaller word widths: proportionally less.
- **Witt ghost_recover**: Allocates `N × N × sizeof(quad_width<W>)` on the stack for the power table. For (N=32, W=uint64, quad_width=uint128_t): 32×32×16 = 16 KiB.
- **`compose`**: Allocates `3 × result_degree × sizeof(W)` on the stack. Result degree max 4095, so up to ~96 KiB for uint64_t.
- **`reciprocal_newton`**: Allocates `4 × NR` limbs on the stack, capped at 65536 bytes by `static_assert`.
- **`DynamicPolynomial`**: Uses `std::vector<W>` — heap-allocated. All functions using it inherit vector allocation costs.