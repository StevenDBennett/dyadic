# Development Notes

## What constexpr proofs caught

The compile-time proofs are not ceremonial — they found real bugs during development:

- `core.h` was missing `artin_schreier()` entirely; `PROOF_AS_KERNEL` and `PROOF_AS_SYMMETRY_8` forced its addition.
- `stirling_1` returns signed values as unsigned wrapping; the proof `s(3,1) = 2` caught a sign error in the expected value.
- `div_2k_adic` was originally a plain right shift; `taylor_to_monomial` produced wrong results for negative 2-adic values, caught by the roundtrip proofs.
- `reversion`'s initial powers were zero-initialised but never computed before the first iteration, producing garbage coefficients.
- Integer promotion in `artin_schreier_symmetry` (`W(1) - x` promoted to `int` for `uint8_t`) produced wrong values for x > 128.

## Integer promotion and UBSan

C++ integer promotion rules silently promote `uint8_t` and `uint16_t` operands to `signed int` before arithmetic. For products like `W(51813) * W(51813)` (with W = `uint16_t`), the result is computed in `int` and overflows INT_MAX — undefined behavior. UBSan flags these at runtime with ASan+UBSan builds.

The fix is to force unsigned arithmetic by multiplying through `1u`: `static_cast<W>(1u * a * a)`. The `1u` pulls the expression into `unsigned int`, where overflow is defined behavior (wrapping modulo 2^32). The final `static_cast<W>` truncates to the target width, matching the original mathematical intent. This pattern is now applied throughout the library.

## The `unsigned __int128` friction points

A deliberate design choice that exposed several GCC-specific rough edges:

- `unsigned __int128` does not satisfy `std::unsigned_integral`, so any standard algorithm or concept constrained on that built-in concept rejects it. Every template constrained with `std::unsigned_integral` in this library silently disallows `__int128`.
- `std::gcd` (libstdc++) refuses to instantiate on `unsigned __int128`. The local `binom_gcd` function works around this.
- This is not a bug in GCC or the library — it is an inherent tension between the C++ standard's unwillingness to acknowledge extended integer types and their necessity for practical 2-adic arithmetic.
- `detail::uint128_t` avoids all these friction points: the software pair satisfies `std::unsigned_integral` (`uint64_t`-based), works with `std::gcd`, and has no `-Wpedantic` warnings.

Despite these friction points, `poly_mul` uses `unsigned __int128` as its inner
accumulator for `W = uint64_t` when `__SIZEOF_INT128__` is defined (guarded by
`#if` preprocessor, not in a `std::unsigned_integral`-constrained context).
The accumulator type is not `W` — it's an internal `accum_t` alias — so the
standard concept constraints don't apply. This gives a ~2.7× speedup at deg=63
while keeping `detail::uint128_t` as the universal fallback for all compilers.

## The pattern for C-style arrays in templates

The original code used raw C arrays sized by template parameters (e.g., `W dp[64][64]`, `W powers[64][64]`). These work in practice but are fragile: they waste stack space (64×64 for small n), the bound is decoupled from the template parameter, and `std::array` provides `.size()`, iterator support, and uniform initialization. The refactor to `std::array` shrinks stack usage (e.g., `powers` went from 64×64 to N×N) and makes the intent explicit.

## O(N³) → O(N²): a case study in template metaprogramming optimization

The original Witt addition recomputed `a[i]^(2^k)` for each ghost sum `j`, yielding O(N³) total operations. The `GhostPowers<N,W>` helper precomputes `a[i]^(2^k)` once per Witt vector in O(N²), and each ghost sum then reads precomputed values in O(N). This is a natural fit for template metaprogramming: the N is known at compile time, so the precomputation buffer is a fixed-size array, and all loops are fully unrollable by the optimizer.

## Complexity wins during refactoring

Several additional hotspots were identified and fixed:

- **`taylor_to_monomial`: O(N³) → O(N²)** — factorial was recomputed from scratch (`for i=2..j`) inside the inner j-loop of a double sum. Precomputed once into a `facts[N]` array.

- **`reversion`: O(N⁴) → O(N³)** — the power matrix `R, R², ..., R^{N-1}` was fully recomputed from scratch each time a new coefficient `R[n]` was set. With incremental updates, only the O(N²) entries affected by the new `R[n]` are propagated through the power hierarchy.

- **Stirling functions: O(N²) space → O(N)** — `stirling_2`, `stirling_1`, and `stirling_1_unsigned` originally allocated full 2D `std::array<std::array<W, MaxN>, MaxN>` tables even for small `n`. Refactored to 1D in-place DP using the backward-iteration pattern `dp[j] = dp[j-1] + j * dp[j]`, reducing stack to `std::array<W, MaxN + 1>`.

## Development history

### Phase A: Compile-time cached Stirling/Pascal tables

Basis conversion functions (`monomial_to_falling`, `falling_to_monomial`, `monomial_to_taylor`, `taylor_to_monomial`) and difference/shift operators (`forward_difference`, `taylor_shift`) previously constructed a new `StirlingCache<N,W>` or `PascalCache<N,W>` on every call. In `constexpr` context the compiler would evaluate this once; at runtime the O(N²) table was rebuilt each invocation.

Fixed by adding `inline constexpr` variable templates:
```cpp
template<int N, std::unsigned_integral W>
inline constexpr StirlingCache<N, W> STIRLING_CACHE{};
```
These are evaluated once per (N,W) at program initialization (compile-time folded when possible) and shared across all calls.

### Phase B: Discrete hedging example

Added `31_discrete_hedging.cpp` to dyadic-examples demonstrating the Δ/Σ operator calculus in a financial context:
- Stock price modelled as polynomial in cents
- `forward_difference` computes daily returns (ΔS)
- Second `forward_difference` gives discrete gamma (Δ²S)
- `indefinite_sum` recovers cumulative P&L (Σ = Δ⁻¹)
- Δ = e^D − I identity connects discrete hedge ratios to continuous delta

### Phase C: Auto-vectorization hints for poly_mul

The inner multiply-accumulate loop of `poly_mul_unsaturated` and the tiled `poly_mul` are textbook dense FMAs — independent iterations with no loop-carried dependencies. Added:
- `DYADIC_RESTRICT` macro (`__restrict__` on GCC/Clang) on all pointer parameters
- `#pragma GCC ivdep` on the inner convolution loops

These hints let the compiler generate SIMD (SSE/AVX2) code for runtime invocations while remaining valid in `constexpr` context (the pragma and restrict are ignored during constexpr evaluation).

### Phase D: Repo consolidation

The `dyadic-core` repository was merged into `dyadic`. Extension headers (`dynamic_polynomial.h`, `pade.h`, `continued_fractions.h`, `matrix.h`) moved from `dyadic-examples` into `dyadic/include/dyadic/`. The umbrella `dyadic.h` at the repo root includes all sub-headers. The `dyadic-examples` repo now contains only example code, with no standalone extension headers.

## Future Directions

### Wider coefficient types

The Taylor basis precision window and Witt vector overflow are both instances of the same problem: 2-adic arithmetic with 64-bit words runs out of headroom when intermediate values exceed 2^64. A natural extension is `dword_t<uint64_t>` → `detail::uint128_t` (already done) and then `W` → `detail::uint128_t` for the coefficient type itself. This would push the precision window from ~2^64/k! to ~2^128/k! — enough to cover all practical polynomial degrees.

### Ghost-map verification beyond N=2

The existing `PROOF_GHOST_HOM_U8_N2` covers all 256 values of each component for N=2, W=uint8. A natural extension verifies larger N and W combinations. The `constexpr` framework scales: the compiler will evaluate `verify_ghost_addition` for N≤8 at compile time.

### Exhaustive small-space verification (partially done)

For `uint8_t`, `PROOF_D_DELTA_U8_N2` covers ALL 65536 degree-1 polynomials for D∘Δ = Δ∘D. `PROOF_GHOST_HOM_U8_N2` covers all 256 values of each component for N=2 Witt vector ghost homomorphism. The remaining gap is full 256⁴ exhaustive ghost homomorphism (exceeds GCC's default constexpr ops budget).

### Cross-compiler support

`detail::uint128_t` already provides a compiler-abstraction layer for 128-bit arithmetic (software pair with optional `__int128` conversion). The remaining `__int128` dependency is confined to `binom()` as an optimization. Clang supports `unsigned __int128` with different constexpr evaluation limits. MSVC support for `__int128` is limited; `detail::uint128_t` works on any C++20 compiler since it is pure software.

### Hardware carry-chain acceleration

The carry chain is a linear scan — O(N) per pass. Modern x86-64 offers `addc` (add with carry) via `_addcarry_u64`, which the carry chain could use directly. The current implementation decomposes into dword accumulation then carry propagation, which is correct but not optimal. A `constexpr` path using `std::array` and a runtime path using BMI2/ADX intrinsics would cover both compile-time and maximum-throughput use cases.
