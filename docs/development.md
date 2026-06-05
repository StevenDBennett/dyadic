# Development Notes

## What constexpr proofs caught

The compile-time proofs are not ceremonial — they found real bugs during development:

- `dyadic.h` was missing `artin_schreier()` entirely; `PROOF_AS_KERNEL` and `PROOF_AS_SYMMETRY_8` forced its addition.
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
- `std::gcd` (libstdc++) refuses to instantiate on `unsigned __int128`. The local `gcd_128` lambda works around this.
- This is not a bug in GCC or the library — it is an inherent tension between the C++ standard's unwillingness to acknowledge extended integer types and their necessity for practical 2-adic arithmetic.
- `detail::uint128_t` avoids all these friction points: the software pair satisfies `std::unsigned_integral` (`uint64_t`-based), works with `std::gcd`, and has no `-Wpedantic` warnings.

## The pattern for C-style arrays in templates

The original code used raw C arrays sized by template parameters (e.g., `W dp[64][64]`, `W powers[64][64]`). These work in practice but are fragile: they waste stack space (64×64 for small n), the bound is decoupled from the template parameter, and `std::array` provides `.size()`, iterator support, and uniform initialization. The refactor to `std::array` shrinks stack usage (e.g., `powers` went from 64×64 to N×N) and makes the intent explicit.

## O(N³) → O(N²): a case study in template metaprogramming optimization

The original Witt addition recomputed `a[i]^(2^k)` for each ghost sum `j`, yielding O(N³) total operations. The `GhostPowers<N,W>` helper precomputes `a[i]^(2^k)` once per Witt vector in O(N²), and each ghost sum then reads precomputed values in O(N). This is a natural fit for template metaprogramming: the N is known at compile time, so the precomputation buffer is a fixed-size array, and all loops are fully unrollable by the optimizer.

## Complexity wins during refactoring

Two additional hotspots were identified and fixed:

- **`taylor_to_monomial`: O(N³) → O(N²)** — factorial was recomputed from scratch (`for i=2..j`) inside the inner j-loop of a double sum. Precomputed once into a `facts[N]` array.

- **`reversion`: O(N⁴) → O(N³)** — the power matrix `R, R², ..., R^{N-1}` was fully recomputed from scratch each time a new coefficient `R[n]` was set. With incremental updates, only the O(N²) entries affected by the new `R[n]` are propagated through the power hierarchy. The m=2 case required special handling because R×R is symmetric (R[n] contributes through both factors with 2× weight).

## Future Directions

### Wider coefficient types

The Taylor basis precision window and Witt vector overflow are both instances of the same problem: 2-adic arithmetic with 64-bit words runs out of headroom when intermediate values exceed 2^64. A natural extension is `dword_t<uint64_t>` → `detail::uint128_t` (already done) and then `W` → `detail::uint128_t` for the coefficient type itself. This would push the precision window from ~2^64/k! to ~2^128/k! — enough to cover all practical polynomial degrees. However, this doubles memory and register pressure and requires a `dword_t<uint128_t>` (there is no `uint256_t` in standard C++).

### Ghost-map verification beyond N=2

The existing `PROOF_GHOST_HOM_U8_N2` covers all 256 values of each component for N=2, W=uint8. A natural extension verifies larger N and W combinations. The `constexpr` framework scales: the compiler will evaluate `verify_ghost_addition` for N≤8 at compile time.

### Exhaustive small-space verification (partially done)

For `uint8_t`, `PROOF_D_DELTA_U8_N2` covers ALL 65536 degree-1 polynomials for D∘Δ = Δ∘D. `PROOF_GHOST_HOM_U8_N2` covers all 256 values of each component for N=2 Witt vector ghost homomorphism. The remaining gap is full 256⁴ exhaustive ghost homomorphism (exceeds GCC's default constexpr ops budget).

### Cross-compiler support

`detail::uint128_t` already provides a compiler-abstraction layer for 128-bit arithmetic (software pair with optional `__int128` conversion). The remaining `__int128` dependency is confined to `binom()` as an optimization. Clang supports `unsigned __int128` with different constexpr evaluation limits. MSVC support for `__int128` is limited; `detail::uint128_t` works on any C++20 compiler since it is pure software.

### Runtime-parameterized N

All polynomials and Witt vectors are templated on N. For applications needing variable-degree polynomials, a type-erased or heap-allocated variant (e.g., `DynamicPolynomial<W>`) would trade compile-time optimization for runtime flexibility. The carry chain, basis conversions, and Witt addition all work with runtime N trivially — only the compile-time proofs depend on template N.

### Hardware carry-chain acceleration

The carry chain is a linear scan — O(N) per pass. Modern x86-64 offers `addc` (add with carry) via `_addcarry_u64`, which the carry chain could use directly. The current implementation decomposes into dword accumulation then carry propagation, which is correct but not optimal. A `constexpr` path using `std::array` and a runtime path using BMI2/ADX intrinsics would cover both compile-time and maximum-throughput use cases.

### Witt vector logarithm and exponential

The library now includes both Witt addition and multiplication, completing the full ring structure. The multiplication uses the same ghost-map-and-recovery strategy: compute the componentwise product of ghost values (via `ghost_j_dword` for exact arithmetic) and recover Witt components by Newton iteration. Ring axioms (distributivity, associativity, commutativity, identity) are verified by compile-time proofs and runtime property tests. Future extensions could include the logarithm and exponential maps, and full Teichmüller lift theory.
