# Future Investigations

Priority order (P1 = most impactful to investigate, P3 = speculative).

## P1 — Concrete Pain Points

- [ ] **Runtime-degree algorithms** — Can Witt, Padé, compose, GCD, Stirling
      be refactored to operate on a `DyadicPolynomial` concept so DynamicPolynomial
      gets full algorithm access? What breaks? What gets slower?
- [ ] **Static bound (N ≤ 32)** — How much would replacing stack arrays with
      `std::vector` in ghost_recover and pade_linear_system cost in complexity
      and constexpr reach? Can the limit be raised to 1024 without blowing
      the constexpr step budget?
- [ ] **Contract macros** — What would a `DYADIC_CONTRACT(cond, msg)` macro
      look like that asserts in debug and optionally throws `std::logic_error`
      in release? How does it interact with the existing fork-based
      precondition tests?

## P2 — Portability & Correctness

- [x] **Hardware uint128** — `#if defined(__SIZEOF_INT128__)` gives 2.7×
      speedup for `poly_mul` on uint64_t (deg 63). Implemented in `core.h`
      via conditional `accum_t` type in the `poly_mul` function. The software
      `detail::uint128_t` remains the default for all non-hot paths.
- [ ] **MSVC proof support** — Can proof instantiations be extracted into a
      single TU so MSVC doesn't need `-fconstexpr-steps`? What's the minimum
      CMake change?
- [ ] **Polynomial private inheritance** — If `: std::array<W,N>` is changed to
      private with using-declarations, what user code breaks? Is it worth the
      migration cost?
- [ ] **Per-TU proof dedup** — Would `extern template` for the 24 proofs
      meaningfully reduce compile times in multi-TU projects?

## P3 — Speculative / Enhancements

- [ ] **uint256 quad_width** — Is a software 4‑limb `uint256` type worth
      implementing so uint64_t Witt ops get a full 4× ghost window, or does
      the 2× present limitation only affect contrived edge cases?
- [ ] **Two-ring API cleanup** — Would renaming `operator*` to `operator^`
      (or adding tag-dispatched `carry_mul` / `cw_mul`) reduce confusion, or
      just break existing code for cosmetic reasons?
- [ ] **DynamicMatrix** — How much of the existing `matrix.h` API can be
      shared with a heap-allocated dynamic-size variant? Bareiss determinant?
- [ ] **Basis-aware iterators** — Could `as_monomial()`, `as_ff()`, `as_taylor()`
      view adaptors be built without runtime overhead? Where would they be
      genuinely useful?
- [ ] **SIMD poly_mul** — Does platform-intrinsic coefficient-wise
      multiplication for large N show measurable speedup, or is the
      bottleneck always the carry chain?
