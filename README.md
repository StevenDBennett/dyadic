# dyadic — Six-Axiom 2-Adic Operator Calculus

Header-only C++20 library for arithmetic in ℤ₂[[t]]: carry chains, formal derivatives, forward differences, Witt vectors, basis conversions, power series composition and reversion. Umbrella header `dyadic.h` includes all sub-headers from `include/dyadic/`.

```cpp
#include <dyadic.h>
#include <cstdio>
using namespace dyadic;

int main() {
    Polynomial<4, uint64_t> p{{1, 2, 3, 4}};

    // Formal derivative: D(P)(t) = dP/dt
    auto dp = formal_derivative(p);     // {2, 6, 12, 0}
    for (auto c : dp) std::printf("%lu ", c);

    // Forward difference: Δ(P)(t) = P(t+1) − P(t)
    auto fp = forward_difference(p);    // {9, 18, 12, 0}

    // Basis conversions (exact roundtrip)
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto back = change_basis<MonomialBasis>(ff);  // == p

    // Evaluate P(5) = 1 + 2·5 + 3·25 + 4·125 = 586
    auto y = p.eval(5);

    // Witt vector ring operations
    WittVector<3, uint64_t> a{{1, 2, 3}}, b{{4, 5, 6}};
    auto sum  = a + b;    // Witt addition (via ghost map)
    auto prod = a * b;    // Witt multiplication

    // Power series reversion (Lagrange inversion)
    Polynomial<6, uint64_t> P{{0, 1, 1}};  // t + t²
    auto Q = reversion(P);
    // Q = t − t² + 2t³ − 5t⁴ + 14t⁵ − 42t⁶ + …
    // (stored as unsigned ℤ₂ values: −1 ≡ 2⁶⁴−1, etc.)
}
```

## Requirements

- **C++20** compiler (GCC 12+, Clang 17+, MSVC 2022+)
- No external dependencies
- `detail::uint128_t` (software 128-bit pair) provides `uint64_t` word support without `unsigned __int128`
- Clang needs `-fconstexpr-steps=50000000` for compile-time proofs; GCC needs `-fconstexpr-ops-limit=200000000`

## Features

| Area | What |
|------|------|
| **2-Adic Primitives** | `v2` (valuation), `modinv_odd`, `div_2k_adic`, `artin_schreier` (℘(x)=x²−x) |
| **Polynomials** | `Polynomial<N,W,Basis>` with `eval`, `+`, `−`, `*`, basis conversion (Monomial / FallingFactorial / Taylor), GCD, resultant, discriminant, pseudo-remainder, divmod, square-free check |
| **Calculus** | `formal_derivative` (D), `forward_difference` (Δ), `taylor_shift`, `indefinite_sum` (Σ = Δ⁻¹) |
| **Witt Vectors** | `WittVector<N,W>` with ghost map, Frobenius, Verschiebung, `+`, `*`, `exp`, `log`, `inverse`, `adams_operation`, `teichmueller_lift` |
| **Carry Chain** | Full-width carry propagation `C = (I−N)⁻¹` — converges in one pass |
| **Compose / Reversion** | Power series composition P(Q(t)) and Lagrange inversion |
| **Compile-Time Proofs** | 22 named `static_assert` proofs (~55 total assertions) verifying ring axioms, basis roundtrips, D∘Δ=Δ∘D, ghost homomorphism, carry idempotence, Witt exp/log roundtrip (see `dyadic/verify.h`) |
| **Combinatorial** | `binom`, `stirling_2`, `stirling_1`, `stirling_1_unsigned` — all `constexpr`, cached at compile time |

## More Examples

### Witt vector ring

```cpp
WittVector<3, uint64_t> a{{1, 2, 3}}, b{{4, 5, 6}};

auto sum  = a + b;
auto prod = a * b;

// Frobenius and Verschiebung
auto fa = a.frobenius();
auto vb = b.verschiebung();
// F(V(a)) == V(F(a))  — verified at compile time

// Adams operation ψⁿ
auto psi = adams_operation(a, 3);

// Teichmüller lift: τ(x) = (x, 0, 0, …)  (N is required)
auto tau = teichmueller_lift<4>(uint64_t{7});
// ghost_j(τ(ab)) = ghost_j(τ(a)) · ghost_j(τ(b))  — verified at compile time
```

### Polynomial GCD, divmod, and ring semantics

```cpp
// NOTE: operator* uses carry-chain (2-adic ring).
// poly_divmod_cw / poly_gcd_cw use coefficient-wise (standard ring).
// Use poly_mul_cw() and verify_divmod_cw() for verification.

Polynomial<3, uint64_t> A{{1, 2, 1}};  // (x+1)²
Polynomial<2, uint64_t> B{{1, 1}};     // (x+1)

auto [Q, R] = poly_divmod_cw(A, B);       // Q=(x+1), R=0
bool ok = verify_divmod_cw(A, Q, B, R);   // true: A=Q·B+R ✓

auto G = poly_gcd_cw(A, B);               // G=(x+1) — divides both
```

### Polynomial composition and reversion

```cpp
Polynomial<6, uint64_t> P{{0, 1, 1}};  // P(t) = t + t²
Polynomial<6, uint64_t> Q{{0, 1, 2}};  // Q(t) = t + 2t²

auto PQ = compose(P, Q);  // P(Q(t)) — degree (N-1)*(M-1)+1 = 26
auto R  = reversion(P);   // Lagrange inverse: P(R(t)) = t
// R(t) = t − t² + 2t³ − 5t⁴ + 14t⁵ − 42t⁶ + …
// (stored in ℤ₂: negative coefficients wrap modulo 2^W)
```

### Compile-time verification

```cpp
#include <dyadic/verify.h>  // triggers all static_asserts at compile time
// If it compiles, all 22 named static_assert proofs passed.
```

Day-to-day compiles include a sampled subset. For the full exhaustive suite (256² cases, 8K+ multiplication cases), define `DYADIC_HEAVY_PROOFS`.

### Precision window checkers

```cpp
Polynomial<6, uint8_t> p{{0, 0, 0, 0, 0, 255}};
if (!check_taylor_roundtrip_precision(p))
    // T₅ = 5! · 255 = 30600 > 256 — high bits lost

WittVector<4, uint32_t> w{{1, 2, 3, 4}};
if (!check_witt_recovery_precision(w))
    // Some rⱼ ≥ 2³²⁻ʲ — ghost recovery inexact
```

### Documentation

| File | What |
|------|------|
| [`docs/precision.md`](docs/precision.md) | Five precision windows unified under a single principle |
| [`docs/theory.md`](docs/theory.md) | Four-domain unified theory: operator calculus ↔ Witt ghosts ↔ Mersenne ghosts ↔ thermodynamic classification |

## Build & Integrate

**As a header collection** — copy `dyadic.h` and the `include/dyadic/` directory into your project's include path, then `#include <dyadic.h>` (umbrella) or `#include <dyadic/core.h>` (individual header).

**With CMake**:
```bash
cmake -B build -DDYADIC_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build
```

Optional: `-DDYADIC_HEAVY_PROOFS=ON` for exhaustive compile-time proofs.

**Via CI script**:
```bash
./ci_compile_check.sh
```

## Tests

| File | What |
|------|------|
| `test_core.cpp` | Core axiom unit tests: six axioms, carry chain, arithmetic, evaluation |
| `test_verify.cpp` | 24 named compile-time proofs (~55 assertions) + runtime verification across 9 (N,W) combos |
| `test_property.cpp` | Randomized property-based tests: 18 invariants × 10 (N,W) combos |
| `test_full.cpp` | 17 functional test groups covering the entire API surface |
| `test_negatives.cpp` | Fork-based assertion verification (6 precondition checks) |
| `benchmark.cpp` | Runtime benchmarks for key operations (build manually: `g++ -O2 -std=c++20 -I. -Iinclude test/benchmark.cpp`) |

All tests pass under GCC 14+ and Clang 17+ with ASan+UBSan. CI covers GCC (light + heavy proofs), Clang, and MSVC on Windows.

## Design Notes

- **Two ring semantics**: `operator*` uses carry-chain `poly_mul` (correct ℤ₂ arithmetic with overflow propagation). Functions with `_cw` suffix (`poly_divmod_cw`, `poly_gcd_cw`, `poly_mul_cw`) operate coefficient-wise (standard ring). These are **different rings** — use `poly_mul_cw()` and `verify_divmod_cw()` to verify divmod/gcd results.
- **`quad_width<W>` alias**: Shorthand for `widen_t<dword_t<W>>` (up to 4× word width), used by ghost-map accumulation and unsaturated polynomial products.
- **`assert` preconditions**: All silent-failure points (`poly_divmod_cw` with zero/even divisor, `witt_exp` with v₂(a₀) < 2, `witt_log`/`witt_inverse` with even a₀, `reversion` with even P[1]) now use `assert()`. Invalid inputs trigger `SIGABRT` in debug builds. A few functions return 0 for invalid inputs by design rather than asserting: `poly_discriminant_cw` returns 0 when degree < 1 or leading coefficient is even; `det_laplace` returns 0 for dim > 6 (use Bareiss via `determinant_cw` instead).
- **Compile-time cached Stirling/Pascal tables**: `detail::STIRLING_CACHE<N,W>` and `detail::PASCAL_CACHE<N,W>` are `inline constexpr` globals — computed once per (N,W) at compile time, shared across all basis conversion and difference calls.
- **Auto-vectorization hints**: `poly_mul_unsaturated` and `poly_mul` use `DYADIC_RESTRICT` (`__restrict__`) and `#pragma GCC ivdep` on their inner multiply-accumulate loops, enabling the compiler to generate SIMD code for runtime invocations without breaking `constexpr`.

## Complexity

See [`docs/complexity.md`](docs/complexity.md) for a complete function-by-function breakdown.

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `v2`, `modinv_odd`, `exact_divide`, `div_2k_adic` | O(1) | O(1) | Single-instruction or fixed-iteration 2-adic primitives |
| `eval` (Monomial) | O(N) | O(1) | Horner's method |
| `eval` (FF / Taylor) | O(N) | O(1) | Incremental falling product / binomial coefficient |
| `formal_derivative` (Monomial) | O(N) | O(N) | Direct coefficient scaling |
| `formal_derivative` (FF / Taylor) | O(N²) | O(N) | Convert → monomial D → convert back |
| `forward_difference` (FF / Taylor) | O(N) | O(N) | FF diagonalizes Δ: scale+shift; Taylor: pure shift |
| `forward_difference` (Monomial) | O(N²) | O(N) | Pascal-cache binomial transform |
| `taylor_shift` | O(N²) | O(N) | Pascal-cache / falling-factorial binomial transform |
| `indefinite_sum` (FF) | O(N) | O(N) | Σ(FFₙ₋₁) = FFₙ / n |
| `indefinite_sum` (Monomial / Taylor) | O(N²) | O(N) | Converts to FF, sums, converts back |
| `poly_exp` / `poly_log` | O(N²) | O(N) | Double-loop recurrence with dword intermediates |
| `poly_mul` / `operator*` | O(N·M) | O(1) chunked | Naive convolution + carry chain; chunked 256-bit buffer |
| `poly_mul_cw` | O(N·M) | O(N+M) | Standard coefficient-wise convolution |
| `compose` | O(N²·M²) | O(N·M) | Naive power-series accumulation; result ≤ 4095 coeffs |
| `reversion` | O(N³) | O(N²) | Incremental Lagrange inversion (not Newton) |
| `change_basis` | O(N²) | O(N) | Stirling-cache triangular matrix multiply; 1–2 conversions per call |
| `witt_add` / `witt_mul` | O(N²) | O(N²) | Ghost-map + Newton recovery; `quad_width` accumulators |
| `witt_log` / `witt_exp` | O(N² + N·T) | O(N²) | T ≈ 2·bit_width (ghost p-adic series terms) |
| `witt_inverse` | O(N²) | O(N²) | Ghost inverse via `modinv_odd` per component |
| `adams_operation` | O(N² + N·n) | O(N²) | Ghost ^ n powering |
| `poly_divmod_cw` | O(N·(N+M)) | O(max(N,M)) | Long division (requires odd lc) |
| `pseudo_remainder_cw` | O(N·(N+M)) | O(N) | Subresultant PRS (no lc requirement) |
| `poly_gcd_cw` | O(K³) worst | O(K) | Euclidean PRS; K = max(N,M) |
| `polynomial_resultant_cw` | O(1) bounded | O(1) | Sylvester matrix (dim ≤ 6) + Laplace det (≤ 6! = 720) |
| `poly_discriminant_cw` | O(N) + O(1) | O(N) | Derivative + resultant |
| `poly_is_square_free_cw` | O(N³) | O(N) | gcd(P, P′) |
| `Matrix::operator*` (M×N × N×P) | O(M·N·P) | O(M·P) | i-k-j loop with zero-skip |
| `Matrix::determinant` | O(M³) | O(M·N) | Bareiss fraction-free elimination |
| `Matrix::inverse` / `solve` / `rref` | O(M³) | O(M²) | Gauss–Jordan with odd-pivot `modinv_odd` |
| `pade_approximant` [M/N] | O(N³ + M·N) | O(N² + M+N) | Gaussian elimination on N×N Toeplitz system |
| `cf_convergent` | O(n²) | O(n) | Classical recurrence; degree grows linearly |
| `div_unsigned` | O(NL·NR) | O(NL+NR) | Knuth Algorithm D (NR ≥ 2) or `divmod_single` (NR=1) |
| `reciprocal_newton` | O(NR²·log NR) | O(NR) | Newton iteration; ceil(log₂((NR+1)·W)) iterations |
| `div_newton` | O(NR²·log NR + NL·NR) | O(NL+NR) | Newton reciprocal division with correction step |
| `binom` / `stirling_2` / `stirling_1` | O(k·log n) / O(n²) | O(1) / O(n) | GCD-multiplicative / DP recurrence |
| `StirlingCache` / `PascalCache` ctor | O(N²) compile-time | O(N²) | `inline constexpr` cached once per (N,W) |

All operations are constexpr; runtime performance matches compile-time complexity bounds. The carry-chain `poly_mul` is auto-vectorized with SIMD hints (`#pragma GCC ivdep`, `DYADIC_RESTRICT`).

## Known Limitations

- **`Polynomial::degree()` renamed to `max_degree()`**: The old name was misleading (it returned `N−1`, the maximum possible degree, not the actual degree). Added `actual_degree()` member function.
- **Taylor basis roundtrip**: `T_k = k! · FF_k` wraps when `FF_k ≥ 2^W / k!`. Use small coefficients for exact roundtrips. FallingFactorial basis has no such limitation.
- **Witt precision window**: Recovery `r_j = (G_j − S_j) / 2^j` requires `r_j < 2^{W−j}`.
- **Witt exp/log term truncation**: Uses valuation-aware dynamic term counting with 2× bit-width budget. Requires `v₂(x) ≥ 2` for exp convergence (mathematical limit — `v₂(x) = 1` stalls at `≤ log₂(n)+1`). Log converges for `v₂(y) ≥ 1` (~135 terms at 128-bit).
- **`detail::uint128_t`** is a software 128-bit pair — no `unsigned __int128` required. `__int128` is used as an optimization in `binom()` (`dyadic/combinatorial.h`) and `poly_mul()` (`dyadic/core.h`, ~2.7× speedup at deg=63 for uint64_t), both guarded by `__SIZEOF_INT128__`.

## License

This project is made available under the terms of the MIT License.
