# dyadic — Six-Axiom 2-Adic Operator Calculus

Single-header C++20 library for arithmetic in ℤ₂[[t]]: carry chains, formal derivatives, forward differences, Witt vectors, basis conversions, power series composition and reversion.

```cpp
#include "dyadic.h"
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
| **Compile-Time Proofs** | 22 named `static_assert` proofs (~55 total assertions) verifying ring axioms, basis roundtrips, D∘Δ=Δ∘D, ghost homomorphism, carry idempotence, Witt exp/log roundtrip (see `dyadic_verify.h`) |

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
#include "dyadic_verify.h"  // triggers all static_asserts at compile time
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

**As a single header** — copy `dyadic.h` into your project and `#include "dyadic.h"`.

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
| `test_verify.cpp` | 22 named compile-time proofs (~55 assertions) + runtime verification across 9 (N,W) combos |
| `test_property.cpp` | Randomized property-based tests: 15 invariants × 10 (N,W) combos |
| `test_full.cpp` | 17 functional test groups covering the entire API surface |
| `test_negatives.cpp` | Fork-based assertion verification (6 precondition checks) |
| `benchmark.cpp` | Runtime benchmarks for key operations (build manually: `g++ -O2 -std=c++20 -I.. benchmark.cpp`) |

All tests pass under GCC 14+ and Clang 17+ with ASan+UBSan. CI covers GCC (light + heavy proofs), Clang, and MSVC on Windows.

## Design Notes

- **Two ring semantics**: `operator*` uses carry-chain `poly_mul` (correct ℤ₂ arithmetic with overflow propagation). Functions with `_cw` suffix (`poly_divmod_cw`, `poly_gcd_cw`, `poly_mul_cw`) operate coefficient-wise (standard ring). These are **different rings** — use `poly_mul_cw()` and `verify_divmod_cw()` to verify divmod/gcd results.
- **`quad_width<W>` alias**: Shorthand for `widen_t<dword_t<W>>` (up to 4× word width), used by ghost-map accumulation and unsaturated polynomial products.
- **`assert` preconditions**: All silent-failure points (`poly_divmod_cw` with zero/even divisor, `witt_exp` with v₂(a₀) < 2, `witt_log`/`witt_inverse` with even a₀, `reversion` with even P[1]) now use `assert()`. Invalid inputs trigger `SIGABRT` in debug builds.
- **Compile-time cached Stirling/Pascal tables**: `detail::STIRLING_CACHE<N,W>` and `detail::PASCAL_CACHE<N,W>` are `inline constexpr` globals — computed once per (N,W) at compile time, shared across all basis conversion and difference calls.
- **Auto-vectorization hints**: `poly_mul_unsaturated` and `poly_mul` use `DYADIC_RESTRICT` (`__restrict__`) and `#pragma GCC ivdep` on their inner multiply-accumulate loops, enabling the compiler to generate SIMD code for runtime invocations without breaking `constexpr`.

## Known Limitations

- **`Polynomial::degree()` renamed to `max_degree()`**: The old name was misleading (it returned `N−1`, the maximum possible degree, not the actual degree). Added `actual_degree()` member function.
- **Taylor basis roundtrip**: `T_k = k! · FF_k` wraps when `FF_k ≥ 2^W / k!`. Use small coefficients for exact roundtrips. FallingFactorial basis has no such limitation.
- **Witt precision window**: Recovery `r_j = (G_j − S_j) / 2^j` requires `r_j < 2^{W−j}`.
- **Witt exp/log term truncation**: Uses valuation-aware dynamic term counting with 2× bit-width budget. Requires `v₂(x) ≥ 2` for exp convergence (mathematical limit — `v₂(x) = 1` stalls at `≤ log₂(n)+1`). Log converges for `v₂(y) ≥ 1` (~135 terms at 128-bit).
- **`detail::uint128_t`** is a software 128-bit pair — no `unsigned __int128` required. `__int128` is used only as an optimization in `binom()`, guarded by feature-test macros.

## License

This project is made available under the terms of the MIT License.
