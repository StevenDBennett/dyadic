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
| **Compile-Time Proofs** | 20+ `static_assert` proofs verifying ring axioms, basis roundtrips, D∘Δ=Δ∘D, ghost homomorphism, carry idempotence, Witt exp/log roundtrip (see `dyadic_verify.h`) |

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
// If it compiles, all 20+ proofs passed.
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
| `test_verify.cpp` | 20+ compile-time proofs + runtime verification across 9 (N,W) combos |
| `test_property.cpp` | Randomized property-based tests: 7 invariants × 10 (N,W) combos |
| `test_full.cpp` | 16 functional test groups covering the entire API surface |
| `benchmark.cpp` | Runtime benchmarks for key operations (build manually: `g++ -O2 -std=c++20 -I.. benchmark.cpp`) |

All tests pass under GCC 14+ and Clang 17+ with ASan+UBSan. CI covers GCC (light + heavy proofs), Clang, and MSVC on Windows.

## Known Limitations

- **Taylor basis roundtrip**: `T_k = k! · FF_k` wraps when `FF_k ≥ 2^W / k!`. Use small coefficients for exact roundtrips. FallingFactorial basis has no such limitation.
- **Witt precision window**: Recovery `r_j = (G_j − S_j) / 2^j` requires `r_j < 2^{W−j}`.
- **Witt exp/log term truncation**: Uses fixed 16-term series; requires `v₂(ghost_j(a)) ≥ 9` for accuracy. This holds automatically for Teichmüller lifts when `v₂(a₀) ≥ 9`. Non-Teichmüller inputs need `v₂(a_j) ≥ 9−j` for j>0.
- **`detail::uint128_t`** is a software 128-bit pair — no `unsigned __int128` required. `__int128` is used only as an optimization in `binom()`, guarded by feature-test macros.

## License

This project is made available under the terms of the MIT License.
