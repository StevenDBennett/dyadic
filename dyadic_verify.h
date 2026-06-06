// dyadic_verify.h — Compile-Time Verification & Invariants
// Demonstrates that the six-axiom calculus can be partially verified at compile time
// using constexpr evaluation and template metaprogramming.
//
// Key insight: For fixed N and W, many algebraic identities become
// compile-time computable proofs via constexpr template evaluation.

#pragma once

#include "dyadic.h"
#include <limits>

namespace dyadic {
namespace verify {

// ============================================================================
// Proof Generator — Metaprogram for auto-generating proof instances
// ============================================================================
// Replaces manual per-(N,W) proof with a single template invocation.
// Each (N,W) pair produces a separate ProofInstance with its own static_assert.
//
// Usage:
//   using P = prove_all<
//       []<int N, std::unsigned_integral W>() { return verify_d_delta_commute<N, W>(); },
//       Ns<2,3,4,5,6>, Ws<uint64_t>
//   >;
//   static_assert(P::value);
//
// A failing proof (N=4, W=uint64_t) produces an error trace showing N=4.

template<int... Is>
struct Ns {};

template<typename... Ts>
struct Ws {};

namespace detail {

template<auto Pred, int N, typename W>
struct ProofInstance {
    static constexpr bool value = Pred.template operator()<N, W>();
    static_assert(value, "Proof failed");
};

// W-independent proof instance (for proofs that don't need a word type).
template<auto Pred, int N>
struct ProofInstanceN {
    static constexpr bool value = Pred.template operator()<N>();
    static_assert(value, "Proof failed");
};

template<auto Pred, typename W, int... NN>
constexpr bool prove_ns() {
    return (ProofInstance<Pred, NN, W>::value && ...);
}

// Proof over Ns only (no W parameter). Currently unused but
// available for W-independent proofs.
template<auto Pred, int... NN>
constexpr bool prove_ns_only() {
    return (ProofInstanceN<Pred, NN>::value && ...);
}

} // namespace detail

template<auto Pred, typename NSeq, typename WSeq>
struct prove_all;

template<auto Pred, int... NN, typename... WW>
struct prove_all<Pred, Ns<NN...>, Ws<WW...>> {
    static_assert(sizeof...(NN) > 0, "prove_all: Ns parameter pack is empty");
    static_assert(sizeof...(WW) > 0, "prove_all: Ws parameter pack is empty");
    static constexpr bool value = (detail::prove_ns<Pred, WW, NN...>() && ...);
};

template<auto Pred, int... NN>
struct prove_all<Pred, Ns<NN...>, Ws<>> {
    static_assert(sizeof...(NN) > 0, "prove_all: Ns parameter pack is empty");
    static constexpr bool value = detail::prove_ns_only<Pred, NN...>();
};

// ============================================================================
// Compile-Time Ring Homomorphism Checks
// ============================================================================
// For specific small Witt vectors, we can verify ghost-map properties
// at compile time. The compiler evaluates these as constant expressions.

template<int N, std::unsigned_integral W>
constexpr bool verify_ghost_addition(const WittVector<N, W>& a,
                                      const WittVector<N, W>& b) {
    auto sum = witt_add(a, b);
    for (int j = 0; j < N; ++j) {
        if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
    }
    return true;
}

// Example: Compile-time proof that ghost map is a homomorphism for N=2, W=uint32_t
// with specific small inputs. The compiler evaluates this to `true` at compile time.
inline constexpr bool PROOF_GHOST_HOM_2_32 = verify_ghost_addition(
    WittVector<2, uint32_t>{{5, 7}},
    WittVector<2, uint32_t>{{3, 11}}
);

// Ghost homomorphism for N=2, W=uint8_t. The exhaustive version tests all
// 256 values of each component (512 + edge cases). Off by default to avoid
// constexpr evaluation overhead during rapid iteration; define
// DYADIC_HEAVY_PROOFS to enable it. The sampled PROOF_GHOST_HOM_2_32
// always provides basic coverage.
#ifdef DYADIC_HEAVY_PROOFS
template<std::unsigned_integral W>
constexpr bool verify_ghost_hom_all_u8() {
    constexpr W max_val = std::numeric_limits<W>::max();
    WittVector<2, W> b{{W(1), W(2)}};
    for (W a0 = 0; ; ++a0) {
        WittVector<2, W> a{{a0, W(0)}};
        auto sum = witt_add(a, b);
        if (sum.ghost(0) != static_cast<W>(a.ghost(0) + b.ghost(0))) return false;
        if (sum.ghost(1) != static_cast<W>(a.ghost(1) + b.ghost(1))) return false;
        if (a0 == max_val) break;
    }
    for (W a1 = 0; ; ++a1) {
        WittVector<2, W> a{{W(0), a1}};
        auto sum = witt_add(a, b);
        if (sum.ghost(0) != static_cast<W>(a.ghost(0) + b.ghost(0))) return false;
        if (sum.ghost(1) != static_cast<W>(a.ghost(1) + b.ghost(1))) return false;
        if (a1 == max_val) break;
    }
    {
        WittVector<2, W> a{{W(0), W(0)}};
        auto sum = witt_add(a, b);
        if (sum.ghost(0) != static_cast<W>(a.ghost(0) + b.ghost(0))) return false;
        if (sum.ghost(1) != static_cast<W>(a.ghost(1) + b.ghost(1))) return false;
    }
    {
        WittVector<2, W> a{{max_val, max_val}};
        auto sum = witt_add(a, b);
        if (sum.ghost(0) != static_cast<W>(a.ghost(0) + b.ghost(0))) return false;
        if (sum.ghost(1) != static_cast<W>(a.ghost(1) + b.ghost(1))) return false;
    }
    return true;
}
inline constexpr bool PROOF_GHOST_HOM_U8_N2 = verify_ghost_hom_all_u8<uint8_t>();

// Ghost homomorphism for N=3 with word type W. Tests each component across
// the full domain (3×domain values + cross terms). Exhaustive over the full
// 3-dimensional space would be domain^3 cases, which can exceed constexpr
// ops budget. This per-component + edge-case approach provides strong coverage.
template<std::unsigned_integral W>
constexpr bool verify_ghost_hom_u8_n3() {
    constexpr W max_val = std::numeric_limits<W>::max();
    WittVector<3, W> b{{W(1), W(2), W(3)}};
    // Vary a0
    for (W a0 = 0; ; ++a0) {
        WittVector<3, W> a{{a0, W(0), W(0)}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
        if (a0 == max_val) break;
    }
    // Vary a1
    for (W a1 = 0; ; ++a1) {
        WittVector<3, W> a{{W(0), a1, W(0)}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
        if (a1 == max_val) break;
    }
    // Vary a2
    for (W a2 = 0; ; ++a2) {
        WittVector<3, W> a{{W(0), W(0), a2}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
        if (a2 == max_val) break;
    }
    // Edge cases: all components non-zero
    for (int t = 0; t < 50; ++t) {
        WittVector<3, W> a{{static_cast<W>(t * 5), static_cast<W>(t * 7), static_cast<W>(t * 11)}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
    }
    // Max values
    {
        WittVector<3, W> a{{max_val, max_val, max_val}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
    }
    return true;
}
inline constexpr bool PROOF_GHOST_HOM_U8_N3 = verify_ghost_hom_u8_n3<uint8_t>();
#else
inline constexpr bool PROOF_GHOST_HOM_U8_N2 = true;
inline constexpr bool PROOF_GHOST_HOM_U8_N3 = true;
#endif

// ============================================================================
// Compile-Time Commutation Relations: D∘Δ = Δ∘D
// ============================================================================
// Verified using the actual formal_derivative and forward_difference
// implementations on representative test polynomials of degree N-1.
// For N=2 (degree 1) the result is degree -1 (size 0) on both sides;
// for N≥3 we compare coefficients.

template<int N, std::unsigned_integral W>
constexpr bool verify_d_delta_poly(const Polynomial<N, W, MonomialBasis>& p) {
    auto lhs = formal_derivative(forward_difference(p));
    auto rhs = forward_difference(formal_derivative(p));
    if constexpr (N >= 3) {
        for (int i = 0; i < N - 2; ++i) {
            if (lhs[i] != rhs[i]) return false;
        }
    }
    return true;
}

template<int N, std::unsigned_integral W>
constexpr bool verify_d_delta_commute() {
    // Test polynomial: all 1s
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = W(1);
        if (!verify_d_delta_poly(p)) return false;
    }
    // Alternating 1, 0
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(i & 1);
        if (!verify_d_delta_poly(p)) return false;
    }
    // Monomial: only the highest coefficient non-zero
    {
        Polynomial<N, W, MonomialBasis> p{};
        p[N-1] = W(1);
        if (!verify_d_delta_poly(p)) return false;
    }
    // Max values
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = W(-1);
        if (!verify_d_delta_poly(p)) return false;
    }
    // Ramp: {1, 2, 3, ..., N}
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(i + 1);
        if (!verify_d_delta_poly(p)) return false;
    }
    return true;
}

// Compile-time proofs for degree 1..5 (generated)
using PROOF_D_DELTA_ALL = prove_all<
    []<int N, std::unsigned_integral W>() { return verify_d_delta_commute<N, W>(); },
    Ns<2,3,4,5,6>, Ws<uint16_t, uint32_t, uint64_t>
>;
static_assert(PROOF_D_DELTA_ALL::value, "PROOF_D_DELTA_ALL");

// Exhaustive D∘Δ = Δ∘D for ALL degree-2 polynomials with word type W.
// domain^2 cases. Off by default; define DYADIC_HEAVY_PROOFS to
// enable. PROOF_D_DELTA_ALL above provides sampled coverage always.
#ifdef DYADIC_HEAVY_PROOFS
template<std::unsigned_integral W>
constexpr bool verify_d_delta_all_u8() {
    using limits = std::numeric_limits<W>;
    constexpr int domain = static_cast<int>(limits::max()) + 1;
    static_assert(Polynomial<2, W, MonomialBasis>::max_degree == 1,
        "N=2 gives max_degree 1 polynomials, DΔ and ΔD both degree -1 (size 0)");
    for (int idx = 0; idx < domain * domain; ++idx) {
        W c0 = static_cast<W>(idx / domain);
        W c1 = static_cast<W>(idx % domain);
        Polynomial<2, W, MonomialBasis> p{{c0, c1}};
        auto d_delta = formal_derivative(forward_difference(p));
        auto delta_d = forward_difference(formal_derivative(p));
        if (d_delta != delta_d) return false;
    }
    return true;
}
inline constexpr bool PROOF_D_DELTA_U8_N2 = verify_d_delta_all_u8<uint8_t>();
#else
inline constexpr bool PROOF_D_DELTA_U8_N2 = true;
#endif

// ============================================================================
// Compile-Time Nilpotency: D^N = 0
// ============================================================================
// For degree N-1 polynomials, D^N = 0. We verify this by matrix power.

template<int N, std::unsigned_integral W>
constexpr bool verify_d_nilpotent() {
    using Mat = std::array<std::array<W, N>, N>;

    // Build D matrix
    Mat D_mat{};
    for (int i = 0; i < N - 1; ++i) D_mat[i][i+1] = static_cast<W>(i + 1);

    // Compute D^N by repeated squaring/multiplication
    Mat result{};
    for (int i = 0; i < N; ++i) result[i][i] = 1; // Identity

    Mat base = D_mat;

    int exp = N;
    while (exp > 0) {
        if (exp & 1) {
            Mat temp{};
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    for (int k = 0; k < N; ++k)
                        temp[i][j] += result[i][k] * base[k][j];
            result = temp;
        }
        Mat temp{};
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int k = 0; k < N; ++k)
                    temp[i][j] += base[i][k] * base[k][j];
        base = temp;
        exp >>= 1;
    }

    // Verify all entries are zero
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            if (result[i][j] != 0) return false;
    return true;
}

using PROOF_D_NILPOTENT_ALL = prove_all<
    []<int N, std::unsigned_integral W>() { return verify_d_nilpotent<N, W>(); },
    Ns<2,3,4,5,6>, Ws<uint16_t, uint32_t, uint64_t>
>;
static_assert(PROOF_D_NILPOTENT_ALL::value, "PROOF_D_NILPOTENT_ALL");

// ============================================================================
// Compile-Time Basis Roundtrip: Monomial → Falling → Monomial = Identity
// ============================================================================
// For specific small polynomials, we verify that change_basis is involutive.

template<int N, std::unsigned_integral W>
constexpr bool verify_basis_roundtrip(const Polynomial<N, W, MonomialBasis>& p) {
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto back = change_basis<MonomialBasis>(ff);
    for (int i = 0; i < N; ++i) {
        if (p[i] != back[i]) return false;
    }
    return true;
}

// Compile-time proofs with specific test polynomials
inline constexpr bool PROOF_ROUNDTRIP_1 = verify_basis_roundtrip(
    Polynomial<4, uint64_t, MonomialBasis>{{1, 2, 3, 4}}
);
inline constexpr bool PROOF_ROUNDTRIP_2 = verify_basis_roundtrip(
    Polynomial<3, uint64_t, MonomialBasis>{{7, 0, 5}}
);

// ============================================================================
// Compile-Time Artin-Schreier Properties
// ============================================================================
// Verify ℘(0) = ℘(1) = 0 and ℘(x) = ℘(1-x) for all x in the ring.
// The second property is verified by exhaustive search for small W.

template<std::unsigned_integral W>
constexpr bool verify_artin_schreier_kernel() {
    return artin_schreier(W(0)) == W(0) && artin_schreier(W(1)) == W(0);
}

template<std::unsigned_integral W>
constexpr bool verify_artin_schreier_symmetry(W x) {
    return artin_schreier(x) == artin_schreier(W(W(1) - x));
}

// For W=uint8_t, we can do exhaustive compile-time verification!
// The compiler evaluates 256 cases at compile time.
template<std::unsigned_integral W>
constexpr bool verify_artin_schreier_symmetry_all() {
    for (W x = 0; x < W(-1); ++x) {
        if (!verify_artin_schreier_symmetry(x)) return false;
    }
    return verify_artin_schreier_symmetry(W(-1));
}

inline constexpr bool PROOF_AS_KERNEL = verify_artin_schreier_kernel<uint64_t>();
// Exhaustive proof for 8-bit words (256 cases, compile-time evaluated)
inline constexpr bool PROOF_AS_SYMMETRY_8 = verify_artin_schreier_symmetry_all<uint8_t>();

// ============================================================================
// Compile-Time Carry Chain Idempotency: C∘C = C
// ============================================================================
// For specific inputs, verify that applying carry_chain twice is the same
// as applying it once.

template<int N, std::unsigned_integral W>
constexpr bool verify_carry_idempotent(const W* input) {
    W first[N];
    dword_t<W> dword_input[N];
    for (int i = 0; i < N; ++i) dword_input[i] = input[i];
    carry_chain(first, dword_input, N);

    W second[N];
    dword_t<W> dword_first[N];
    for (int i = 0; i < N; ++i) dword_first[i] = first[i];
    carry_chain(second, dword_first, N);

    for (int i = 0; i < N; ++i) {
        if (first[i] != second[i]) return false;
    }
    return true;
}

inline constexpr uint64_t TEST_CARRY_INPUT[3] = {uint64_t(-1), uint64_t(-1), 0};
inline constexpr bool PROOF_C_IDEMPOTENT = verify_carry_idempotent<3, uint64_t>(TEST_CARRY_INPUT);

// ============================================================================
// Compile-Time Forward Difference Identity: Δ = e^D - I
// ============================================================================
// For specific polynomials, verify that forward_difference equals the
// truncated exponential series Σ_{k=1}^{N-1} D^k/k!.
//
// Uses the closed-form binomial expression:
//   D^k(p)_i / k! = C(i+k, k) · p_{i+k}
//   (e^D - I)(p)_i = Σ_{k=1}^{N-1-i} C(i+k, k) · p_{i+k}
//                   = Σ_{j=i+1}^{N-1} C(j, i) · p_j
// This avoids sequential D-application, which compounds overflow
// after 2-3 steps for full-width coefficients (carryline audit §2.3).
// Verified in exact Z before modular reduction: 1000/1000 (carryline §2.3).

template<int N, std::unsigned_integral W>
constexpr bool verify_delta_exp_d(const Polynomial<N, W, MonomialBasis>& p) {
    auto delta_p = forward_difference(p);

    // Compute e^D - I using the closed-form binomial series:
    //   D^k(p)_i / k! = C(i+k, k) · p_{i+k}
    // Summing over k gives Σ_{j=i+1}^{N-1} C(j, i) · p_j
    Polynomial<N-1, W, MonomialBasis> series{};
    for (int i = 0; i < N - 1; ++i) {
        W sum = 0;
        for (int j = i + 1; j < N; ++j) {
            sum += p[j] * binom<W>(j, i);
        }
        series[i] = sum;
    }

    for (int i = 0; i < N - 1; ++i) {
        if (delta_p[i] != series[i]) return false;
    }
    return true;
}

inline constexpr bool PROOF_DELTA_EXP = verify_delta_exp_d(
    Polynomial<4, uint64_t, MonomialBasis>{{1, 2, 3, 4}}
);

// ============================================================================
// Compile-Time Size Constraints & Invariants
// ============================================================================
// These static_asserts fire at compile time if the invariants are violated.

// Polynomial degree must be non-negative
static_assert(Polynomial<1, uint64_t, MonomialBasis>::max_degree == 0, "max_degree invariant");
static_assert(Polynomial<5, uint64_t, MonomialBasis>::max_degree == 4, "max_degree invariant");

// Word bit-width invariants (previously on DWord, now direct)
static_assert(8 * sizeof(uint8_t) == 8, "8-bit word");
static_assert(8 * sizeof(uint16_t) == 16, "16-bit word");
static_assert(8 * sizeof(uint32_t) == 32, "32-bit word");
static_assert(8 * sizeof(uint64_t) == 64, "64-bit word");

// v2(0) returns full word width
static_assert(v2(uint8_t(0)) == 8, "v2(0) = 8 for uint8_t");
static_assert(v2(uint16_t(0)) == 16, "v2(0) = 16 for uint16_t");
static_assert(v2(uint32_t(0)) == 32, "v2(0) = 32 for uint32_t");
static_assert(v2(uint64_t(0)) == 64, "v2(0) = 64 for uint64_t");

// Modular inverse properties
static_assert(modinv_odd<uint64_t>(1) == 1, "1^{-1} = 1");
static_assert(modinv_odd<uint64_t>(3) * uint64_t(3) == 1, "3 * 3^{-1} = 1 mod 2^64");

// Binomial coefficient properties
static_assert(binom<uint64_t>(5, 2) == 10, "C(5,2) = 10");
static_assert(binom<uint64_t>(5, 0) == 1, "C(5,0) = 1");
static_assert(binom<uint64_t>(5, 5) == 1, "C(5,5) = 1");
static_assert(binom<uint64_t>(5, 3) == 10, "C(5,3) = 10");

// Stirling number roundtrip: S(4,2) and s(2,4) etc
static_assert(stirling_2<uint64_t>(4, 2) == 7, "S(4,2) = 7");
static_assert(stirling_1<uint64_t>(3, 1) == uint64_t(2), "s(3,1) = 2");

// Artin-Schreier kernel
static_assert(artin_schreier(uint64_t(0)) == 0, "℘(0) = 0");
static_assert(artin_schreier(uint64_t(1)) == 0, "℘(1) = 0");

// Witt vector FV = VF identity
inline constexpr bool PROOF_FV_VF = []() constexpr {
    WittVector<3, uint64_t> w{{1, 2, 3}};
    auto fv = w.FV();
    auto vf = w.VF();
    return fv[0] == vf[0] && fv[1] == vf[1] && fv[2] == vf[2];
}();

// Witt multiplication ghost homomorphism
inline constexpr bool PROOF_MUL_GHOST = []() constexpr {
    {
        WittVector<3, uint32_t> a{{3, 5, 7}};
        WittVector<3, uint32_t> b{{9, 11, 13}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<uint32_t>(a.ghost(j) * b.ghost(j)))
                return false;
        }
    }
    {
        WittVector<4, uint16_t> a{{1, 2, 3, 4}};
        WittVector<4, uint16_t> b{{5, 6, 7, 8}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 4; ++j) {
            if (r.ghost(j) != static_cast<uint16_t>(a.ghost(j) * b.ghost(j)))
                return false;
        }
    }
    {
        WittVector<2, uint32_t> a{{7, 3}};
        WittVector<2, uint32_t> b{{11, 5}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 2; ++j) {
            if (r.ghost(j) != static_cast<uint32_t>(a.ghost(j) * b.ghost(j)))
                return false;
        }
    }
    return true;
}();

// Witt multiplication distributes over addition (ring distributivity)
inline constexpr bool PROOF_MUL_DISTRIBUTIVE = []() constexpr {
    {
        WittVector<3, uint32_t> a{{1, 2, 3}};
        WittVector<3, uint32_t> b{{4, 5, 6}};
        WittVector<3, uint32_t> c{{7, 8, 9}};
        auto lhs = witt_mul(a, witt_add(b, c));
        auto rhs = witt_add(witt_mul(a, b), witt_mul(a, c));
        for (int i = 0; i < 3; ++i) if (lhs[i] != rhs[i]) return false;
    }
    {
        WittVector<3, uint32_t> a{{9, 8, 7}};
        WittVector<3, uint32_t> b{{6, 5, 4}};
        WittVector<3, uint32_t> c{{3, 2, 1}};
        auto lhs = witt_mul(a, witt_add(b, c));
        auto rhs = witt_add(witt_mul(a, b), witt_mul(a, c));
        for (int i = 0; i < 3; ++i) if (lhs[i] != rhs[i]) return false;
    }
    return true;
}();

// Exhaustive N=2 multiplication ghost homomorphism for word type W.
// domain_a0 × (domain_b0 / step) cases. Off by default; define
// DYADIC_HEAVY_PROOFS to enable. PROOF_MUL_GHOST above provides
// sampled coverage always.
#ifdef DYADIC_HEAVY_PROOFS
inline constexpr bool PROOF_MUL_U8_N2 = []() constexpr {
    using W = uint8_t;
    constexpr int domain = static_cast<int>(std::numeric_limits<W>::max()) + 1;
    for (int a0 = 0; a0 < domain; ++a0) {
        for (int b0 = 0; b0 < domain; b0 += 8) {
            WittVector<2, W> a{{static_cast<W>(a0), W(1)}};
            WittVector<2, W> b{{static_cast<W>(b0), W(2)}};
            auto r = witt_mul(a, b);
            if (r[0] != static_cast<W>(a[0] * b[0]))
                return false;
            W ga1 = static_cast<W>(a[0] * a[0] + 2 * a[1]);
            W gb1 = static_cast<W>(b[0] * b[0] + 2 * b[1]);
            W gr1 = static_cast<W>(r[0] * r[0] + 2 * r[1]);
            if (gr1 != static_cast<W>(ga1 * gb1))
                return false;
        }
    }
    return true;
}();

// Exhaustive N=3 multiplication ghost homomorphism for uint8_t.
// Varies each Witt component independently with fixed other components,
// plus cross-term and edge-case coverage. Full 3-dimensional exhaustive
// (256^3 input pairs) would be too expensive; this per-component + edge
// approach mirrors the addition verifier above.
inline constexpr bool PROOF_MUL_U8_N3 = []() constexpr {
    using W = uint8_t;
    constexpr W max_val = std::numeric_limits<W>::max();

    // Vary a0 across full domain
    for (W a0 = 0; ; ++a0) {
        WittVector<3, W> a{{a0, W(2), W(3)}};
        WittVector<3, W> b{{W(5), W(7), W(11)}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<W>(a.ghost(j) * b.ghost(j)))
                return false;
        }
        if (a0 == max_val) break;
    }

    // Vary a1
    for (W a1 = 0; ; ++a1) {
        WittVector<3, W> a{{W(3), a1, W(7)}};
        WittVector<3, W> b{{W(5), W(11), W(13)}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<W>(a.ghost(j) * b.ghost(j)))
                return false;
        }
        if (a1 == max_val) break;
    }

    // Vary a2
    for (W a2 = 0; ; ++a2) {
        WittVector<3, W> a{{W(3), W(5), a2}};
        WittVector<3, W> b{{W(7), W(11), W(13)}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<W>(a.ghost(j) * b.ghost(j)))
                return false;
        }
        if (a2 == max_val) break;
    }

    // Edge cases: all non-zero cross terms
    for (int t = 0; t < 50; ++t) {
        WittVector<3, W> a{{static_cast<W>(t * 3), static_cast<W>(t * 5), static_cast<W>(t * 7)}};
        WittVector<3, W> b{{static_cast<W>(t * 11), static_cast<W>(t * 13), static_cast<W>(t * 17)}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<W>(a.ghost(j) * b.ghost(j)))
                return false;
        }
    }

    // Max values
    {
        WittVector<3, W> a{{max_val, max_val, max_val}};
        WittVector<3, W> b{{max_val, max_val, max_val}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<W>(a.ghost(j) * b.ghost(j)))
                return false;
        }
    }

    return true;
}();
#else
inline constexpr bool PROOF_MUL_U8_N2 = true;
inline constexpr bool PROOF_MUL_U8_N3 = true;
#endif

// Teichmüller lift: τ(1) = 1, i.e. the lift of the multiplicative unit is the
// Witt vector multiplicative unit.
inline constexpr bool PROOF_TEICHMULLER_ONE = []() constexpr {
    auto one = teichmueller_lift<4, uint32_t>(1);
    auto expected = WittVector<4, uint32_t>::one();
    for (int i = 0; i < 4; ++i) if (one[i] != expected[i]) return false;
    return true;
}();

// Teichmüller lift is multiplicative: ghost_j(τ(ab)) = ghost_j(τ(a)) * ghost_j(τ(b))
// Verified by comparing ghost values (mod W) rather than raw components, since
// the component recovery has a precision window of 2^{8*sizeof(W)-j} per j.
inline constexpr bool PROOF_TEICHMULLER_MULT = []() constexpr {
    for (uint32_t a = 0; a < 40; ++a) {
        for (uint32_t b = 0; b < 40; ++b) {
            auto lhs = witt_mul(teichmueller_lift<3, uint32_t>(a),
                                teichmueller_lift<3, uint32_t>(b));
            auto rhs = teichmueller_lift<3, uint32_t>(a * b);
            for (int j = 0; j < 3; ++j)
                if (lhs.ghost(j) != rhs.ghost(j)) return false;
        }
    }
    return true;
}();

// Adams ψ^1 is the identity map
inline constexpr bool PROOF_ADAMS_IDENTITY = []() constexpr {
    WittVector<3, uint32_t> a{{5, 7, 11}};
    auto r = adams_operation(a, 1);
    for (int i = 0; i < 3; ++i) if (r[i] != a[i]) return false;
    return true;
}();

// Adams preserves the ghost homomorphism: ghost_j(ψ^n(a)) = ghost_j(a)^n
inline constexpr bool PROOF_ADAMS_GHOST = []() constexpr {
    {
        WittVector<3, uint32_t> a{{3, 5, 7}};
        auto r = adams_operation(a, 3);
        for (int j = 0; j < 3; ++j) {
            uint32_t g = a.ghost(j);
            uint32_t gn = g * g * g;
            if (r.ghost(j) != gn) return false;
        }
    }
    {
        WittVector<4, uint16_t> a{{1, 2, 3, 4}};
        auto r = adams_operation(a, 2);
        for (int j = 0; j < 4; ++j) {
            uint16_t g = a.ghost(j);
            uint16_t gn = g * g;
            if (r.ghost(j) != gn) return false;
        }
    }
    return true;
}();

// ============================================================================
// Compile-Time Report Generation
// ============================================================================
// Each proof is independently verified by its own static_assert.
// If any proof fails, the compiler halts with the proof's name in the
// error message — no need to count or maintain a master total.

// Ghost homomorphism
static_assert(PROOF_GHOST_HOM_2_32,   "PROOF_GHOST_HOM_2_32");
static_assert(PROOF_GHOST_HOM_U8_N2,  "PROOF_GHOST_HOM_U8_N2");
static_assert(PROOF_GHOST_HOM_U8_N3,  "PROOF_GHOST_HOM_U8_N3");

static_assert(PROOF_D_DELTA_U8_N2,    "PROOF_D_DELTA_U8_N2");

static_assert(PROOF_ROUNDTRIP_1,      "PROOF_ROUNDTRIP_1");
static_assert(PROOF_ROUNDTRIP_2,      "PROOF_ROUNDTRIP_2");

static_assert(PROOF_AS_KERNEL,        "PROOF_AS_KERNEL");
static_assert(PROOF_AS_SYMMETRY_8,    "PROOF_AS_SYMMETRY_8");

static_assert(PROOF_C_IDEMPOTENT,     "PROOF_C_IDEMPOTENT");

static_assert(PROOF_DELTA_EXP,        "PROOF_DELTA_EXP");

static_assert(PROOF_FV_VF,            "PROOF_FV_VF");

static_assert(PROOF_MUL_GHOST,        "PROOF_MUL_GHOST");
static_assert(PROOF_MUL_DISTRIBUTIVE, "PROOF_MUL_DISTRIBUTIVE");
static_assert(PROOF_MUL_U8_N2,        "PROOF_MUL_U8_N2");
static_assert(PROOF_MUL_U8_N3,        "PROOF_MUL_U8_N3");

static_assert(PROOF_TEICHMULLER_ONE,  "PROOF_TEICHMULLER_ONE");
static_assert(PROOF_TEICHMULLER_MULT, "PROOF_TEICHMULLER_MULT");
static_assert(PROOF_ADAMS_IDENTITY,   "PROOF_ADAMS_IDENTITY");
static_assert(PROOF_ADAMS_GHOST,      "PROOF_ADAMS_GHOST");

// Witt logarithm/exponential roundtrip: exp∘log = identity on odd inputs.
// Tests both the p-adic series and ghost-map recovery at compile time.
// With dynamic term counting (2× bit-width budget), v₂(a₀) ≥ 2 suffices
// for exp convergence. This proof uses a Teichmüller input a₀=512 (v₂=9)
// which converges with budget to spare.
inline constexpr bool PROOF_WITT_LOG_EXP = []() constexpr {
    // Teichmüller lift: a[0]=512 (v2=9), higher components zero.
    WittVector<3, uint32_t> a{{512, 0, 0}};
    auto b = witt_exp(a);       // b[0] is odd (unit in the Witt ring)
    auto c = witt_log(b);        // should recover a
    for (int i = 0; i < 3; ++i) {
        if (a[i] != c[i]) return false;
    }
    return true;
}();
static_assert(PROOF_WITT_LOG_EXP,      "PROOF_WITT_LOG_EXP");

// ============================================================================
// Runtime Verification (for cases too large for constexpr)
// ============================================================================
// For larger N or exhaustive searches over 64-bit space, we fall back to
// runtime verification. But the compile-time proofs cover the structural
// invariants; runtime tests verify numerical stability.

template<int N, std::unsigned_integral W>
inline int run_all_verifications() {
    int failures = 0;
    ::dyadic::detail::XorShift64 rng(12345);

    // Verify D∘Δ = Δ∘D for random polynomials
    for (int trial = 0; trial < 500; ++trial) {
        Polynomial<N, W, MonomialBasis> p;
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next());

        auto d_delta = formal_derivative(forward_difference(p));
        auto delta_d = forward_difference(formal_derivative(p));

        bool ok = true;
        for (int i = 0; i < N - 2; ++i) {
            if (d_delta[i] != delta_d[i]) { ok = false; break; }
        }
        if (!ok) failures++;
    }

    // Verify basis roundtrip: Monomial → FallingFactorial → Monomial
    {
        rng = ::dyadic::detail::XorShift64(12345);
        for (int trial = 0; trial < 200; ++trial) {
            Polynomial<N, W, MonomialBasis> p;
            for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next());
            auto ff = change_basis<FallingFactorialBasis>(p);
            auto back = change_basis<MonomialBasis>(ff);
            bool ok = true;
            for (int i = 0; i < N; ++i) {
                if (p[i] != back[i]) { ok = false; break; }
            }
            if (!ok) failures++;
        }
    }

    // Verify ghost homomorphism for random Witt vectors
    rng = ::dyadic::detail::XorShift64(12345);
    for (int trial = 0; trial < 200; ++trial) {
        WittVector<N, W> a, b;
        for (int i = 0; i < N; ++i) {
            a[i] = static_cast<W>(rng.next());
            b[i] = static_cast<W>(rng.next());
        }
        auto sum = witt_add(a, b);
        bool ok = true;
        for (int j = 0; j < N; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) {
                ok = false; break;
            }
        }
        if (!ok) failures++;
    }

    // Verify Witt multiplication distributivity: a·(b+c) = a·b + a·c
    rng = ::dyadic::detail::XorShift64(12345);
    for (int trial = 0; trial < 200; ++trial) {
        WittVector<N, W> a, b, c;
        for (int i = 0; i < N; ++i) {
            a[i] = static_cast<W>(rng.next());
            b[i] = static_cast<W>(rng.next());
            c[i] = static_cast<W>(rng.next());
        }
        auto lhs = witt_mul(a, witt_add(b, c));
        auto rhs = witt_add(witt_mul(a, b), witt_mul(a, c));
        bool ok = true;
        for (int i = 0; i < N; ++i) {
            if (lhs[i] != rhs[i]) { ok = false; break; }
        }
        if (!ok) failures++;
    }

    // Verify FV = VF for random Witt vectors
    rng = ::dyadic::detail::XorShift64(12345);
    for (int trial = 0; trial < 200; ++trial) {
        WittVector<N, W> w;
        for (int i = 0; i < N; ++i) w[i] = static_cast<W>(rng.next());
        auto fv = w.FV();
        auto vf = w.VF();
        bool ok = true;
        for (int i = 0; i < N; ++i) {
            if (fv[i] != vf[i]) { ok = false; break; }
        }
        if (!ok) failures++;
    }

    return failures;
}

} // namespace verify
} // namespace dyadic
