// dyadic/verify.h — Compile-Time Verification & Invariants
// Demonstrates that the six-axiom calculus can be partially verified at
// compile time using constexpr evaluation and template metaprogramming.

#pragma once

#include <dyadic/core.h>
#include <dyadic/basis.h>
#include <dyadic/witt.h>
#include <dyadic/prng.h>
#include <dyadic/combinatorial.h>

#include <limits>

namespace dyadic {
namespace verify {

// ============================================================================
// Taylor basis roundtrip precision check
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr bool check_taylor_roundtrip_precision(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    auto taylor = change_basis<TaylorBasis>(p);
    auto back = change_basis<MonomialBasis>(taylor);
    for (int i = 0; i < N; ++i) {
        if (p[i] != back[i]) return false;
    }
    return true;
}

// ============================================================================
// Proof Generator — Metaprogram for auto-generating proof instances
// ============================================================================

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

template<auto Pred, int N>
struct ProofInstanceN {
    static constexpr bool value = Pred.template operator()<N>();
    static_assert(value, "Proof failed");
};

template<auto Pred, typename W, int... NN>
constexpr bool prove_ns() {
    return (ProofInstance<Pred, NN, W>::value && ...);
}

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

template<int N, std::unsigned_integral W>
constexpr bool verify_ghost_addition(const WittVector<N, W>& a,
                                      const WittVector<N, W>& b) {
    auto sum = witt_add(a, b);
    for (int j = 0; j < N; ++j) {
        if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
    }
    return true;
}

inline constexpr bool PROOF_GHOST_HOM_2_32 = verify_ghost_addition(
    WittVector<2, uint32_t>{{5, 7}},
    WittVector<2, uint32_t>{{3, 11}}
);

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

template<std::unsigned_integral W>
constexpr bool verify_ghost_hom_u8_n3() {
    constexpr W max_val = std::numeric_limits<W>::max();
    WittVector<3, W> b{{W(1), W(2), W(3)}};
    for (W a0 = 0; ; ++a0) {
        WittVector<3, W> a{{a0, W(0), W(0)}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
        if (a0 == max_val) break;
    }
    for (W a1 = 0; ; ++a1) {
        WittVector<3, W> a{{W(0), a1, W(0)}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
        if (a1 == max_val) break;
    }
    for (W a2 = 0; ; ++a2) {
        WittVector<3, W> a{{W(0), W(0), a2}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
        if (a2 == max_val) break;
    }
    for (int t = 0; t < 50; ++t) {
        WittVector<3, W> a{{static_cast<W>(t * 5), static_cast<W>(t * 7), static_cast<W>(t * 11)}};
        auto sum = witt_add(a, b);
        for (int j = 0; j < 3; ++j) {
            if (sum.ghost(j) != static_cast<W>(a.ghost(j) + b.ghost(j))) return false;
        }
    }
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
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = W(1);
        if (!verify_d_delta_poly(p)) return false;
    }
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(i & 1);
        if (!verify_d_delta_poly(p)) return false;
    }
    {
        Polynomial<N, W, MonomialBasis> p{};
        p[N-1] = W(1);
        if (!verify_d_delta_poly(p)) return false;
    }
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = W(-1);
        if (!verify_d_delta_poly(p)) return false;
    }
    {
        Polynomial<N, W, MonomialBasis> p{};
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(i + 1);
        if (!verify_d_delta_poly(p)) return false;
    }
    return true;
}

using PROOF_D_DELTA_ALL = prove_all<
    []<int N, std::unsigned_integral W>() { return verify_d_delta_commute<N, W>(); },
    Ns<2,3,4,5,6>, Ws<uint16_t, uint32_t, uint64_t>
>;
static_assert(PROOF_D_DELTA_ALL::value, "PROOF_D_DELTA_ALL");

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

template<int N, std::unsigned_integral W>
constexpr bool verify_d_nilpotent() {
    using Mat = std::array<std::array<W, N>, N>;

    Mat D_mat{};
    for (int i = 0; i < N - 1; ++i) D_mat[i][i+1] = static_cast<W>(i + 1);

    Mat result{};
    for (int i = 0; i < N; ++i) result[i][i] = 1;

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
// Compile-Time Basis Roundtrip: Monomial -> Falling -> Monomial = Identity
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr bool verify_basis_roundtrip(const Polynomial<N, W, MonomialBasis>& p) {
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto back = change_basis<MonomialBasis>(ff);
    for (int i = 0; i < N; ++i) {
        if (p[i] != back[i]) return false;
    }
    return true;
}

inline constexpr bool PROOF_ROUNDTRIP_1 = verify_basis_roundtrip(
    Polynomial<4, uint64_t, MonomialBasis>{{1, 2, 3, 4}}
);
inline constexpr bool PROOF_ROUNDTRIP_2 = verify_basis_roundtrip(
    Polynomial<3, uint64_t, MonomialBasis>{{7, 0, 5}}
);

// ============================================================================
// Compile-Time Artin-Schreier Properties
// ============================================================================

template<std::unsigned_integral W>
constexpr bool verify_artin_schreier_kernel() {
    return artin_schreier(W(0)) == W(0) && artin_schreier(W(1)) == W(0);
}

template<std::unsigned_integral W>
constexpr bool verify_artin_schreier_symmetry(W x) {
    return artin_schreier(x) == artin_schreier(W(W(1) - x));
}

template<std::unsigned_integral W>
constexpr bool verify_artin_schreier_symmetry_all() {
    for (W x = 0; x < W(-1); ++x) {
        if (!verify_artin_schreier_symmetry(x)) return false;
    }
    return verify_artin_schreier_symmetry(W(-1));
}

inline constexpr bool PROOF_AS_KERNEL = verify_artin_schreier_kernel<uint64_t>();
inline constexpr bool PROOF_AS_SYMMETRY_8 = verify_artin_schreier_symmetry_all<uint8_t>();

// ============================================================================
// Compile-Time Carry Chain Idempotency: C∘C = C
// ============================================================================

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

template<int N, std::unsigned_integral W>
constexpr bool verify_delta_exp_d(const Polynomial<N, W, MonomialBasis>& p) {
    auto delta_p = forward_difference(p);

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

static_assert(Polynomial<1, uint64_t, MonomialBasis>::max_degree == 0, "max_degree invariant");
static_assert(Polynomial<5, uint64_t, MonomialBasis>::max_degree == 4, "max_degree invariant");

static_assert(8 * sizeof(uint8_t) == 8, "8-bit word");
static_assert(8 * sizeof(uint16_t) == 16, "16-bit word");
static_assert(8 * sizeof(uint32_t) == 32, "32-bit word");
static_assert(8 * sizeof(uint64_t) == 64, "64-bit word");

static_assert(v2(uint8_t(0)) == 8, "v2(0) = 8 for uint8_t");
static_assert(v2(uint16_t(0)) == 16, "v2(0) = 16 for uint16_t");
static_assert(v2(uint32_t(0)) == 32, "v2(0) = 32 for uint32_t");
static_assert(v2(uint64_t(0)) == 64, "v2(0) = 64 for uint64_t");

static_assert(modinv_odd<uint64_t>(1) == 1, "1^{-1} = 1");
static_assert(modinv_odd<uint64_t>(3) * uint64_t(3) == 1, "3 * 3^{-1} = 1 mod 2^64");

static_assert(binom<uint64_t>(5, 2) == 10, "C(5,2) = 10");
static_assert(binom<uint64_t>(5, 0) == 1, "C(5,0) = 1");
static_assert(binom<uint64_t>(5, 5) == 1, "C(5,5) = 1");
static_assert(binom<uint64_t>(5, 3) == 10, "C(5,3) = 10");

static_assert(stirling_2<uint64_t>(4, 2) == 7, "S(4,2) = 7");
static_assert(stirling_1<uint64_t>(3, 1) == uint64_t(2), "s(3,1) = 2");

static_assert(artin_schreier(uint64_t(0)) == 0, "P(0) = 0");
static_assert(artin_schreier(uint64_t(1)) == 0, "P(1) = 0");

// ============================================================================
// Witt vector proofs
// ============================================================================

inline constexpr bool PROOF_FV_VF = []() constexpr {
    WittVector<3, uint64_t> w{{1, 2, 3}};
    auto fv = w.FV();
    auto vf = w.VF();
    return fv[0] == vf[0] && fv[1] == vf[1] && fv[2] == vf[2];
}();

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

inline constexpr bool PROOF_MUL_U8_N3 = []() constexpr {
    using W = uint8_t;
    constexpr W max_val = std::numeric_limits<W>::max();

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

    for (int t = 0; t < 50; ++t) {
        WittVector<3, W> a{{static_cast<W>(t * 3), static_cast<W>(t * 5), static_cast<W>(t * 7)}};
        WittVector<3, W> b{{static_cast<W>(t * 11), static_cast<W>(t * 13), static_cast<W>(t * 17)}};
        auto r = witt_mul(a, b);
        for (int j = 0; j < 3; ++j) {
            if (r.ghost(j) != static_cast<W>(a.ghost(j) * b.ghost(j)))
                return false;
        }
    }

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

inline constexpr bool PROOF_TEICHMULLER_ONE = []() constexpr {
    auto one = teichmueller_lift<4, uint32_t>(1);
    auto expected = WittVector<4, uint32_t>::one();
    for (int i = 0; i < 4; ++i) if (one[i] != expected[i]) return false;
    return true;
}();

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

inline constexpr bool PROOF_ADAMS_IDENTITY = []() constexpr {
    WittVector<3, uint32_t> a{{5, 7, 11}};
    auto r = adams_operation(a, 1);
    for (int i = 0; i < 3; ++i) if (r[i] != a[i]) return false;
    return true;
}();

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

inline constexpr bool PROOF_WITT_LOG_EXP = []() constexpr {
    WittVector<5, uint32_t> a{{512, 0, 0, 0, 0}};
    auto b = witt_exp(a);
    auto c = witt_log(b);
    for (int i = 0; i < 5; ++i) {
        if (a[i] != c[i]) return false;
    }
    return true;
}();
static_assert(PROOF_WITT_LOG_EXP,      "PROOF_WITT_LOG_EXP");

inline constexpr bool PROOF_WITT_INVERSE = []() constexpr {
    for (uint32_t a0 : {1u, 3u, 5u, 7u, 9u, 65535u}) {
        WittVector<3, uint32_t> w{{a0, 2, 3}};
        if (!check_witt_inverse(w)) return false;
    }
    return true;
}();
static_assert(PROOF_WITT_INVERSE,       "PROOF_WITT_INVERSE");

inline constexpr bool PROOF_WITT_EXP_LOG_ROUNDTRIP = []() constexpr {
    for (uint32_t a0 : {1u, 5u, 9u, 13u, 65533u}) {
        WittVector<3, uint32_t> w{{a0, 2, 3}};
        if (!check_witt_exp_log_roundtrip(w)) return false;
    }
    return true;
}();
static_assert(PROOF_WITT_EXP_LOG_ROUNDTRIP, "PROOF_WITT_EXP_LOG_ROUNDTRIP");

// ============================================================================
// Runtime Verification
// ============================================================================

template<int N, std::unsigned_integral W>
inline int run_all_verifications(uint64_t base_seed = 0xDEADBEEF) {
    int failures = 0;
    auto section_seed = [base_seed](uint64_t tag) {
        return base_seed ? (base_seed ^ tag) : tag;
    };

    {
    ::dyadic::detail::XorShift64 rng(section_seed(0xD1FFACE1));
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
    }

    {
    ::dyadic::detail::XorShift64 rng(section_seed(0xCAFEBABE));
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

    {
    ::dyadic::detail::XorShift64 rng(section_seed(0xFEEDFACE));
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
    }

    {
    ::dyadic::detail::XorShift64 rng(section_seed(0xD15C0DE5));
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
    }

    {
    ::dyadic::detail::XorShift64 rng(section_seed(0xF1FFE));
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
    }

    return failures;
}

} // namespace verify
} // namespace dyadic
