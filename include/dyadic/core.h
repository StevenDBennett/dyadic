// dyadic/core.h — Six-Axiom 2-Adic Operator Calculus (Core)
// The self-contained minimum: word types → carry operator → polynomials → D/Δ
//
// Six axioms:
//   1. synthetic_divide  — carry propagation as (I−N)⁻¹
//   2. formal_derivative — D
//   3. forward_difference — Δ
//   4. exact_forward_difference — exact Δ
//   5. exact_formal_derivative — exact D
//   6. [change_basis — requires dyadic/basis.h]
//
// Extensions (basis conversion, Witt vectors, compose/reversion, etc.)
// are in separate dyadic/*.h headers. This core has zero dependencies
// beyond the C++20 standard library.

#pragma once

#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <concepts>
#include <type_traits>
#include <utility>

namespace dyadic {

// ============================================================================
// 1. Type Tags and Ring Categories
// ============================================================================

struct Exact {};
template<typename W> struct Modular {};
struct Unnormalized {};

template<typename T> concept RingType = requires {
    typename T::word_type;
    { T::value } -> std::convertible_to<typename T::word_type>;
};

template<std::unsigned_integral W, typename Tag>
struct Ring {
    W value{};
    using word_type = W;
    using tag_type = Tag;
};

template<std::unsigned_integral W>
using ExactRing_t = Ring<W, Exact>;

template<std::unsigned_integral W>
using NormalizedRing_t = Ring<W, Modular<W>>;

template<std::unsigned_integral W>
using UnnormalizedRing_t = Ring<W, Unnormalized>;

// ============================================================================
// 2. Double-Width Type Mapping
// ============================================================================
// Consistent 2x mapping for all word sizes: 8->16, 16->32, 32->64, 64->128.
// For uint64_t the 128-bit pair is a software type (detail::uint128_t) that
// works on all compilers. Used by carry_chain, poly_mul_unsaturated, and
// ghost computations. The secondary widen_t<dword_t<W>> provides 4x for
// unsaturated accumulation.

namespace detail { struct uint128_t; }

template<std::unsigned_integral T>
struct double_width {};
template<> struct double_width<uint8_t>  { using type = uint16_t; };
template<> struct double_width<uint16_t> { using type = uint32_t; };
template<> struct double_width<uint32_t> { using type = uint64_t; };
template<> struct double_width<uint64_t> { using type = detail::uint128_t; };

template<std::unsigned_integral W>
using dword_t = typename double_width<W>::type;

namespace detail {

// Software 128-bit unsigned pair. Provides constexpr arithmetic for
// carry-chain, ghost-map, and polynomial accumulation. Convertible
// to/from unsigned __int128 where available (GCC/Clang).
struct uint128_t {
    static constexpr int word_bits  = 8 * sizeof(uint64_t);
    static constexpr int half_bits  = word_bits / 2;
    static constexpr int full_bits  = 2 * word_bits;
    static constexpr uint64_t low_mask = ~uint64_t(0) >> half_bits;

    uint64_t lo;
    uint64_t hi;

    constexpr uint128_t() : lo(0), hi(0) {}
    constexpr uint128_t(uint64_t v) : lo(v), hi(0) {}
    constexpr uint128_t(uint64_t hi_, uint64_t lo_) : lo(lo_), hi(hi_) {}
    constexpr explicit operator uint64_t() const { return lo; }
    template<std::unsigned_integral U>
    constexpr explicit operator U() const { return static_cast<U>(lo); }

#if defined(__SIZEOF_INT128__) && !defined(__STRICT_ANSI__)
    constexpr explicit operator unsigned __int128() const {
        return (static_cast<unsigned __int128>(hi) << word_bits) | lo;
    }
#endif

    friend constexpr uint128_t operator+(uint128_t a, uint128_t b) {
        uint128_t r;
        r.lo = a.lo + b.lo;
        r.hi = a.hi + b.hi + (r.lo < a.lo ? 1 : 0);
        return r;
    }

    friend constexpr uint128_t operator-(uint128_t a, uint128_t b) {
        uint128_t r;
        r.lo = a.lo - b.lo;
        r.hi = a.hi - b.hi - (r.lo > a.lo ? 1 : 0);
        return r;
    }

    // 64x64 -> 128 multiplication decomposed into half-width limbs.
    friend constexpr uint128_t operator*(uint128_t a, uint128_t b) {
        auto mul64 = [](uint64_t x, uint64_t y) -> uint128_t {
            uint64_t xl = x & low_mask, xh = x >> half_bits;
            uint64_t yl = y & low_mask, yh = y >> half_bits;
            uint64_t p00 = xl * yl;
            uint64_t p01 = xl * yh;
            uint64_t p10 = xh * yl;
            uint64_t p11 = xh * yh;
            uint64_t carry = (p00 >> half_bits) + (p01 & low_mask) + (p10 & low_mask);
            uint128_t r;
            r.lo = (p00 & low_mask) | ((carry & low_mask) << half_bits);
            r.hi = p11 + (p01 >> half_bits) + (p10 >> half_bits) + (carry >> half_bits);
            return r;
        };
        auto p00 = mul64(a.lo, b.lo);
        auto p01 = mul64(a.lo, b.hi);
        auto p10 = mul64(a.hi, b.lo);
        uint64_t lo = p00.lo;
        uint64_t hi = p00.hi + p01.lo + p10.lo;
        return uint128_t(hi, lo);
    }

    friend constexpr uint128_t operator<<(uint128_t a, int shift) {
        if (shift >= full_bits) return uint128_t(0, 0);
        if (shift <= 0) return a;
        if (shift >= word_bits) return uint128_t(a.lo << (shift - word_bits), 0);
        return uint128_t((a.hi << shift) | (a.lo >> (word_bits - shift)), a.lo << shift);
    }

    friend constexpr uint128_t operator>>(uint128_t a, int shift) {
        if (shift >= full_bits) return uint128_t(0, 0);
        if (shift <= 0) return a;
        if (shift >= word_bits) return uint128_t(0, a.hi >> (shift - word_bits));
        return uint128_t(a.hi >> shift, (a.lo >> shift) | (a.hi << (word_bits - shift)));
    }

    friend constexpr bool operator==(uint128_t a, uint128_t b) { return a.lo == b.lo && a.hi == b.hi; }
    friend constexpr bool operator!=(uint128_t a, uint128_t b) { return !(a == b); }
    friend constexpr bool operator<(uint128_t a, uint128_t b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
    friend constexpr bool operator>(uint128_t a, uint128_t b) { return b < a; }
    friend constexpr bool operator<=(uint128_t a, uint128_t b) { return !(b < a); }
    friend constexpr bool operator>=(uint128_t a, uint128_t b) { return !(a < b); }

    friend constexpr uint128_t operator|(uint128_t a, uint128_t b) { return {a.hi | b.hi, a.lo | b.lo}; }
    friend constexpr uint128_t operator&(uint128_t a, uint128_t b) { return {a.hi & b.hi, a.lo & b.lo}; }
    friend constexpr uint128_t operator^(uint128_t a, uint128_t b) { return {a.hi ^ b.hi, a.lo ^ b.lo}; }
    friend constexpr uint128_t operator~(uint128_t a) { return {~a.hi, ~a.lo}; }

    friend constexpr uint128_t operator/(uint128_t a, uint128_t b) {
        assert(b != uint128_t(0) && "uint128_t division by zero");
        uint128_t q = 0, r = 0;
        for (int i = full_bits - 1; i >= 0; --i) {
            r = r << 1;
            r.lo |= (a >> i).lo & 1;
            if (r >= b) { r = r - b; q.lo |= (uint64_t(1) << (i % word_bits)); if (i >= word_bits) q.hi |= uint64_t(1) << (i - word_bits); }
        }
        return q;
    }

    friend constexpr uint128_t operator%(uint128_t a, uint128_t b) {
        assert(b != uint128_t(0) && "uint128_t modulo by zero");
        uint128_t r = 0;
        for (int i = full_bits - 1; i >= 0; --i) {
            r = r << 1;
            r.lo |= (a >> i).lo & 1;
            if (r >= b) r = r - b;
        }
        return r;
    }

    static constexpr uint128_t all_ones() { return {~0ULL, ~0ULL}; }
    static constexpr int bit_width = full_bits;

    constexpr uint128_t& operator+=(uint128_t v) { *this = *this + v; return *this; }
    constexpr uint128_t& operator-=(uint128_t v) { *this = *this - v; return *this; }
    constexpr uint128_t& operator*=(uint128_t v) { *this = *this * v; return *this; }
    constexpr uint128_t& operator<<=(int s) { *this = *this << s; return *this; }
    constexpr uint128_t& operator>>=(int s) { *this = *this >> s; return *this; }
    constexpr uint128_t& operator%=(uint128_t v) { *this = *this % v; return *this; }
};

// SFINAE-friendly widen: yields double_width<T>::type when available, else T.
template<typename T, typename = void>
struct widen { using type = T; };
template<typename T>
struct widen<T, std::void_t<typename double_width<T>::type>> {
    using type = typename double_width<T>::type;
};
template<typename T> using widen_t = typename widen<T>::type;

// Compile-time all-ones mask for any integer type.
template<typename U>
struct all_ones_mask { static constexpr U value = U(-1); };
template<>
struct all_ones_mask<uint128_t> { static constexpr uint128_t value = uint128_t::all_ones(); };
template<typename U>
inline constexpr U all_ones_mask_v = all_ones_mask<U>::value;

} // namespace detail

// 4x width alias for the common widen_t<dword_t<W>> pattern.
template<std::unsigned_integral W>
using quad_width = detail::widen_t<dword_t<W>>;

// ============================================================================
// 3. 2-Adic Primitives
// ============================================================================

template<std::unsigned_integral W>
constexpr int v2(W x) noexcept {
    if (x == 0) return 8 * sizeof(W);
    return std::countr_zero(x);
}

template<std::unsigned_integral W>
constexpr W modinv_odd(W a) noexcept {
    W x = 1;
    constexpr int iterations = std::bit_width(8 * sizeof(W)) - 1;
    for (int i = 0; i < iterations; ++i) {
        x = x * (W(2) - a * x);
    }
    return x;
}

template<std::unsigned_integral W>
constexpr W exact_divide(W x, W d) noexcept {
    return x * modinv_odd(d);
}

// Exact 2-adic division by 2^k: sign-extended right shift.
template<typename U>
constexpr U div_2k_adic(U val, int k) noexcept {
    constexpr int B = 8 * sizeof(U);
    if (k <= 0) return val;
    if (k >= B) return 0;
    if ((val >> (B - 1)) == 0) return val >> k;
    return (val >> k) | (detail::all_ones_mask_v<U> << (B - k));
}

// Artin-Schreier operator P(x) = x^2 - x.
// Kernel: P(0) = P(1) = 0. Symmetry: P(x) = P(1 - x).
template<std::unsigned_integral W>
constexpr W artin_schreier(W x) noexcept {
    return x * x - x;
}

// ============================================================================
// 4. Basis Tags
// ============================================================================

struct MonomialBasis {};
struct FallingFactorialBasis {};
struct TaylorBasis {};

// ============================================================================
// 5. Combinatorial Caches
// ============================================================================

namespace detail {

template<int N, std::unsigned_integral W>
struct StirlingCache {
    std::array<std::array<W, N>, N> S2{};
    std::array<std::array<W, N>, N> s1{};

    constexpr StirlingCache() noexcept {
        if constexpr (N > 0) {
            S2[0][0] = 1;
            s1[0][0] = 1;
        }
        for (int i = 1; i < N; ++i) {
            for (int j = 1; j <= i; ++j) {
                S2[i][j] = S2[i-1][j-1] + static_cast<W>(j) * S2[i-1][j];
                s1[i][j] = s1[i-1][j-1] - static_cast<W>(i - 1) * s1[i-1][j];
            }
        }
    }
};

template<int N, std::unsigned_integral W>
struct PascalCache {
    std::array<std::array<W, N>, N> C{};

    constexpr PascalCache() noexcept {
        if constexpr (N > 0) {
            C[0][0] = 1;
        }
        for (int n = 1; n < N; ++n) {
            C[n][0] = C[n][n] = 1;
            for (int k = 1; k < n; ++k) {
                C[n][k] = C[n-1][k-1] + C[n-1][k];
            }
        }
    }
};

template<int N, std::unsigned_integral W>
inline constexpr StirlingCache<N, W> STIRLING_CACHE{};

template<int N, std::unsigned_integral W>
inline constexpr PascalCache<N, W> PASCAL_CACHE{};

} // namespace detail

// ============================================================================
// 6. Polynomial
// ============================================================================

template<int N, std::unsigned_integral W, typename Basis = MonomialBasis>
struct Polynomial : std::array<W, N> {
    using word_type = W;
    using basis_type = Basis;
    static constexpr int max_degree = N - 1;

    constexpr int actual_degree() const noexcept {
        for (int i = N - 1; i >= 0; --i) {
            if ((*this)[i] != W(0)) return i;
        }
        return -1;
    }

    constexpr Polynomial() noexcept = default;
    constexpr Polynomial(std::array<W, N> coeffs) noexcept : std::array<W, N>(coeffs) {}

    constexpr W eval(W x) const noexcept {
        static_assert(std::is_same_v<Basis, MonomialBasis>,
            "Polynomial::eval: non-monomial bases require dyadic/basis.h");
        W r = 0;
        for (int i = N - 1; i >= 0; --i) {
            r = r * x + (*this)[i];
        }
        return r;
    }

    constexpr Polynomial operator+(const Polynomial& o) const noexcept {
        Polynomial r;
        for (int i = 0; i < N; ++i) r[i] = (*this)[i] + o[i];
        return r;
    }

    constexpr Polynomial operator-(const Polynomial& o) const noexcept {
        Polynomial r;
        for (int i = 0; i < N; ++i) r[i] = (*this)[i] - o[i];
        return r;
    }
};

// ============================================================================
// 7. [change_basis — in dyadic/basis.h]
// ============================================================================

// ============================================================================
// 8. synthetic_divide (Axiom 1)
// ============================================================================

template<std::unsigned_integral W>
constexpr Ring<W, Modular<W>> synthetic_divide(Ring<W, Unnormalized> x) noexcept {
    return Ring<W, Modular<W>>{x.value};
}

template<int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N, W, Basis> synthetic_divide(Polynomial<N, W, Basis> p) noexcept {
    constexpr int shift = (8 * sizeof(W)) / 2;
    constexpr W mask = (W(1) << shift) - 1;

    W carry = 0;
    for (int i = 0; i < N; ++i) {
        W old = p[i];
        W val = old + carry;
        bool overflow = val < old;
        p[i] = val & mask;
        carry = (val >> shift) | (overflow ? (W(1) << shift) : W(0));
    }
    return p;
}

// ============================================================================
// 9. carry_normalize
// ============================================================================

template<std::unsigned_integral W>
constexpr Ring<W, Modular<W>> carry_normalize(Ring<W, Unnormalized> x) noexcept {
    return Ring<W, Modular<W>>{x.value};
}

template<std::unsigned_integral W>
constexpr Ring<W, Modular<W>> carry_normalize(Ring<W, Modular<W>> x) noexcept {
    return x;
}

template<std::unsigned_integral W>
constexpr Ring<W, Modular<W>> carry_normalize(Ring<W, Exact> x) noexcept {
    return Ring<W, Modular<W>>{x.value};
}

// ============================================================================
// 10. Formal Derivative — Monomial Basis (Axiom 2)
// ============================================================================
// FF and Taylor overloads require dyadic/basis.h

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, MonomialBasis> formal_derivative(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    Polynomial<N-1, W, MonomialBasis> r{};
    for (int i = 1; i < N; ++i) {
        r[i-1] = static_cast<W>(i) * p[i];
    }
    return r;
}

// ============================================================================
// 11. Forward Difference — Monomial Basis (Axiom 3)
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, MonomialBasis> forward_difference(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    constexpr auto& C = detail::PASCAL_CACHE<N, W>;
    Polynomial<N-1, W, MonomialBasis> r{};
    for (int i = 0; i < N - 1; ++i) {
        W sum = 0;
        for (int j = i + 1; j < N; ++j) {
            sum += p[j] * C.C[j][i];
        }
        r[i] = sum;
    }
    return r;
}

// ============================================================================
// 12. Exact Variants — Monomial Basis (Axioms 4 & 5)
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, MonomialBasis> exact_forward_difference(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    return forward_difference(p);
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, MonomialBasis> exact_formal_derivative(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    return formal_derivative(p);
}

// ============================================================================
// 13. Carry Chain
// ============================================================================

namespace detail {

// Add with carry: returns (sum, carry_out) where carry_in/carry_out are 0 or 1.
template<std::unsigned_integral W>
constexpr std::pair<W, W> adc(W a, W b, W carry_in) noexcept {
    using dw_t = dword_t<W>;
    dw_t sum = static_cast<dw_t>(a) + static_cast<dw_t>(b) + static_cast<dw_t>(carry_in);
    return {static_cast<W>(sum), static_cast<W>(sum >> (8 * sizeof(W)))};
}

// ADX-accelerated 64-bit add with carry.
#if defined(__x86_64__) && defined(__ADX__)
#if defined(_MSC_VER)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
template<>
inline std::pair<uint64_t, uint64_t> adc<uint64_t>(uint64_t a, uint64_t b, uint64_t carry_in) noexcept {
    unsigned char cf = static_cast<unsigned char>(carry_in);
    uint64_t result;
    cf = _addcarryx_u64(cf, a, b, &result);
    return {result, static_cast<uint64_t>(cf)};
}
#endif

// Overflow-checked addition (constexpr on GCC 5+/Clang 3.8+).
template<std::unsigned_integral T>
constexpr bool add_overflow(T a, T b, T* result) noexcept {
#if (defined(__GNUC__) || defined(__clang__)) && !defined(__STRICT_ANSI__)
    return __builtin_add_overflow(a, b, result);
#else
    *result = a + b;
    return *result < a;
#endif
}

// Overflow-checked multiplication (constexpr on GCC 5+/Clang 3.8+).
template<std::unsigned_integral T>
constexpr bool mul_overflow(T a, T b, T* result) noexcept {
#if (defined(__GNUC__) || defined(__clang__)) && !defined(__STRICT_ANSI__)
    return __builtin_mul_overflow(a, b, result);
#else
    *result = a * b;
    return a != 0 && *result / a != b;
#endif
}

} // namespace detail

template<std::unsigned_integral W, typename Accum = dword_t<W>>
constexpr Accum carry_chain(W* r, const Accum* v, int N) noexcept {
    Accum carry{};
    for (int i = 0; i < N; ++i) {
        Accum sum = v[i] + carry;
        r[i] = static_cast<W>(sum);
        carry = sum >> (8 * sizeof(W));
    }
    return carry;
}

template<std::unsigned_integral W>
constexpr W carry_chain_word(W hi, W lo) noexcept {
    dword_t<W> sum = (static_cast<dword_t<W>>(hi) << (8 * sizeof(W))) + lo;
    return static_cast<W>(sum);
}

// Restrict pointer annotation for auto-vectorization.
#if defined(__GNUC__) || defined(__clang__)
#define DYADIC_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define DYADIC_RESTRICT __restrict
#else
#define DYADIC_RESTRICT
#endif

// ============================================================================
// 14. Polynomial Multiplication
// ============================================================================

namespace detail {

template<std::unsigned_integral W, typename Accum = quad_width<W>>
constexpr void poly_mul_unsaturated(Accum* DYADIC_RESTRICT r,
    const W* DYADIC_RESTRICT a, int na,
    const W* DYADIC_RESTRICT b, int nb) noexcept {
    int nr = na + nb - 1;
    for (int i = 0; i < nr; ++i) r[i] = 0;
    for (int i = 0; i < na; ++i) {
#if defined(__clang__)
        #pragma clang loop vectorize(enable)
#elif defined(__GNUC__)
        #pragma GCC ivdep
#endif
        for (int j = 0; j < nb; ++j) {
            r[i + j] += static_cast<Accum>(a[i]) * static_cast<Accum>(b[j]);
        }
    }
}

} // namespace detail

template<std::unsigned_integral W>
constexpr void poly_mul(W* DYADIC_RESTRICT r, const W* DYADIC_RESTRICT a, int na,
                        const W* DYADIC_RESTRICT b, int nb) noexcept {
    using accum_t = quad_width<W>;
    constexpr int CHUNK_COUNT = 256 / sizeof(accum_t);
    int nr = na + nb - 1;
    if (nr <= CHUNK_COUNT) {
        std::array<accum_t, CHUNK_COUNT> buf{};
        detail::poly_mul_unsaturated(buf.data(), a, na, b, nb);
        carry_chain(r, buf.data(), nr);
    } else {
        std::array<accum_t, CHUNK_COUNT> buf{};
        accum_t carry{};
        int pos = 0;

        while (pos < nr) {
            int end = pos + CHUNK_COUNT;
            if (end > nr) end = nr;
            int chunk_len = end - pos;
            for (int i = 0; i < chunk_len; ++i) buf[i] = 0;

            for (int i = 0; i < na; ++i) {
                int j_start = (pos > i) ? pos - i : 0;
                int j_end = nb;
                if (i + j_end > end) j_end = end - i;
#if defined(__clang__)
                #pragma clang loop vectorize(enable)
#elif defined(__GNUC__)
                #pragma GCC ivdep
#endif
                for (int j = j_start; j < j_end; ++j) {
                    buf[i + j - pos] += static_cast<accum_t>(a[i]) * static_cast<accum_t>(b[j]);
                }
            }

            for (int i = 0; i < chunk_len; ++i) {
                accum_t sum = buf[i] + carry;
                r[pos + i] = static_cast<W>(sum);
                carry = sum >> (8 * sizeof(W));
            }
            pos = end;
        }
    }
}

template<int N, int M, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N+M-1, W, Basis>
operator*(const Polynomial<N, W, Basis>& p, const Polynomial<M, W, Basis>& q) noexcept {
    constexpr int NR = N + M - 1;
    Polynomial<NR, W, Basis> r{};
    poly_mul(r.data(), p.data(), N, q.data(), M);
    return r;
}

// Coefficient-wise polynomial multiplication (standard ring convolution).
// Use this to verify poly_divmod results (operator* is carry-chain, not standard ring).
template<int N, int M, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N+M-1, W, Basis>
poly_mul_cw(const Polynomial<N, W, Basis>& a, const Polynomial<M, W, Basis>& b) noexcept {
    constexpr int NR = N + M - 1;
    Polynomial<NR, W, Basis> r{};
    for (int i = 0; i < N; ++i) {
        if (a[i] == W(0)) continue;
        for (int j = 0; j < M; ++j) {
            r[i+j] += a[i] * b[j];
        }
    }
    return r;
}

// Verify A = Q*B + R in the coefficient-wise (standard) ring.
template<int N, int M, std::unsigned_integral W>
constexpr bool verify_divmod_cw(
    const Polynomial<N, W, MonomialBasis>& A,
    const Polynomial<N, W, MonomialBasis>& Q,
    const Polynomial<M, W, MonomialBasis>& B,
    const Polynomial<M, W, MonomialBasis>& R) noexcept
{
    auto QB = poly_mul_cw(Q, B);
    for (int i = 0; i < N; ++i) {
        W rhs = QB[i];
        if (i < M) rhs += R[i];
        if (A[i] != rhs) return false;
    }
    return true;
}

} // namespace dyadic
