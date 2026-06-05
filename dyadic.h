// dyadic — Six-Axiom 2-Adic Operator Calculus
// Single-header C++20 library — Sections 1-23

#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <concepts>
#include <cstdio>

namespace dyadic {

// ============================================================================
// 1. Type Tags and Ring Categories
// ============================================================================

struct Exact {};
template<typename W> struct Modular {};
struct Unnormalized {};

template<typename T> struct RingCategory;

template<typename T> concept IsExact = std::is_same_v<typename RingCategory<T>::type, Exact>;
template<typename T> concept IsUnnormalized = std::is_same_v<typename RingCategory<T>::type, Unnormalized>;
template<typename T> concept IsModular = std::is_same_v<typename RingCategory<T>::type, Modular<typename T::word_type>>;
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

template<std::unsigned_integral W, typename Tag>
struct RingCategory<Ring<W, Tag>> {
    using type = Tag;
};

template<std::unsigned_integral W>
using ExactRing_t = Ring<W, Exact>;

template<std::unsigned_integral W>
using NormalizedRing_t = Ring<W, Modular<W>>;

template<std::unsigned_integral W>
using UnnormalizedRing_t = Ring<W, Unnormalized>;

// ============================================================================
// 2. double_width / dword_t — Generic Double-Width Type Mapping
// ============================================================================
// Consistent 2× mapping for all word sizes: 8→16, 16→32, 32→64, 64→128.
// For uint64_t the 128-bit pair is a software type (detail::uint128_t) that
// works on all compilers; on GCC/Clang with __int128 it also supports
// implicit conversion for __int128-optimized paths.
// Used by carry_chain, poly_mul_unsaturated, and ghost computations.
// The secondary widen_t<dword_t<W>> provides 4× for ghost-map accumulation
// and unsaturated polynomial products where extra headroom is needed.

// Forward declaration for double_width<uint64_t> specialization.
namespace detail { struct uint128_t; }

template<std::unsigned_integral T>
struct double_width {};

template<> struct double_width<uint8_t>  { using type = uint16_t; };   // 2×
template<> struct double_width<uint16_t> { using type = uint32_t; };
template<> struct double_width<uint32_t> { using type = uint64_t; };
template<> struct double_width<uint64_t> { using type = detail::uint128_t; };

template<std::unsigned_integral W>
using dword_t = typename double_width<W>::type;

namespace detail {

// Software 128-bit unsigned pair. Provides constexpr arithmetic sufficient
// for carry-chain, ghost-map, and polynomial accumulation. Convertible
// to/from unsigned __int128 where available (GCC/Clang).
struct uint128_t {
    // Widths derived from the constituent word type so all limb constants
    // track sizeof(uint64_t) automatically.
    static constexpr int word_bits  = 8 * sizeof(uint64_t);  // 64
    static constexpr int half_bits  = word_bits / 2;          // 32
    static constexpr int full_bits  = 2 * word_bits;           // 128
    static constexpr uint64_t low_mask = ~uint64_t(0) >> half_bits;  // 0xFFFFFFFF

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

    // 64×64 → 128 multiplication decomposed into half-width limbs.
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
        // 128×128 → 128: low 128 bits of (a.hi*2^64 + a.lo) * (b.hi*2^64 + b.lo)
        auto p00 = mul64(a.lo, b.lo);
        auto p01 = mul64(a.lo, b.hi);
        auto p10 = mul64(a.hi, b.lo);
        // a.hi * b.hi contributes only to bits ≥128, discarded
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
    static constexpr uint128_t all_ones() { return {~0ULL, ~0ULL}; }
    static constexpr int bit_width = full_bits;

    constexpr uint128_t& operator+=(uint128_t v) { *this = *this + v; return *this; }
    constexpr uint128_t& operator-=(uint128_t v) { *this = *this - v; return *this; }
    constexpr uint128_t& operator*=(uint128_t v) { *this = *this * v; return *this; }
    constexpr uint128_t& operator<<=(int s) { *this = *this << s; return *this; }
    constexpr uint128_t& operator>>=(int s) { *this = *this >> s; return *this; }
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
    // Newton doubles precision each iteration. Starting from 1 bit,
    // after k iterations we have 2^k bits. Need 2^k ≥ bits => k = ceil(log2(bits)).
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
// Works on any unsigned integer type including software multi-word types.
template<typename U>
constexpr U div_2k_adic(U val, int k) noexcept {
    constexpr int B = 8 * sizeof(U);
    if (k <= 0) return val;
    if (k >= B) return 0;
    if ((val >> (B - 1)) == 0) return val >> k;
    return (val >> k) | (detail::all_ones_mask_v<U> << (B - k));
}

// Artin-Schreier operator ℘(x) = x^2 - x.
// Kernel: ℘(0) = ℘(1) = 0. Symmetry: ℘(x) = ℘(1 - x).
template<std::unsigned_integral W>
constexpr W artin_schreier(W x) noexcept {
    return x * x - x;
}

// ============================================================================
// 5. Bases
// ============================================================================

struct MonomialBasis {};
struct FallingFactorialBasis {};
struct TaylorBasis {};

// ============================================================================
// 6. Combinatorial Functions
// ============================================================================

namespace detail {

template<typename T>
constexpr T binom_gcd(T a, T b) {
    while (b) { T t = b; b = a % b; a = t; }
    return a;
}

} // namespace detail

template<std::unsigned_integral W, int MaxN = 67>
constexpr W binom(int n, int k) noexcept {
    if (k < 0 || k > n || n > MaxN) return W(0);
    if (k == 0 || k == n) return W(1);
    if (k > n - k) k = n - k;

    static_assert(MaxN <= 67,
        "C(67,33) exceeds uint64_t max. Reduce MaxN or ensure __int128 is available.");

#if defined(__SIZEOF_INT128__) && !defined(__STRICT_ANSI__)
    using acc_t = unsigned __int128;
    static_assert(MaxN <= 100,
        "C(n, n/2) exceeds unsigned __int128 for n > 100. Reduce MaxN.");
#else
    using acc_t = uint64_t;
#endif

    // Dual-accumulator approach: keep numerator and denominator in
    // reduced form throughout to avoid intermediate overflow.
    // For 64-bit fallback, C(60,30) ≈ 1.18e18 fits in uint64_t, but
    // the unreduced multiplication by (n-k+i) can exceed 2^64. The
    // interleaved gcd reduction prevents this.
    acc_t num = 1, den = 1;
    for (int i = 1; i <= k; ++i) {
        acc_t ni = static_cast<acc_t>(n - k + i);
        acc_t idi = static_cast<acc_t>(i);
        acc_t g;

        // Cancel gcd(ni, idi) before multiplying to keep intermediates small
        g = detail::binom_gcd(ni, idi);
        ni /= g;
        idi /= g;

        // Cancel gcd(num, idi)
        g = detail::binom_gcd(num, idi);
        num /= g;
        idi /= g;

        // Cancel gcd(den, ni)
        g = detail::binom_gcd(den, ni);
        den /= g;
        ni /= g;

        num *= ni;
        den *= idi;
    }
    return static_cast<W>(num / den);
}

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

} // namespace detail

template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_2(int n, int k) noexcept {
    if (n < 0 || k < 0 || n > MaxN) return W(0);
    if (n == 0 && k == 0) return W(1);
    if (n == 0 || k == 0) return W(0);
    if (k > n) return W(0);

    std::array<std::array<W, MaxN>, MaxN> dp{};
    dp[0][0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= i; ++j) {
            dp[i][j] = dp[i-1][j-1] + static_cast<W>(j) * dp[i-1][j];
        }
    }
    return dp[n][k];
}

// Signed Stirling numbers of the first kind s(n,k).
// The recurrence s(n,k) = s(n-1,k-1) - (n-1)*s(n-1,k) produces signed
// values. Returned as unsigned W (2-adic wrapping): e.g. s(2,1) = -1
// becomes W(-1) = 2^W - 1. Use with div_2k_adic for exact 2-adic division
// when a factorial denominator is present (see taylor_to_monomial).
template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_1(int n, int k) noexcept {
    if (n < 0 || k < 0 || n > MaxN) return W(0);
    if (n == 0 && k == 0) return W(1);
    if (n == 0 || k == 0) return W(0);
    if (k > n) return W(0);

    std::array<std::array<W, MaxN>, MaxN> dp{};
    dp[0][0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= i; ++j) {
            W term = static_cast<W>(i - 1) * dp[i-1][j];
            dp[i][j] = dp[i-1][j-1] - term;
        }
    }
    return dp[n][k];
}

template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_1_unsigned(int n, int k) noexcept {
    if (n < 0 || k < 0 || n > MaxN) return W(0);
    if (n == 0 && k == 0) return W(1);
    if (n == 0 || k == 0) return W(0);
    if (k > n) return W(0);

    std::array<std::array<W, MaxN>, MaxN> dp{};
    dp[0][0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= i; ++j) {
            dp[i][j] = dp[i-1][j-1] + static_cast<W>(i - 1) * dp[i-1][j];
        }
    }
    return dp[n][k];
}

// ============================================================================
// 7. Polynomial
// ============================================================================

template<int N, std::unsigned_integral W, typename Basis = MonomialBasis>
struct Polynomial : std::array<W, N> {
    using word_type = W;
    using basis_type = Basis;
    static constexpr int degree = N - 1;

    constexpr Polynomial() = default;
    constexpr Polynomial(std::array<W, N> coeffs) noexcept : std::array<W, N>(coeffs) {}

    constexpr W eval(W x) const noexcept {
        if constexpr (std::is_same_v<Basis, MonomialBasis>) {
            W r = 0;
            for (int i = N - 1; i >= 0; --i) {
                r = r * x + (*this)[i];
            }
            return r;
        } else if constexpr (std::is_same_v<Basis, FallingFactorialBasis>) {
            W sum = 0;
            W falling = 1;
            for (int i = 0; i < N; ++i) {
                if (i > 0) {
                    falling = falling * (x - static_cast<W>(i - 1));
                }
                sum += (*this)[i] * falling;
            }
            return sum;
        } else if constexpr (std::is_same_v<Basis, TaylorBasis>) {
            W sum = 0;
            W b = 1;
            for (int i = 0; i < N; ++i) {
                sum += (*this)[i] * b;
                W denom = static_cast<W>(i + 1);
                W num = x - static_cast<W>(i);
                if (denom % 2 == 1) {
                    b = b * num * modinv_odd(denom);
                } else {
                    int shift = v2(denom);
                    W odd_part = denom >> shift;
                    b = div_2k_adic(b * num, shift) * modinv_odd(odd_part);
                }
            }
            return sum;
        } else {
            static_assert(std::is_same_v<Basis, MonomialBasis> ||
                          std::is_same_v<Basis, FallingFactorialBasis> ||
                          std::is_same_v<Basis, TaylorBasis>,
                "Polynomial::eval: unsupported basis type");
            return W{};
        }
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
// 8. change_basis (Axiom 6)
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, FallingFactorialBasis> monomial_to_falling(const Polynomial<N, W, MonomialBasis>& p) noexcept {
    constexpr auto cache = detail::StirlingCache<N, W>{};
    Polynomial<N, W, FallingFactorialBasis> r{};
    for (int k = 0; k < N; ++k) {
        W sum = 0;
        for (int n = k; n < N; ++n) {
            sum += p[n] * cache.S2[n][k];
        }
        r[k] = sum;
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> falling_to_monomial(const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    constexpr auto cache = detail::StirlingCache<N, W>{};
    Polynomial<N, W, MonomialBasis> r{};
    for (int k = 0; k < N; ++k) {
        W sum = 0;
        for (int n = k; n < N; ++n) {
            sum += 1u * p[n] * cache.s1[n][k];
        }
        r[k] = sum;
    }
    return r;
}

// Taylor basis coefficients T_k = k! * FF_k (k! times falling factorial coef).
// T_k wraps 2^W when FF_k ≥ 2^W / k! — see taylor_to_monomial limitation.
template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, TaylorBasis> monomial_to_taylor(const Polynomial<N, W, MonomialBasis>& p) noexcept {
    constexpr auto cache = detail::StirlingCache<N, W>{};
    Polynomial<N, W, TaylorBasis> r{};
    W fact = 1;
    for (int k = 0; k < N; ++k) {
        W sum = 0;
        if (k > 0) fact *= static_cast<W>(k);
        for (int n = k; n < N; ++n) {
            sum += p[n] * cache.S2[n][k] * fact;
        }
        r[k] = sum;
    }
    return r;
}

// Compute T_j * s(j,k) / j! using double-width intermediates.
// LIMITATION: T_j = j! * FF_j may wrap 2^W when FF_j ≥ 2^W / j!. Once
// wrapped, the division by j! cannot recover the exact result. This is a
// precision window overflow analogous to Witt vector ghost-map arithmetic.
// Use small coefficients (FF_j < 2^W / j!) for exact roundtrips.
template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> taylor_to_monomial(const Polynomial<N, W, TaylorBasis>& p) noexcept {
    constexpr auto cache = detail::StirlingCache<N, W>{};
    using accum_t = detail::widen_t<dword_t<W>>;
    // Precompute factorials 0!..(N-1)!
    W facts[N];
    facts[0] = 1;
    for (int i = 1; i < N; ++i) facts[i] = facts[i-1] * static_cast<W>(i);

    Polynomial<N, W, MonomialBasis> r{};
    for (int k = 0; k < N; ++k) {
        accum_t sum{};
        for (int j = k; j < N; ++j) {
            W fact = facts[j];
            accum_t val = static_cast<accum_t>(p[j]) *
                          static_cast<accum_t>(cache.s1[j][k]);
            if (fact % 2 == 1) {
                accum_t t = val * static_cast<accum_t>(modinv_odd(fact));
                sum += t;
            } else {
                int shift = v2(fact);
                W odd = fact >> shift;
                accum_t shifted = div_2k_adic(val, shift);
                accum_t t = shifted * static_cast<accum_t>(modinv_odd(odd));
                sum += t;
            }
        }
        r[k] = static_cast<W>(sum);
    }
    return r;
}

template<typename ToBasis, int N, std::unsigned_integral W, typename FromBasis>
constexpr Polynomial<N, W, ToBasis> change_basis(const Polynomial<N, W, FromBasis>& p) noexcept {
    if constexpr (std::is_same_v<ToBasis, FromBasis>) {
        return p;
    } else if constexpr (std::is_same_v<FromBasis, MonomialBasis> && std::is_same_v<ToBasis, FallingFactorialBasis>) {
        return monomial_to_falling(p);
    } else if constexpr (std::is_same_v<FromBasis, FallingFactorialBasis> && std::is_same_v<ToBasis, MonomialBasis>) {
        return falling_to_monomial(p);
    } else if constexpr (std::is_same_v<FromBasis, MonomialBasis> && std::is_same_v<ToBasis, TaylorBasis>) {
        return monomial_to_taylor(p);
    } else if constexpr (std::is_same_v<FromBasis, TaylorBasis> && std::is_same_v<ToBasis, MonomialBasis>) {
        return taylor_to_monomial(p);
    } else {
        static_assert(
            std::is_same_v<FromBasis, FallingFactorialBasis> ||
            std::is_same_v<FromBasis, TaylorBasis>,
            "No conversion path to MonomialBasis from this basis type. "
            "Define a monomial_to_* and *_to_monomial function pair.");
        return change_basis<ToBasis>(change_basis<MonomialBasis>(p));
    }
}

// ============================================================================
// 9. synthetic_divide (Axiom 1 — The Only Carry Primitive)
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
// 10. carry_normalize
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
// 11. Formal Derivative (Axiom 2)
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, MonomialBasis> formal_derivative(const Polynomial<N, W, MonomialBasis>& p) noexcept {
    Polynomial<N-1, W, MonomialBasis> r{};
    for (int i = 1; i < N; ++i) {
        r[i-1] = static_cast<W>(i) * p[i];
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, FallingFactorialBasis> formal_derivative(const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    auto mono = change_basis<MonomialBasis>(p);
    auto dmono = formal_derivative(mono);
    return change_basis<FallingFactorialBasis>(dmono);
}

// ============================================================================
// 12. Forward Difference (Axiom 3)
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, MonomialBasis> forward_difference(const Polynomial<N, W, MonomialBasis>& p) noexcept {
    constexpr auto C = detail::PascalCache<N, W>{};
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

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, FallingFactorialBasis> forward_difference(const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    Polynomial<N-1, W, FallingFactorialBasis> r{};
    for (int i = 1; i < N; ++i) {
        r[i-1] = static_cast<W>(i) * p[i];
    }
    return r;
}

// ============================================================================
// 13. Exact Variants (Axioms 4 & 5)
// ============================================================================

template<int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N-1, W, Basis> exact_forward_difference(const Polynomial<N, W, Basis>& p) noexcept {
    return forward_difference(p);
}

template<int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N-1, W, Basis> exact_formal_derivative(const Polynomial<N, W, Basis>& p) noexcept {
    return formal_derivative(p);
}

// ============================================================================
// 14. Witt Vectors
// ============================================================================

namespace detail {

// Precomputed powers a[i]^(2^k) for the ghost map.
// Reduces ghost sum from O(N³) to O(N²).
template<int N, std::unsigned_integral W>
struct GhostPowers {
    std::array<std::array<W, N>, N> pow{};

    constexpr void init(const W* a) noexcept {
        for (int i = 0; i < N; ++i) build_row(i, a[i]);
    }

    constexpr void set(int i, W val) noexcept {
        build_row(i, val);
    }

    constexpr W get(int i, int k) const noexcept { return pow[i][k]; }

  private:
    constexpr void build_row(int i, W val) noexcept {
        pow[i][0] = val;
        for (int k = 1; k < N - i; ++k) {
            pow[i][k] = static_cast<W>(1u * pow[i][k-1] * pow[i][k-1]);
        }
    }
};

// Ghost sum using widen_t<dword_t<W>> for extra headroom.
// For W=uint8_t this is uint32_t (4×); for larger types it's the next
// level in the 2× chain (dword_t widened once more).
template<int N, std::unsigned_integral W>
constexpr auto ghost_j_widen(const GhostPowers<N, W>& pows, int j) noexcept {
    using gw_t = detail::widen_t<dword_t<W>>;
    gw_t sum{};
    for (int i = 0; i <= j; ++i) {
        sum += (gw_t(uint64_t(1)) << i) *
               static_cast<gw_t>(pows.get(i, j - i));
    }
    return sum;
}

// Backward-compatible wrapper returning dword_t<W>.
// Used by WittVector::ghost() and ghost_vector() where truncation to W
// is acceptable. Internally widens then narrows.
template<int N, std::unsigned_integral W>
constexpr dword_t<W> ghost_j_dword(const GhostPowers<N, W>& pows, int j) noexcept {
    return static_cast<dword_t<W>>(ghost_j_widen(pows, j));
}

} // namespace detail

template<int N, std::unsigned_integral W>
struct WittVector {
    std::array<W, N> a{};

    constexpr WittVector() = default;
    constexpr WittVector(std::array<W, N> c) : a(c) {}

    constexpr W& operator[](int i) noexcept { return a[i]; }
    constexpr W operator[](int i) const noexcept { return a[i]; }

    constexpr W ghost(int j) const noexcept {
        detail::GhostPowers<N, W> pows;
        pows.init(a.data());
        return static_cast<W>(detail::ghost_j_dword(pows, j));
    }

    constexpr std::array<W, N> ghost_vector() const noexcept {
        detail::GhostPowers<N, W> pows;
        pows.init(a.data());
        std::array<W, N> g{};
        for (int j = 0; j < N; ++j) {
            g[j] = static_cast<W>(detail::ghost_j_dword(pows, j));
        }
        return g;
    }

    constexpr WittVector frobenius() const noexcept {
        WittVector r;
        for (int i = 0; i < N; ++i) r.a[i] = static_cast<W>(1u * a[i] * a[i]);
        return r;
    }

    constexpr WittVector verschiebung() const noexcept {
        WittVector r;
        r.a[0] = 0;
        for (int i = 1; i < N; ++i) r.a[i] = a[i-1];
        return r;
    }

    constexpr WittVector FV() const noexcept {
        return frobenius().verschiebung();
    }

    constexpr WittVector VF() const noexcept {
        return verschiebung().frobenius();
    }

    constexpr friend WittVector operator+(const WittVector& a, const WittVector& b) noexcept {
        return witt_add(a, b);
    }

    constexpr friend WittVector operator*(const WittVector& a, const WittVector& b) noexcept {
        return witt_mul(a, b);
    }

    static constexpr WittVector zero() noexcept { return WittVector{}; }
    static constexpr WittVector one() noexcept { WittVector r; r[0] = W(1); return r; }
};

template<int N, std::unsigned_integral W>
constexpr bool check_ghost_ring(const WittVector<N,W>& a, const WittVector<N,W>& b) noexcept {
    auto ga = a.ghost_vector();
    auto gb = b.ghost_vector();
    auto gc = witt_add(a, b).ghost_vector();
    for (int j = 0; j < N; ++j) {
        if (gc[j] != ga[j] + gb[j]) return false;
    }
    return true;
}

// ============================================================================
// 15. Carry Chain (Full-Width C = (I−N)^{−1})
// ============================================================================
// The carry chain C = (I−N)^{−1}, operating on raw dword arrays.
// Completes in one pass for any input. No headroom limitation.

template<std::unsigned_integral W, typename Accum = dword_t<W>>
constexpr Accum carry_chain(W* r, const Accum* v, int N) noexcept {
    Accum carry{};
    for (int i = 0; i < N; ++i) {
        Accum sum = v[i] + carry;
        r[i] = static_cast<W>(sum);
        carry = sum >> (8 * sizeof(W));
    }
    return carry; // final carry; non-zero means overflow
}

template<std::unsigned_integral W>
constexpr W carry_chain_word(W hi, W lo) noexcept {
    dword_t<W> sum = (static_cast<dword_t<W>>(hi) << (8 * sizeof(W))) + lo;
    return static_cast<W>(sum);
}


// ============================================================================
// 16. Polynomial Multiplication
// ============================================================================
// Two-phase multiplication: unsaturated dword accumulation + carry chain.
// For results larger than the per-specialization CHUNK, uses a tiled carry
// chain that processes the unsaturated product in fixed-size chunks while
// propagating the carry between chunks. CHUNK targets ~256 bytes of accum_t
// entries, so it scales with the word size.

template<std::unsigned_integral W, typename Accum = detail::widen_t<dword_t<W>>>
constexpr void poly_mul_unsaturated(Accum* r, const W* a, int na,
                                    const W* b, int nb) noexcept {
    int nr = na + nb - 1;
    for (int i = 0; i < nr; ++i) r[i] = 0;
    for (int i = 0; i < na; ++i) {
        for (int j = 0; j < nb; ++j) {
            r[i + j] += static_cast<Accum>(a[i]) * static_cast<Accum>(b[j]);
        }
    }
}

template<std::unsigned_integral W>
constexpr void poly_mul(W* r, const W* a, int na, const W* b, int nb) noexcept {
    using accum_t = detail::widen_t<dword_t<W>>;
    constexpr int CHUNK = 256 / sizeof(accum_t);
    int nr = na + nb - 1;
    if (nr <= CHUNK) {
        std::array<accum_t, CHUNK> buf{};
        poly_mul_unsaturated(buf.data(), a, na, b, nb);
        carry_chain(r, buf.data(), nr);
    } else {
        std::array<accum_t, CHUNK> buf{};
        accum_t carry{};
        int pos = 0;

        while (pos < nr) {
            int end = pos + CHUNK;
            if (end > nr) end = nr;
            int chunk_len = end - pos;
            for (int i = 0; i < chunk_len; ++i) buf[i] = 0;

            for (int i = 0; i < na; ++i) {
                int j_start = (pos > i) ? pos - i : 0;
                int j_end = nb;
                if (i + j_end > end) j_end = end - i;
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

// ============================================================================
// 17. Taylor Shift
// ============================================================================
// P(t) → P(t + δ) using closed-form binomial transform.
// Avoids compound overflow from sequential D-application (DRAFT §6.3).

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> taylor_shift(const Polynomial<N, W, MonomialBasis>& p, W delta) noexcept {
    constexpr auto C = detail::PascalCache<N, W>{};
    Polynomial<N, W, MonomialBasis> r{};
    std::array<W, N> delta_pow{};
    delta_pow[0] = 1;
    for (int i = 1; i < N; ++i) delta_pow[i] = delta_pow[i-1] * delta;

    for (int k = 0; k < N; ++k) {
        W sum = 0;
        for (int n = k; n < N; ++n) {
            sum += p[n] * C.C[n][k] * delta_pow[n - k];
        }
        r[k] = sum;
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, FallingFactorialBasis> taylor_shift(const Polynomial<N, W, FallingFactorialBasis>& p, W delta) noexcept {
    constexpr auto C = detail::PascalCache<N, W>{};
    Polynomial<N, W, FallingFactorialBasis> r{};

    std::array<W, N> fall{};
    fall[0] = 1;
    for (int i = 1; i < N; ++i) {
        fall[i] = fall[i-1] * (delta - static_cast<W>(i - 1));
    }

    for (int i = 0; i < N; ++i) {
        W sum = 0;
        for (int k = i; k < N; ++k) {
            sum += p[k] * C.C[k][i] * fall[k - i];
        }
        r[i] = sum;
    }
    return r;
}

// ============================================================================
// 18. Indefinite Sum Σ (Inverse of Δ)
// ============================================================================
// In falling factorial basis: Σ((t)_{n-1}) = (t)_n / n
// Division by n: modinv_odd when n odd, dword right-shift when n even.
// Even n require the dividend to be divisible by 2^{v2(n)} for exactness.

template<int N, std::unsigned_integral W>
constexpr Polynomial<N+1, W, FallingFactorialBasis> indefinite_sum(const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    Polynomial<N+1, W, FallingFactorialBasis> r{};
    r[0] = 0; // constant of integration (kernel of Δ)
    for (int i = 0; i < N; ++i) {
        W denom = static_cast<W>(i + 1);
        if (denom % 2 == 1) {
            r[i+1] = p[i] * modinv_odd(denom);
        } else {
            int shift = v2(denom);
            W odd_part = denom >> shift;
            // Exact 2-adic division: divide by 2^shift (sign-extended shift)
            // then invert odd part. Requires p[i] divisible by 2^shift.
            r[i+1] = div_2k_adic(p[i], shift) * modinv_odd(odd_part);
        }
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N+1, W, MonomialBasis> indefinite_sum(const Polynomial<N, W, MonomialBasis>& p) noexcept {
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto sum_ff = indefinite_sum(ff);
    return change_basis<MonomialBasis>(sum_ff);
}

// ============================================================================
// 19–20. Witt Addition & Multiplication (ghost-map + Newton recovery)
// ============================================================================
// Algorithm: compute ghost values G_j = combine(ghost(a)_j, ghost(b)_j)
// (sum for addition, product for multiplication), then recover Witt
// components via Newton iteration.
//
// LIMITATION: Recovery succeeds only when r_j < 2^{8*sizeof(W)-j}. For larger
// components, the division by 2^j loses high bits (precision window overflow).
// This is a fundamental constraint of ghost-map arithmetic with word size W.

namespace detail {

// Newton recovery from ghost values: given G[j] (widened precision), recover
// Witt components r[j] via r_j = (G_j - S_j) / 2^j where
// S_j = Σ_{i<j} 2^i · r_i^{2^{j-i}}.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> ghost_recover(const std::array<detail::widen_t<dword_t<W>>, N>& G) noexcept {
    using gw_t = detail::widen_t<dword_t<W>>;
    WittVector<N, W> r;
    r[0] = static_cast<W>(G[0]);

    GhostPowers<N, W> r_pows;
    r_pows.set(0, r[0]);

    for (int j = 1; j < N; ++j) {
        gw_t S{};
        for (int i = 0; i < j; ++i) {
            S += (gw_t(uint64_t(1)) << i) *
                 static_cast<gw_t>(r_pows.get(i, j - i));
        }

        gw_t diff = G[j] - S;
        r[j] = static_cast<W>(diff >> j);
        r_pows.set(j, r[j]);
    }

    return r;
}

template<int N, std::unsigned_integral W, typename Combine>
constexpr WittVector<N, W> ghost_op(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b,
                                    Combine combine) noexcept {
    using gw_t = detail::widen_t<dword_t<W>>;
    GhostPowers<N, W> a_pows, b_pows;
    a_pows.init(a.a.data());
    b_pows.init(b.a.data());

    std::array<gw_t, N> G{};
    for (int j = 0; j < N; ++j) {
        G[j] = combine(ghost_j_widen(a_pows, j),
                       ghost_j_widen(b_pows, j));
    }

    return ghost_recover<N, W>(G);
}

} // namespace detail

template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_add(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b) noexcept {
    using gw_t = detail::widen_t<dword_t<W>>;
    return detail::ghost_op<N>(a, b, std::plus<gw_t>{});
}

template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_mul(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b) noexcept {
    using gw_t = detail::widen_t<dword_t<W>>;
    return detail::ghost_op<N>(a, b, std::multiplies<gw_t>{});
}

// Adams operation ψ^n: ghost_j(ψ^n(a)) = ghost_j(a)^n (componentwise power).
// A ring endomorphism on Witt vectors; ψ^p = Frobenius for prime p.
// ψ^{mn} = ψ^m ∘ ψ^n, ψ^1 = identity. Requires n ≥ 1.
// Ghost power uses widen_t<dword_t<W>> for extra headroom (4× for uint8_t,
// 2× of dword_t for others). For W=uint64, dword_t=uint128_t and no wider
// type exists, so overflow may occur for large ghost values with n ≥ 2.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> adams_operation(const WittVector<N, W>& a, int n) noexcept {
    if (n <= 1) return a;

    using gw_t = detail::widen_t<dword_t<W>>;
    detail::GhostPowers<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> G{};
    for (int j = 0; j < N; ++j) {
        gw_t g = detail::ghost_j_widen<N>(pows, j);
        gw_t gn = gw_t(uint64_t(1));
        for (int i = 0; i < n; ++i) gn = gn * g;
        G[j] = gn;
    }
    return detail::ghost_recover<N, W>(G);
}

// Teichmüller lift: τ(x) = (x, 0, ..., 0).
// Embeds the residue field into the Witt ring; τ(xy) = τ(x) × τ(y).
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> teichmueller_lift(W x) noexcept {
    WittVector<N, W> r{};
    r[0] = x;
    return r;
}

// ============================================================================
// 21. Power Series Composition and Reversion
// ============================================================================

template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<(N-1)*(M-1)+1, W, MonomialBasis>
compose(const Polynomial<N, W, MonomialBasis>& P,
        const Polynomial<M, W, MonomialBasis>& Q) noexcept {
    constexpr int R_SIZE = (N - 1) * (M - 1) + 1;

    std::array<W, R_SIZE> result{};
    std::array<W, R_SIZE> power{};
    power[0] = 1;
    int power_deg = 0;

    for (int k = 0; k < N; ++k) {
        for (int i = 0; i <= power_deg; ++i) {
            result[i] += P[k] * power[i];
        }

        if (k < N - 1) {
            int next_deg = power_deg + M - 1;
            if (next_deg >= R_SIZE) next_deg = R_SIZE - 1;

            std::array<W, R_SIZE> new_power{};
            for (int i = 0; i <= power_deg; ++i) {
                if (power[i] == 0) continue;
                for (int j = 0; j < M && i + j < R_SIZE; ++j) {
                    new_power[i + j] += power[i] * Q[j];
                }
            }
            power = new_power;
            power_deg = next_deg;
        }
    }

    return Polynomial<R_SIZE, W, MonomialBasis>(result);
}

template<int N, std::unsigned_integral W>
    requires (N >= 2)
constexpr Polynomial<N, W, MonomialBasis>
reversion(const Polynomial<N, W, MonomialBasis>& P) noexcept {
    // Precondition: P[1] is odd (invertible mod 2^{8*sizeof(W)}). P[0] is assumed 0.
    if (P[1] % 2 == 0) return Polynomial<N, W, MonomialBasis>{};
    W p1_inv = modinv_odd(P[1]);

    Polynomial<N, W, MonomialBasis> R{};
    R[0] = 0;
    R[1] = p1_inv;

    // powers[m][d] = coeff of x^d in R^m, for m = 1..N-1
    std::array<std::array<W, N>, N> powers{};
    for (int d = 0; d < N; ++d) powers[1][d] = R[d];
    for (int m = 2; m < N; ++m) {
        for (int d = 0; d < N; ++d) powers[m][d] = 0;
        for (int i = 0; i < N; ++i) {
            if (powers[m-1][i] == 0) continue;
            for (int j = 0; j < N && i + j < N; ++j) {
                powers[m][i + j] += powers[m-1][i] * R[j];
            }
        }
    }

    for (int n = 2; n < N; ++n) {
        W sum = 0;
        for (int k = 2; k <= n && k < N; ++k) {
            sum += P[k] * powers[k][n];
        }
        // P[1]*R[n] + sum = 0  =>  R[n] = -sum * p1_inv
        R[n] = (W(0) - sum) * p1_inv;

        // Incremental power update: new R[n] contributes to all higher powers.
        // Previously R[n] was 0; adding it changes every R^m where R[n] enters
        // as a factor. O(N²) per n, O(N³) total vs O(N⁴) for full recompute.
        powers[1][n] = R[n];
        // m = 2: R*R is symmetric — R[n] contributes through both factors
        if (2 * n < N) {
            powers[2][2 * n] += R[n] * R[n];
        }
        for (int i = 0; i < N - n; ++i) {
            W p = powers[1][i];
            if (p == 0) continue;
            if (i != n) powers[2][i + n] += 2 * p * R[n];
        }
        // m ≥ 3: only the second factor (R itself) changes directly
        for (int m = 3; m < N; ++m) {
            for (int i = 0; i < N - n; ++i) {
                W p = powers[m - 1][i];
                if (p == 0) continue;
                powers[m][i + n] += p * R[n];
            }
        }
    }

    return R;
}

// ============================================================================
// 22. PRNG
// ============================================================================

struct XorShift64 {
    uint64_t state;
    constexpr explicit XorShift64(uint64_t seed) : state(seed) {}
    constexpr uint64_t next() noexcept {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }
};

// ============================================================================
// 23. Verification Helpers
// ============================================================================

enum TestResult : int { PASS = 0, FAIL = 1 };

inline int report(const char* name, bool ok) noexcept {
    std::printf("%s  %s\n", name, ok ? "PASS" : "FAIL");
    return ok ? PASS : FAIL;
}

// Precision window checks — return true when recovery will be exact.
//
// Taylor basis: monomial→taylor→monomial roundtrip requires
//   T_j = j! * FF_j < 2^{8*sizeof(W)}  for all j.
// This check performs the roundtrip with widened arithmetic and compares.
template<int N, std::unsigned_integral W>
constexpr bool check_taylor_roundtrip_precision(const Polynomial<N, W, MonomialBasis>& p) noexcept {
    auto taylor = change_basis<TaylorBasis>(p);
    auto back = change_basis<MonomialBasis>(taylor);
    for (int i = 0; i < N; ++i) {
        if (p[i] != back[i]) return false;
    }
    return true;
}

// Witt vector ghost recovery: recovery r_j = (G_j - S_j) / 2^j succeeds
// exactly when r_j < 2^{8*sizeof(W)-j} for all j. This check verifies by
// computing ghost_vector, recovering, and comparing ghosts.
template<int N, std::unsigned_integral W>
constexpr bool check_witt_recovery_precision(const WittVector<N, W>& w) noexcept {
    auto gv = w.ghost_vector();
    std::array<detail::widen_t<dword_t<W>>, N> G{};
    for (int j = 0; j < N; ++j) G[j] = static_cast<detail::widen_t<dword_t<W>>>(gv[j]);
    auto recovered = detail::ghost_recover<N, W>(G);
    auto rgv = recovered.ghost_vector();
    for (int j = 0; j < N; ++j) {
        if (gv[j] != rgv[j]) return false;
    }
    return true;
}

} // namespace dyadic
