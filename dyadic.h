// dyadic — Six-Axiom 2-Adic Operator Calculus
// Single-header C++20 library — Sections 1-24

#pragma once

#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <functional>

#include <type_traits>
#include <concepts>
#include <cstdio>
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

    // Binary long division returning quotient.
    friend constexpr uint128_t operator/(uint128_t a, uint128_t b) {
        assert(b != uint128_t(0) && "uint128_t division by zero");
        uint128_t q = 0, r = 0;
        for (int i = full_bits - 1; i >= 0; --i) {
            r = r << 1;
            r.lo |= (a >> i).lo & 1;
            if (r >= b) {
                r = r - b;
                q.lo |= (uint64_t(1) << (i % word_bits));
                if (i >= word_bits) q.hi |= uint64_t(1) << (i - word_bits);
            }
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

// 4× width alias for the common widen_t<dword_t<W>> pattern.
// Used by ghost-map accumulation, unsaturated polynomial products.
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
// 4. Bases
// ============================================================================

struct MonomialBasis {};
struct FallingFactorialBasis {};
struct TaylorBasis {};

// ============================================================================
// 5. Combinatorial Functions
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

#if defined(__SIZEOF_INT128__) && !defined(__STRICT_ANSI__)
    using acc_t = unsigned __int128;
    static_assert(MaxN <= 100,
        "C(n, n/2) exceeds unsigned __int128 for n > 100. Reduce MaxN.");
#else
    static_assert(MaxN <= 67,
        "C(67,33) exceeds uint64_t max. Reduce MaxN.");
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

    // 1D in-place DP: S(n,k) = S(n-1,k-1) + k * S(n-1,k)
    // Iterate backwards to reuse the same row.
    std::array<W, MaxN + 1> dp{};
    dp[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j >= 1; --j) {
            dp[j] = dp[j-1] + static_cast<W>(j) * dp[j];
        }
        dp[0] = 0;
    }
    return dp[k];
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

    std::array<W, MaxN + 1> dp{};
    dp[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j >= 1; --j) {
            dp[j] = dp[j-1] - static_cast<W>(i - 1) * dp[j];
        }
        dp[0] = 0;
    }
    return dp[k];
}

template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_1_unsigned(int n, int k) noexcept {
    if (n < 0 || k < 0 || n > MaxN) return W(0);
    if (n == 0 && k == 0) return W(1);
    if (n == 0 || k == 0) return W(0);
    if (k > n) return W(0);

    std::array<W, MaxN + 1> dp{};
    dp[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j >= 1; --j) {
            dp[j] = dp[j-1] + static_cast<W>(i - 1) * dp[j];
        }
        dp[0] = 0;
    }
    return dp[k];
}

// ============================================================================
// 6. Polynomial
// ============================================================================

template<int N, std::unsigned_integral W, typename Basis = MonomialBasis>
struct Polynomial : std::array<W, N> {
    using word_type = W;
    using basis_type = Basis;
    static constexpr int max_degree = N - 1;

    // Actual degree: index of last non-zero coefficient. Returns -1 for zero poly.
    constexpr int actual_degree() const noexcept {
        for (int i = N - 1; i >= 0; --i) {
            if ((*this)[i] != W(0)) return i;
        }
        return -1;
    }

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
// 7. change_basis (Axiom 6)
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
    using accum_t = quad_width<W>;
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
// 8. synthetic_divide (Axiom 1 — The Only Carry Primitive)
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
// 10. Formal Derivative (Axiom 2)
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
// 11. Forward Difference (Axiom 3)
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
// 12. Exact Variants (Axioms 4 & 5)
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
// 13. Witt Vectors
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

// GhostPowersFull — gw_t-precision precomputed powers for ghost_map.
// Like GhostPowers but stores at gw_t (= widen_t<dword_t<W>>) precision
// Single O(N²) init.  Stores at gw_t (= quad_width<W>) precision to avoid
// truncation in log/exp/inverse recovery.
template<int N, std::unsigned_integral W>
struct GhostPowersFull {
    using gw_t = quad_width<W>;
    std::array<std::array<gw_t, N>, N> pow{};

    constexpr void init(const W* a) noexcept {
        for (int i = 0; i < N; ++i) {
            pow[i][0] = static_cast<gw_t>(a[i]);
            for (int k = 1; k < N - i; ++k) {
                pow[i][k] = pow[i][k-1] * pow[i][k-1];
            }
        }
    }

    constexpr gw_t ghost(int j) const noexcept {
        gw_t sum{};
        for (int i = 0; i <= j; ++i) {
            sum += (gw_t(uint64_t(1)) << i) * pow[i][j - i];
        }
        return sum;
    }
};

// Ghost sum using quad_width<W> (widen_t<dword_t<W>>) for extra headroom.
// For W=uint8_t this is uint32_t (4×); for larger types it's the next
// level in the 2× chain (dword_t widened once more).
template<int N, std::unsigned_integral W>
constexpr auto ghost_j_widen(const GhostPowers<N, W>& pows, int j) noexcept {
    using gw_t = quad_width<W>;
    gw_t sum{};
    for (int i = 0; i <= j; ++i) {
        sum += (gw_t(uint64_t(1)) << i) *
               static_cast<gw_t>(pows.get(i, j - i));
    }
    return sum;
}

// Wrapper from quad_width<W> → dword_t<W>.
// Ghost_j_widen accumulates at quad_width precision (4× headroom for uint8_t);
// this wrapper narrows to dword_t<W> for the WittVector::ghost() API,
// which further narrows to W on return. The dword_t intermediate prevents
// spurious overflow in ghost-compare operations.
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
// 14. Carry Chain & Hardware-Accelerated Arithmetic
// ============================================================================
// The carry chain C = (I−N)^{−1}, operating on raw dword arrays.
// Completes in one pass for any input. No headroom limitation.
//
// The adc/add_overflow/mul_overflow helpers provide constexpr-safe
// platform-optimized arithmetic, using __builtin_add_overflow and
// ADX/BMI2 intrinsics where available.

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
#include <x86intrin.h>
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
    (void)a; (void)b;
    return false;
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


// ============================================================================
// 15. Polynomial Multiplication
// ============================================================================
// Two-phase multiplication: unsaturated dword accumulation + carry chain.
// For results larger than the per-specialization CHUNK_COUNT, uses a tiled carry
// chain that processes the unsaturated product in fixed-size chunks while
// propagating the carry between chunks. CHUNK_COUNT targets ~256 bytes of accum_t
// entries, so it scales with the word size.

template<std::unsigned_integral W, typename Accum = quad_width<W>>
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
    using accum_t = quad_width<W>;
    constexpr int CHUNK_COUNT = 256 / sizeof(accum_t);  // ~256 bytes per chunk
    int nr = na + nb - 1;
    if (nr <= CHUNK_COUNT) {
        std::array<accum_t, CHUNK_COUNT> buf{};
        poly_mul_unsaturated(buf.data(), a, na, b, nb);
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
// Use this to verify poly_divmod_cw/poly_gcd_cw results.
// NOTE: This is NOT the same ring as operator* (which uses carry-chain poly_mul).
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

// Verify A = Q*B + R in the coefficient-wise ring.
// This is the CORRECT way to verify poly_divmod_cw / poly_gcd_cw results.
// (Using operator* for verification gives wrong answers — it's carry-chain.)
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

// ============================================================================
// 15b. Polynomial GCD, Resultant, and Discriminant
// ============================================================================
// Ported from carryline (BigInteger/extras.hpp). All operations are in
// monomial basis. Polynomial GCD uses pseudo-remainder PRS to avoid
// division by zero divisors (even coefficients in Z/2^W Z).
//
// RING SEMANTICS:
//   - Functions with `_cw` suffix operate coefficient-wise (standard ring).
//   - operator* uses carry-chain poly_mul (2-adic ring with carry propagation).
//   - These are DIFFERENT rings. Use poly_mul_cw() to verify _cw results.
//   - verify_divmod_cw() is provided as a correctness checker.
//
// Key constraints:
//   - poly_divmod_cw requires an odd leading coefficient (unit mod 2^W)
//   - pseudo_remainder_cw works for any coefficients (no division needed)
//   - polynomial_resultant_cw requires deg(A)+deg(B) ≤ 6 (Laplace O(n!))

namespace detail {

// Actual degree of a polynomial (index of last non-zero coeff).
// Returns -1 for zero polynomial, 0..N-1 otherwise.
// Delegates to the member function; kept as detail helper for convenience.
template<int N, std::unsigned_integral W, typename Basis>
constexpr int poly_actual_degree(const Polynomial<N, W, Basis>& p) noexcept {
    return p.actual_degree();
}

// Formal derivative in monomial basis returning std::array for internal use.
template<int N, std::unsigned_integral W>
constexpr std::array<W, N> poly_diff_arr(const std::array<W, N>& a, int deg) noexcept {
    std::array<W, N> r{};
    for (int i = 0; i < deg; ++i) {
        r[i] = static_cast<W>(static_cast<dword_t<W>>(i + 1) * static_cast<dword_t<W>>(a[i + 1]));
    }
    return r;
}

} // namespace detail

// Pseudo-remainder: prem(A, B) = lc(B)^(deg(A)-deg(B)+1) * A mod B.
// Works over any ring (no division by leading coefficient required).
template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis>
pseudo_remainder_cw(const Polynomial<N, W, MonomialBasis>& A,
                 const Polynomial<M, W, MonomialBasis>& B) noexcept {
    int da = A.actual_degree();
    int db = B.actual_degree();
    if (db < 0) return Polynomial<N, W, MonomialBasis>{}; // B = 0
    if (da < db) return A;

    // Copy A into a mutable std::array
    std::array<W, N> R{};
    for (int i = 0; i < N; ++i) R[i] = A[i];
    int nr = da;

    W lcB = B[db];
    int e = nr - db;

    for (int k = e; k >= 0; --k) {
        int pos = db + k;
        if (pos >= N) continue;
        W ck = R[pos];
        if (ck == W(0)) continue;

        // Multiply all of R by lcB
        for (int i = 0; i <= nr; ++i) R[i] = R[i] * lcB;

        // Subtract ck * x^k * B
        for (int j = 0; j <= db && k + j < N; ++j) {
            R[k + j] = R[k + j] - ck * B[j];
        }
    }

    Polynomial<N, W, MonomialBasis> result{};
    for (int i = 0; i < N; ++i) result[i] = R[i];
    return result;
}

// Polynomial division: A = Q*B + R.
// Requires B's leading coefficient to be odd (unit modulo 2^W).
// Returns (Q, R) where deg(R) < deg(B).
template<int N, int M, std::unsigned_integral W>
constexpr std::pair<Polynomial<N, W, MonomialBasis>, Polynomial<M, W, MonomialBasis>>
poly_divmod_cw(const Polynomial<N, W, MonomialBasis>& A,
             const Polynomial<M, W, MonomialBasis>& B) noexcept {
    int da = A.actual_degree();
    int db = B.actual_degree();
    assert(db >= 0 && "poly_divmod_cw: division by zero polynomial");
    assert((B[db] & W(1)) != W(0) && "poly_divmod_cw: leading coefficient must be odd");

    // Build M-sized remainder (may be smaller than N)
    auto make_remainder = [&](const std::array<W, N>& src) {
        Polynomial<M, W, MonomialBasis> r{};
        int copy_n = (M < N) ? M : N;
        for (int i = 0; i < copy_n; ++i) r[i] = src[i];
        return r;
    };

    if (da < db) {
        return {Polynomial<N, W, MonomialBasis>{}, make_remainder(A)};
    }

    W inv_lcB = modinv_odd(B[db]);

    std::array<W, N> R{};
    for (int i = 0; i < N; ++i) R[i] = A[i];

    Polynomial<N, W, MonomialBasis> Q{};

    for (int i = da - db; i >= 0; --i) {
        int pos = db + i;
        if (pos >= N) continue;
        W c = R[pos];
        if (c == W(0)) continue;
        W q = c * inv_lcB;
        Q[i] = q;
        for (int j = 0; j <= db && i + j < N; ++j) {
            R[i + j] = R[i + j] - q * B[j];
        }
    }

    return {Q, make_remainder(R)};
}

// Polynomial GCD via pseudo-remainder PRS (Euclidean algorithm).
// Works for any coefficients (even zero divisors) because it uses
// pseudo_remainder_cw instead of ordinary polynomial remainder.
// Returns the GCD normalized to have leading coefficient 1.
// NOTE: Internally converts to K=max(N,M)-sized arrays so the
// Euclidean step can swap degrees without type conflicts.
template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<(N > M ? N : M), W, MonomialBasis>
poly_gcd_cw(const Polynomial<N, W, MonomialBasis>& A_,
         const Polynomial<M, W, MonomialBasis>& B_) noexcept {
    constexpr int K = (N > M ? N : M);
    using Poly = Polynomial<K, W, MonomialBasis>;

    // Copy into K-sized internal arrays
    Poly A{}, B{};
    for (int i = 0; i < N; ++i) A[i] = A_[i];
    for (int i = 0; i < M; ++i) B[i] = B_[i];

    int da = A.actual_degree();
    int db = B.actual_degree();

    if (db < 0) {
        Poly result{};
        for (int i = 0; i <= da; ++i) result[i] = A[i];
        if (da >= 0 && result[da] != W(1) && (result[da] & W(1)) != W(0)) {
            W inv_lc = modinv_odd(result[da]);
            for (int i = 0; i <= da; ++i) result[i] = result[i] * inv_lc;
        }
        return result;
    }
    if (da < db) {
        Poly tmp = A; A = B; B = tmp;
        int t = da; da = db; db = t;
    }

    while (db >= 0) {
        auto R = pseudo_remainder_cw(
            Polynomial<K, W, MonomialBasis>(A),
            Polynomial<K, W, MonomialBasis>(B));

        A = B;
        B = Poly{};
        for (int i = 0; i < K; ++i) B[i] = R[i];

        da = A.actual_degree();
        db = B.actual_degree();
    }

    Poly result{};
    for (int i = 0; i <= da; ++i) result[i] = A[i];
    if (da >= 0 && result[da] != W(1) && (result[da] & W(1)) != W(0)) {
        W inv_lc = modinv_odd(result[da]);
        for (int i = 0; i <= da; ++i) result[i] = result[i] * inv_lc;
    }
    return result;
}

// Determinant via Laplace expansion for small matrices (dim ≤ 6).
// Fully constexpr-friendly with no heap allocation. For larger matrices,
// returns 0 and the caller should use a different method.
template<std::unsigned_integral W>
constexpr W det_laplace(const std::array<std::array<W, 6>, 6>& M, int n) noexcept {
    if (n == 0) return W(1);
    if (n == 1) return M[0][0];
    if (n == 2) return M[0][0] * M[1][1] - M[0][1] * M[1][0];

    W det = W(0);
    for (int j = 0; j < n; ++j) {
        if (M[0][j] == W(0)) continue;
        // Build submatrix (remove row 0, column j)
        std::array<std::array<W, 6>, 6> sub{};
        for (int r = 1; r < n; ++r) {
            int col = 0;
            for (int c = 0; c < n; ++c) {
                if (c == j) continue;
                sub[r - 1][col] = M[r][c];
                ++col;
            }
        }
        W cofactor = det_laplace(sub, n - 1);
        // Sign: (-1)^{0+j} = 1 if j even, -1 if j odd
        if (j % 2 == 0) {
            det += M[0][j] * cofactor;
        } else {
            det -= M[0][j] * cofactor;
        }
    }
    return det;
}

// Polynomial resultant: det(Sylvester(A, B)).
// The Sylvester matrix is (deg(A)+deg(B)) × (deg(A)+deg(B)).
// Uses Laplace expansion (O(n!)), limited to dim ≤ 6.
// For dyadic's typical small N (2..6), this covers all cases.
template<int N, int M, std::unsigned_integral W>
constexpr W polynomial_resultant_cw(const Polynomial<N, W, MonomialBasis>& A,
                                 const Polynomial<M, W, MonomialBasis>& B) noexcept {
    int m = A.actual_degree();
    int n = B.actual_degree();

    if (m < 0 || n < 0) return W(0);
    if (m == 0 && n == 0) return W(1);

    int dim = m + n;
    if (dim > 6) return W(0); // Laplace expansion limit

    // Build Sylvester matrix: dim × dim
    std::array<std::array<W, 6>, 6> Mtx{};

    // First n rows: coefficients of A (descending by index), shifted right
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j <= m; ++j) {
            Mtx[k][k + j] = A[m - j];
        }
    }
    // Last m rows: coefficients of B (descending by index), shifted right
    for (int k = 0; k < m; ++k) {
        for (int j = 0; j <= n; ++j) {
            Mtx[n + k][k + j] = B[n - j];
        }
    }

    return det_laplace<W>(Mtx, dim);
}

// Polynomial discriminant: (-1)^{d(d-1)/2} * res(P, P') / lc(P).
// Returns 0 if the leading coefficient is even (cannot divide reliably).
template<int N, std::unsigned_integral W>
constexpr W poly_discriminant_cw(const Polynomial<N, W, MonomialBasis>& P) noexcept {
    int d = P.actual_degree();
    if (d < 1) return W(0);
    if (d == 1) return W(1);

    W lc = P[d];
    if ((lc & W(1)) == W(0)) return W(0);

    // Compute P'
    auto P_prime = formal_derivative(P);

    W res = polynomial_resultant_cw(P, P_prime);
    int sign_exp = d * (d - 1) / 2;
    W sign = (sign_exp & 1) ? W(W(0) - W(1)) : W(1);
    return sign * res * modinv_odd(lc);
}

// Square-free check: gcd(P, P') of degree 0 means no repeated factors.
// Works even when the leading coefficient is even (unlike discriminant).
template<int N, std::unsigned_integral W>
constexpr bool poly_is_square_free_cw(const Polynomial<N, W, MonomialBasis>& P) noexcept {
    int d = P.actual_degree();
    if (d < 1) return false;

    auto P_prime = formal_derivative(P);
    auto g = poly_gcd_cw(P, P_prime);
    return g.actual_degree() == 0;
}

// ============================================================================
// 16. Taylor Shift
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
// 17. Indefinite Sum Σ (Inverse of Δ)
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
// 18–19. Witt Addition & Multiplication (ghost-map + Newton recovery)
// ============================================================================
// Algorithm: compute ghost values G_j = combine(ghost(a)_j, ghost(b)_j)
// (sum for addition, product for multiplication), then recover Witt
// components via Newton iteration.
//
// LIMITATION: Recovery succeeds only when r_j < 2^{8*sizeof(W)-j}. For larger
// components, the division by 2^j loses high bits (precision window overflow).
// This is a fundamental constraint of ghost-map arithmetic with word size W.

namespace detail {

// Compute ghost_j at full gw_t precision from raw Witt components.
// Unlike ghost_j_widen (which uses W-precision GhostPowers), this computes
// each a[i]^(2^{j-i}) at gw_t precision, avoiding truncation. Required when
// Witt components are large (e.g., after inverse/log/exp recovery).
template<int N, std::unsigned_integral W>
constexpr auto ghost_j_full(const W* a, int j) noexcept {
    using gw_t = quad_width<W>;
    gw_t sum{};
    for (int i = 0; i <= j; ++i) {
        gw_t pow = static_cast<gw_t>(a[i]);
        for (int k = 0; k < j - i; ++k) {
            pow = pow * pow;
        }
        sum += (gw_t(uint64_t(1)) << i) * pow;
    }
    return sum;
}

// Newton recovery from ghost values: given G[j] (widened precision), recover
// Witt components r[j] via r_j = (G_j - S_j) / 2^j where
// S_j = Σ_{i<j} 2^i · r_i^{2^{j-i}}.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> ghost_recover(const std::array<quad_width<W>, N>& G) noexcept {
    using gw_t = quad_width<W>;
    WittVector<N, W> r;
    r[0] = static_cast<W>(G[0]);

    // Precompute r[i]^(2^k) at gw_t precision to avoid truncation
    std::array<std::array<gw_t, N>, N> pows{};
    pows[0][0] = static_cast<gw_t>(r[0]);
    for (int k = 1; k < N; ++k) {
        pows[0][k] = pows[0][k-1] * pows[0][k-1];
    }

    for (int j = 1; j < N; ++j) {
        gw_t S{};
        for (int i = 0; i < j; ++i) {
            S += (gw_t(uint64_t(1)) << i) * pows[i][j - i];
        }

        gw_t diff = G[j] - S;
        r[j] = static_cast<W>(diff >> j);

        pows[j][0] = static_cast<gw_t>(r[j]);
        for (int k = 1; k < N - j; ++k) {
            pows[j][k] = pows[j][k-1] * pows[j][k-1];
        }
    }

    return r;
}

template<int N, std::unsigned_integral W, typename Combine>
constexpr WittVector<N, W> ghost_op(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b,
                                    Combine combine) noexcept {
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> a_pows, b_pows;
    a_pows.init(a.a.data());
    b_pows.init(b.a.data());

    std::array<gw_t, N> G{};
    for (int j = 0; j < N; ++j) {
        G[j] = combine(a_pows.ghost(j), b_pows.ghost(j));
    }

    return ghost_recover<N, W>(G);
}

} // namespace detail

template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_add(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b) noexcept {
    using gw_t = quad_width<W>;
    return detail::ghost_op<N>(a, b, std::plus<gw_t>{});
}

template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_mul(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b) noexcept {
    using gw_t = quad_width<W>;
    return detail::ghost_op<N>(a, b, std::multiplies<gw_t>{});
}

// Adams operation ψ^n: ghost_j(ψ^n(a)) = ghost_j(a)^n (componentwise power).
// A ring endomorphism on Witt vectors; ψ^p = Frobenius for prime p.
// ψ^{mn} = ψ^m ∘ ψ^n, ψ^1 = identity. Requires n ≥ 1.
// Ghost power uses quad_width<W> for extra headroom (4× for uint8_t,
// 2× of dword_t for others). For W=uint64, dword_t=uint128_t and no wider
// type exists, so overflow may occur for large ghost values with n ≥ 2.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> adams_operation(const WittVector<N, W>& a, int n) noexcept {
    if (n <= 1) return a;

    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> G{};
    for (int j = 0; j < N; ++j) {
        gw_t g = pows.ghost(j);
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
// 20. Power Series Composition and Reversion
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
    assert((P[1] & W(1)) != W(0) && "reversion: requires P[1] odd");
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
// 21. PRNG
// ============================================================================

struct XorShift64 {
    uint64_t state;
    constexpr explicit XorShift64(uint64_t seed) : state(seed) {
        // XorShift64 has a fixed point at 0; require non-zero seed.
        // All seeds used in the test suite are non-zero.
    }
    constexpr uint64_t next() noexcept {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }
    static constexpr bool is_valid_seed(uint64_t s) { return s != 0; }
};

// ============================================================================
// 22. Verification Helpers
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
    std::array<quad_width<W>, N> G{};
    for (int j = 0; j < N; ++j) G[j] = static_cast<quad_width<W>>(gv[j]);
    auto recovered = detail::ghost_recover<N, W>(G);
    auto rgv = recovered.ghost_vector();
    for (int j = 0; j < N; ++j) {
        if (gv[j] != rgv[j]) return false;
    }
    return true;
}

// ============================================================================
// 23. Witt Logarithm, Exponential, and Inverse
// ============================================================================
// The ghost map ghost: W(R) → R^N is a ring isomorphism onto its image.
// Componentwise log/exp/inverse of ghost values followed by ghost recovery
// gives the corresponding Witt vector operations.
//   - witt_log(a): requires a[0] odd (unit in the Witt ring)
//   - witt_exp(a): requires a[0] ≡ 0 (mod 4) for convergence in ℤ₂
//   - witt_inverse(a): requires a[0] odd
//
// Ghost map operations run at gw_t (= quad_width<W>) precision, which
// for uint64_t is uint128_t. The 2-adic primitives (v2, modinv_odd, div_2k_adic)
// are overloaded for uint128_t in the detail namespace.

namespace detail {

// 2-adic primitives for uint128_t (not std::unsigned_integral)
// Named with _128 suffix to avoid hiding dyadic::v2, modinv_odd, div_2k_adic.
constexpr int v2_128(uint128_t x) noexcept {
    if (x.lo != 0) return std::countr_zero(x.lo);
    if (x.hi != 0) return 64 + std::countr_zero(x.hi);
    return 128;
}

constexpr uint128_t modinv_odd_128(uint128_t a) noexcept {
    uint128_t x = 1;
    x = x * (uint128_t(2) - a * x);
    x = x * (uint128_t(2) - a * x);
    x = x * (uint128_t(2) - a * x);
    x = x * (uint128_t(2) - a * x);
    x = x * (uint128_t(2) - a * x);
    x = x * (uint128_t(2) - a * x);
    x = x * (uint128_t(2) - a * x);
    return x;
}

constexpr uint128_t div_2k_adic_128(uint128_t val, int k) noexcept {
    if (k <= 0) return val;
    bool negative = (val.hi >> 63) != 0;
    val = val >> k;
    if (negative) {
        uint128_t mask = ~(uint128_t::all_ones() >> k);
        val = val | mask;
    }
    return val;
}

// Divide a ℤ₂ value (generic type T) by integer n using 2-adic division.
template<typename T>
constexpr T div_by_int_generic(T val, int n) noexcept {
    T n_t = static_cast<T>(static_cast<uint64_t>(n));
    int shift;
    T odd_part;
    if constexpr (std::is_same_v<T, uint128_t>) {
        shift = v2_128(n_t);
        odd_part = n_t >> shift;
        return div_2k_adic_128(val, shift) * modinv_odd_128(odd_part);
    } else {
        shift = dyadic::v2(n_t);
        odd_part = n_t >> shift;
        return dyadic::div_2k_adic(val, shift) * dyadic::modinv_odd(odd_part);
    }
}

// 2-adic valuation for the generic types used in p-adic series.
template<typename T>
constexpr int p_adic_val(T x) noexcept {
    if constexpr (std::is_same_v<T, uint128_t>) {
        return v2_128(x);
    } else {
        return dyadic::v2(x);
    }
}

// ℤ₂ logarithm: log(1+y) = Σ_{n=1}^{∞} (-1)^{n+1} y^n / n
// Converges when v₂(y) ≥ 1 (y even). Terms have valuation n·v₂(y) − v₂(n),
// vanishing when this reaches the bit width of T.
// For v₂(y) = 1, needs ~bits + log₂(bits) terms (135 for 128-bit).
// Uses a 2× budget to cover the v₂(y) = 1 case.
// NOTE: This computes log(1+y), NOT log(y). For log(y), pass y-1.
template<typename T>
constexpr T p_adic_log_1plus_impl(T y) noexcept {
    if (y == T(0)) return T(0);
    int bits = 8 * int(sizeof(T));
    int max_terms = 2 * bits;  // generous budget for v₂=1 inputs
    T sum = y;
    T term = y;
    for (int n = 2; n <= max_terms; ++n) {
        term = term * y;
        term = term * T(static_cast<uint64_t>(n - 1));
        term = div_by_int_generic(term, n);
        if (term == T(0)) break;
        sum = (n & 1) ? (sum + term) : (sum - term);
    }
    return sum;
}

// ℤ₂ exponential: exp(x) = Σ_{n=0}^{∞} x^n / n!
// Converges when v₂(x) ≥ 2 (x ≡ 0 mod 4). Terms have valuation
// n·(v₂(x)−1) + s₂(n), vanishing when this reaches the bit width of T.
// For v₂(x) = 1, the series never converges at finite precision
// (s₂(n) ≤ log₂(n)+1). Uses a 2× budget and early-exits on vanishing terms.
template<typename T>
constexpr T p_adic_exp_impl(T x) noexcept {
    int bits = 8 * int(sizeof(T));
    int v = p_adic_val(x);
    int max_terms = (v < 2) ? (2 * bits) : ((bits + v - 3) / (v - 1) + 1);
    if (max_terms > 2 * bits) max_terms = 2 * bits;
    T sum = T(1);
    T term = T(1);
    for (int n = 1; n <= max_terms; ++n) {
        term = term * x;
        term = div_by_int_generic(term, n);
        if (term == T(0)) break;
        sum = sum + term;
    }
    return sum;
}

} // namespace detail

// Witt vector logarithm: applies ℤ₂ log componentwise to ghost values.
// Defined when a[0] is odd (unit in the Witt ring).
// Returns b such that witt_exp(b) = a (when both are defined).
//
// p_adic_log_1plus_impl(y) computes log(1+y). We want log(ghost_j(a)),
// we pass ghost_j(a) - 1 so that log(1 + (ghost_j(a)-1)) = log(ghost_j(a)).
// Convergence requires ghost_j(a) ≡ 1 (mod 2), which holds when a[0] is odd
// (a[0]^{2^j} ≡ 1 mod 2 for odd a[0]).
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_log(const WittVector<N, W>& a) noexcept {
    assert((a[0] & W(1)) != W(0) && "witt_log: requires a[0] odd (unit in Witt ring)");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> H{};
    for (int j = 0; j < N; ++j) {
        H[j] = detail::p_adic_log_1plus_impl(pows.ghost(j) - gw_t(1));
    }
    return detail::ghost_recover<N, W>(H);
}

// Witt vector exponential: applies ℤ₂ exp componentwise to ghost values.
// Requires a[0] ≡ 0 (mod 4) for convergence in the 2-adic series.
// Returns b such that witt_log(b) = a (when both are defined).
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_exp(const WittVector<N, W>& a) noexcept {
    assert((a[0] & W(3)) == W(0) && "witt_exp: requires a[0] ≡ 0 (mod 4), v₂(a[0]) ≥ 2");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> H{};
    for (int j = 0; j < N; ++j) {
        H[j] = detail::p_adic_exp_impl(pows.ghost(j));
    }
    return detail::ghost_recover<N, W>(H);
}

// Witt vector inverse: a^{-1} in the Witt ring.
// Defined when a[0] is odd (unit). Uses componentwise ghost inversion.
// Satisfies a * witt_inverse(a) = τ(1) = {1, 0, ..., 0}.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_inverse(const WittVector<N, W>& a) noexcept {
    assert((a[0] & W(1)) != W(0) && "witt_inverse: requires a[0] odd (unit in Witt ring)");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> H{};
    for (int j = 0; j < N; ++j) {
        gw_t Gj = pows.ghost(j);
        if constexpr (std::is_same_v<gw_t, detail::uint128_t>) {
            H[j] = detail::modinv_odd_128(Gj);
        } else {
            H[j] = dyadic::modinv_odd(Gj);
        }
    }
    return detail::ghost_recover<N, W>(H);
}

// ============================================================================
// 24. Verification Helpers for Witt Log/Exp/Inverse
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr bool check_witt_inverse(const WittVector<N, W>& a) noexcept {
    auto inv = witt_inverse(a);
    auto prod = a * inv;
    for (int i = 0; i < N; ++i) {
        if (i == 0 && prod[i] != 1) return false;
        if (i > 0 && prod[i] != 0) return false;
    }
    return true;
}

// Check exp∘log = identity on the multiplicative domain (odd a[0]).
// For the additive roundtrip (log∘exp = id on v₂(a₀) ≥ 2), see PROOF_WITT_LOG_EXP.
template<int N, std::unsigned_integral W>
constexpr bool check_witt_exp_log_roundtrip(const WittVector<N, W>& a) noexcept {
    auto b = witt_exp(witt_log(a));
    for (int i = 0; i < N; ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

} // namespace dyadic
