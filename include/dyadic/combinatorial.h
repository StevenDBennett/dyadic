// dyadic/combinatorial.h — Binomial coefficients, Stirling numbers
// Free-standing combinatorial functions. No dyadic-core dependencies
// beyond the C++20 standard library.
//
// Requires: <array>, <concepts>, <cstdint>

#pragma once

#include <array>
#include <concepts>
#include <cstdint>

namespace dyadic {

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
    acc_t num = 1, den = 1;
    for (int i = 1; i <= k; ++i) {
        acc_t ni = static_cast<acc_t>(n - k + i);
        acc_t idi = static_cast<acc_t>(i);
        acc_t g = detail::binom_gcd(ni, idi);
        ni /= g; idi /= g;
        g = detail::binom_gcd(num, idi);
        num /= g; idi /= g;
        g = detail::binom_gcd(den, ni);
        den /= g; ni /= g;
        num *= ni; den *= idi;
    }
    return static_cast<W>(num / den);
}

template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_2(int n, int k) noexcept {
    if (n < 0 || k < 0 || n > MaxN) return W(0);
    if (n == 0 && k == 0) return W(1);
    if (n == 0 || k == 0) return W(0);
    if (k > n) return W(0);
    std::array<W, MaxN + 1> dp{};
    dp[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j >= 1; --j)
            dp[j] = dp[j - 1] + static_cast<W>(j) * dp[j];
        dp[0] = 0;
    }
    return dp[k];
}

template<std::unsigned_integral W, int MaxN = 64>
constexpr W stirling_1(int n, int k) noexcept {
    if (n < 0 || k < 0 || n > MaxN) return W(0);
    if (n == 0 && k == 0) return W(1);
    if (n == 0 || k == 0) return W(0);
    if (k > n) return W(0);
    std::array<W, MaxN + 1> dp{};
    dp[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = i; j >= 1; --j)
            dp[j] = dp[j - 1] - static_cast<W>(i - 1) * dp[j];
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
        for (int j = i; j >= 1; --j)
            dp[j] = dp[j - 1] + static_cast<W>(i - 1) * dp[j];
        dp[0] = 0;
    }
    return dp[k];
}

} // namespace dyadic
