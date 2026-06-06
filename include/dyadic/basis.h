// dyadic/basis.h — Basis Conversion (Axiom 6) + FF/Taylor Polynomial Operations
// Requires: dyadic/core.h
//
// Provides change_basis between Monomial, FallingFactorial, and Taylor bases,
// plus FF/Taylor overloads for eval, formal_derivative, forward_difference.

#pragma once

#include <dyadic/core.h>

namespace dyadic {

// ============================================================================
// Basis Conversion (Axiom 6)
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, FallingFactorialBasis> monomial_to_falling(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    constexpr auto& cache = detail::STIRLING_CACHE<N, W>;
    using dw_t = dword_t<W>;
    Polynomial<N, W, FallingFactorialBasis> r{};
    for (int k = 0; k < N; ++k) {
        dw_t sum = 0;
        for (int n = k; n < N; ++n) {
            sum += static_cast<dw_t>(p[n]) * static_cast<dw_t>(cache.S2[n][k]);
        }
        r[k] = static_cast<W>(sum);
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> falling_to_monomial(
    const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    constexpr auto& cache = detail::STIRLING_CACHE<N, W>;
    using dw_t = dword_t<W>;
    Polynomial<N, W, MonomialBasis> r{};
    for (int k = 0; k < N; ++k) {
        dw_t sum = 0;
        for (int n = k; n < N; ++n) {
            sum += static_cast<dw_t>(p[n]) * static_cast<dw_t>(cache.s1[n][k]);
        }
        r[k] = static_cast<W>(sum);
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, TaylorBasis> monomial_to_taylor(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    constexpr auto& cache = detail::STIRLING_CACHE<N, W>;
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

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> taylor_to_monomial(
    const Polynomial<N, W, TaylorBasis>& p) noexcept {
    constexpr auto& cache = detail::STIRLING_CACHE<N, W>;
    using accum_t = quad_width<W>;
    std::array<W, N> facts{};
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
constexpr Polynomial<N, W, ToBasis> change_basis(
    const Polynomial<N, W, FromBasis>& p) noexcept {
    if constexpr (std::is_same_v<ToBasis, FromBasis>) {
        return p;
    } else if constexpr (std::is_same_v<FromBasis, MonomialBasis> &&
                         std::is_same_v<ToBasis, FallingFactorialBasis>) {
        return monomial_to_falling(p);
    } else if constexpr (std::is_same_v<FromBasis, FallingFactorialBasis> &&
                         std::is_same_v<ToBasis, MonomialBasis>) {
        return falling_to_monomial(p);
    } else if constexpr (std::is_same_v<FromBasis, MonomialBasis> &&
                         std::is_same_v<ToBasis, TaylorBasis>) {
        return monomial_to_taylor(p);
    } else if constexpr (std::is_same_v<FromBasis, TaylorBasis> &&
                         std::is_same_v<ToBasis, MonomialBasis>) {
        return taylor_to_monomial(p);
    } else {
        return change_basis<ToBasis>(change_basis<MonomialBasis>(p));
    }
}

// ============================================================================
// Polynomial eval — FF and Taylor bases
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr W eval(const Polynomial<N, W, FallingFactorialBasis>& p, W x) noexcept {
    W sum = 0;
    W falling = 1;
    for (int i = 0; i < N; ++i) {
        if (i > 0) falling = falling * (x - static_cast<W>(i - 1));
        sum += p[i] * falling;
    }
    return sum;
}

template<int N, std::unsigned_integral W>
constexpr W eval(const Polynomial<N, W, TaylorBasis>& p, W x) noexcept {
    W sum = 0;
    W b = 1;
    for (int i = 0; i < N; ++i) {
        sum += p[i] * b;
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
}

// ============================================================================
// Formal derivative — FF and Taylor bases
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, FallingFactorialBasis> formal_derivative(
    const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    auto mono = change_basis<MonomialBasis>(p);
    auto dmono = formal_derivative(mono);
    return change_basis<FallingFactorialBasis>(dmono);
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, TaylorBasis> formal_derivative(
    const Polynomial<N, W, TaylorBasis>& p) noexcept {
    auto mono = change_basis<MonomialBasis>(p);
    auto dmono = formal_derivative(mono);
    return change_basis<TaylorBasis>(dmono);
}

// ============================================================================
// Forward difference — FF and Taylor bases
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, FallingFactorialBasis> forward_difference(
    const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    Polynomial<N-1, W, FallingFactorialBasis> r{};
    for (int i = 1; i < N; ++i) {
        r[i-1] = static_cast<W>(i) * p[i];
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N-1, W, TaylorBasis> forward_difference(
    const Polynomial<N, W, TaylorBasis>& p) noexcept {
    auto mono = change_basis<MonomialBasis>(p);
    auto dmono = forward_difference(mono);
    return change_basis<TaylorBasis>(dmono);
}

} // namespace dyadic
