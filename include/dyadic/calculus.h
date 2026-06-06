// dyadic/calculus.h — Taylor Shift and Indefinite Sum
// Requires: dyadic/basis.h

#pragma once

#include <dyadic/basis.h>

namespace dyadic {

// ============================================================================
// Taylor Shift — P(t) -> P(t + delta) using closed-form binomial transform
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis> taylor_shift(
    const Polynomial<N, W, MonomialBasis>& p, W delta) noexcept {
    constexpr auto& C = detail::PASCAL_CACHE<N, W>;
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
constexpr Polynomial<N, W, FallingFactorialBasis> taylor_shift(
    const Polynomial<N, W, FallingFactorialBasis>& p, W delta) noexcept {
    constexpr auto& C = detail::PASCAL_CACHE<N, W>;
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

template<int N, std::unsigned_integral W>
constexpr Polynomial<N, W, TaylorBasis> taylor_shift(
    const Polynomial<N, W, TaylorBasis>& p, W delta) noexcept {
    auto mono = change_basis<MonomialBasis>(p);
    auto shifted = taylor_shift(mono, delta);
    return change_basis<TaylorBasis>(shifted);
}

// ============================================================================
// Indefinite Sum Σ (Inverse of Delta)
// ============================================================================
// In falling factorial basis: Σ((t)_{n-1}) = (t)_n / n

template<int N, std::unsigned_integral W>
constexpr Polynomial<N+1, W, FallingFactorialBasis> indefinite_sum(
    const Polynomial<N, W, FallingFactorialBasis>& p) noexcept {
    Polynomial<N+1, W, FallingFactorialBasis> r{};
    r[0] = 0;
    for (int i = 0; i < N; ++i) {
        W denom = static_cast<W>(i + 1);
        if (denom % 2 == 1) {
            r[i+1] = p[i] * modinv_odd(denom);
        } else {
            int shift = v2(denom);
            W odd_part = denom >> shift;
            r[i+1] = div_2k_adic(p[i], shift) * modinv_odd(odd_part);
        }
    }
    return r;
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N+1, W, MonomialBasis> indefinite_sum(
    const Polynomial<N, W, MonomialBasis>& p) noexcept {
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto sum_ff = indefinite_sum(ff);
    return change_basis<MonomialBasis>(sum_ff);
}

template<int N, std::unsigned_integral W>
constexpr Polynomial<N+1, W, TaylorBasis> indefinite_sum(
    const Polynomial<N, W, TaylorBasis>& p) noexcept {
    auto mono = change_basis<MonomialBasis>(p);
    auto sum = indefinite_sum(mono);
    return change_basis<TaylorBasis>(sum);
}

} // namespace dyadic
