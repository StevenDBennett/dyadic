// dyadic/compose.h — Power series composition and Lagrange reversion
// Requires: dyadic/core.h

#pragma once

#include <dyadic/core.h>

namespace dyadic {

// ============================================================================
// Power Series Composition
// ============================================================================
// compose(P, Q) = P(Q(t)) truncated to degree (N-1)*(M-1).
// Complexity: O(N^2 * M^2) naive power accumulation.

template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<(N-1)*(M-1)+1, W, MonomialBasis>
compose(const Polynomial<N, W, MonomialBasis>& P,
        const Polynomial<M, W, MonomialBasis>& Q) noexcept {
    constexpr int R_SIZE = (N - 1) * (M - 1) + 1;
    static_assert(R_SIZE <= 4096, "compose: result too large, stack overflow risk");

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

// ============================================================================
// Lagrange Reversion
// ============================================================================
// Given P(t) = t + a_2 t^2 + ... with P[1] odd, find R(t) such that
// P(R(t)) = t. Uses Lagrange inversion, O(N^3).

template<int N, std::unsigned_integral W>
    requires (N >= 2)
constexpr Polynomial<N, W, MonomialBasis>
reversion(const Polynomial<N, W, MonomialBasis>& P) noexcept {
    assert((P[1] & W(1)) != W(0) && "reversion: requires P[1] odd");
    W p1_inv = modinv_odd(P[1]);

    Polynomial<N, W, MonomialBasis> R{};
    R[0] = 0;
    R[1] = p1_inv;

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
        R[n] = (W(0) - sum) * p1_inv;

        powers[1][n] = R[n];
        if (2 * n < N) {
            powers[2][2 * n] += R[n] * R[n];
        }
        for (int i = 0; i < N - n; ++i) {
            W p = powers[1][i];
            if (p == 0) continue;
            if (i != n) powers[2][i + n] += 2 * p * R[n];
        }
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

} // namespace dyadic
