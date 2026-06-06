// pade.h — Pade Approximants in Z_2[[t]] (Dyadic extension)
// Requires: dyadic/core.h

#pragma once

#include <dyadic/core.h>

#include <array>
#include <utility>

namespace dyadic {

template<int M, int N, std::unsigned_integral W>
constexpr std::pair<
    Polynomial<M+1, W, MonomialBasis>,
    Polynomial<N+1, W, MonomialBasis>
>
pade_approximant(const Polynomial<M+N+1, W, MonomialBasis>& series) noexcept {
    using P = Polynomial<M+1, W, MonomialBasis>;
    using Q = Polynomial<N+1, W, MonomialBasis>;

    if constexpr (N == 0) {
        std::array<W, M+1> p{};
        for (int i = 0; i <= M; ++i) p[i] = series[i];
        return {p, std::array<W, 1>{W(1)}};
    }

    std::array<std::array<W, N + 1>, N> mat{};
    for (int i = 0; i < N; ++i) {
        mat[i][N] = W(0) - series[M + 1 + i];
        for (int j = 0; j < N; ++j) {
            int idx = M + i - j;
            mat[i][j] = (idx >= 0) ? series[idx] : W(0);
        }
    }

    for (int col = 0; col < N; ++col) {
        int pivot = col;
        while (pivot < N && mat[pivot][col] % 2 == 0) ++pivot;
        if (pivot == N) {
            return {P{}, Q{}};
        }
        if (pivot != col) {
            for (int j = col; j <= N; ++j) std::swap(mat[col][j], mat[pivot][j]);
        }
        W inv_pivot = modinv_odd(mat[col][col]);
        for (int i = col + 1; i < N; ++i) {
            if (mat[i][col] == 0) continue;
            W factor = mat[i][col] * inv_pivot;
            for (int j = col; j <= N; ++j) {
                mat[i][j] -= factor * mat[col][j];
            }
        }
    }

    std::array<W, N+1> q{};
    q[0] = W(1);
    for (int i = N - 1; i >= 0; --i) {
        W sum = 0;
        for (int j = i + 1; j < N; ++j) {
            sum += mat[i][j] * q[j + 1];
        }
        q[i + 1] = (mat[i][N] - sum) * modinv_odd(mat[i][i]);
    }

    std::array<W, M+1> p{};
    for (int k = 0; k <= M; ++k) {
        W sum = series[k];
        for (int j = 1; j <= N && j <= k; ++j) {
            sum += q[j] * series[k - j];
        }
        p[k] = sum;
    }

    return {p, q};
}

} // namespace dyadic
