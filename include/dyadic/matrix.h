// matrix.h — Matrix / Linear Algebra over Z_2 (Dyadic extension)
// Requires: dyadic/core.h

#pragma once

#include <dyadic/core.h>

#include <array>
#include <tuple>

namespace dyadic {

template<int M, int N, std::unsigned_integral W>
struct Matrix {
    std::array<std::array<W, N>, M> data{};

    constexpr std::array<W, N>& operator[](int i) { return data[i]; }
    constexpr const std::array<W, N>& operator[](int i) const { return data[i]; }

    static constexpr Matrix identity() noexcept requires (M == N) {
        Matrix r;
        for (int i = 0; i < M; ++i) r[i][i] = 1;
        return r;
    }

    static constexpr Matrix zero() noexcept {
        return {};
    }

    constexpr Matrix operator+(const Matrix& o) const noexcept {
        Matrix r;
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                r[i][j] = data[i][j] + o[i][j];
        return r;
    }

    constexpr Matrix operator-(const Matrix& o) const noexcept {
        Matrix r;
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                r[i][j] = data[i][j] - o[i][j];
        return r;
    }

    constexpr bool operator==(const Matrix& o) const noexcept {
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                if (data[i][j] != o[i][j]) return false;
        return true;
    }

    constexpr bool operator!=(const Matrix& o) const noexcept {
        return !(*this == o);
    }

    constexpr bool is_singular() const noexcept requires (M == N) {
        return determinant() % 2 == 0;
    }

    template<int P>
    constexpr Matrix<M, P, W> operator*(const Matrix<N, P, W>& o) const noexcept {
        Matrix<M, P, W> r;
        for (int i = 0; i < M; ++i)
            for (int k = 0; k < N; ++k)
                if (data[i][k] != 0)
                    for (int j = 0; j < P; ++j)
                        r[i][j] += data[i][k] * o[k][j];
        return r;
    }

    constexpr Matrix<N, M, W> transpose() const noexcept {
        Matrix<N, M, W> r;
        for (int i = 0; i < M; ++i)
            for (int j = 0; j < N; ++j)
                r[j][i] = data[i][j];
        return r;
    }

    constexpr W trace() const noexcept requires (M == N) {
        W t = 0;
        for (int i = 0; i < M; ++i) t += data[i][i];
        return t;
    }

    constexpr std::tuple<Matrix, int, W> rref() const noexcept {
        Matrix A = *this;
        W det_sign = 1;
        int rk = 0;
        auto pivot_col = [&](int col) -> bool {
            int p = rk;
            while (p < M && A[p][col] % 2 == 0) ++p;
            if (p == M) return false;
            if (p != rk) {
                std::swap(A[p], A[rk]);
                det_sign = W(0) - det_sign;
            }
            W inv = modinv_odd(A[rk][col]);
            for (int i = 0; i < M; ++i) {
                if (i == rk || A[i][col] == 0) continue;
                W factor = A[i][col] * inv;
                for (int j = col; j < N; ++j)
                    A[i][j] -= factor * A[rk][j];
            }
            ++rk;
            return true;
        };
        for (int j = 0; j < N && rk < M; ++j)
            pivot_col(j);
        return {A, rk, det_sign};
    }

    constexpr int rank() const noexcept {
        std::array<W, N> rows[M];
        for (int i = 0; i < M; ++i) rows[i] = data[i];
        int rk = 0;
        for (int col = 0; col < N && rk < M; ++col) {
            int p = rk;
            while (p < M && (rows[p][col] & 1) == 0) ++p;
            if (p == M) continue;
            if (p != rk) std::swap(rows[p], rows[rk]);
            for (int i = rk + 1; i < M; ++i) {
                if ((rows[i][col] & 1) == 0) continue;
                for (int j = col; j < N; ++j)
                    rows[i][j] -= rows[rk][j];
            }
            ++rk;
        }
        return rk;
    }

    constexpr W determinant() const noexcept requires (M == N) {
        if constexpr (M == 0) return W(1);
        if constexpr (M == 1) return data[0][0];
        if constexpr (M == 2) return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        using row_t = std::array<W, N>;
        row_t A[M];
        for (int i = 0; i < M; ++i) A[i] = data[i];
        W prev = 1;
        int parity = 0;
        for (int k = 0; k < M - 1; ++k) {
            if (A[k][k] % 2 == 0) {
                int p = k + 1;
                while (p < M && A[p][k] % 2 == 0) ++p;
                if (p == M) return 0;
                std::swap(A[p], A[k]);
                parity ^= 1;
            }
            int shift = 0;
            W odd_part = 1;
            bool need_div = (prev != 1);
            if (need_div) {
                shift = dyadic::v2(prev);
                odd_part = prev >> shift;
            }
            for (int i = k + 1; i < M; ++i) {
                for (int j = k + 1; j < M; ++j) {
                    W num = A[i][j] * A[k][k] - A[i][k] * A[k][j];
                    if (!need_div) {
                        A[i][j] = num;
                    } else {
                        A[i][j] = dyadic::div_2k_adic(num, shift) * dyadic::modinv_odd(odd_part);
                    }
                }
            }
            prev = A[k][k];
        }
        W det = A[M-1][M-1];
        return (parity == 0) ? det : (W(0) - det);
    }

    constexpr Matrix inverse() const noexcept requires (M == N) {
        std::array<W, 2 * M> rows[M]{};
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) rows[i][j] = data[i][j];
            rows[i][M + i] = 1;
        }
        for (int col = 0; col < M; ++col) {
            int p = col;
            while (p < M && rows[p][col] % 2 == 0) ++p;
            if (p == M) return {};
            if (p != col) std::swap(rows[p], rows[col]);
            W inv = modinv_odd(rows[col][col]);
            for (int i = 0; i < M; ++i) {
                if (i == col || rows[i][col] == 0) continue;
                W factor = rows[i][col] * inv;
                for (int j = col; j < 2 * M; ++j)
                    rows[i][j] -= factor * rows[col][j];
            }
        }
        Matrix r;
        for (int i = 0; i < M; ++i) {
            W inv = modinv_odd(rows[i][i]);
            for (int j = 0; j < M; ++j)
                r[i][j] = rows[i][M + j] * inv;
        }
        return r;
    }

    constexpr std::array<W, M> solve(const std::array<W, M>& b) const noexcept requires (M == N) {
        std::array<W, M + 1> rows[M]{};
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < M; ++j) rows[i][j] = data[i][j];
            rows[i][M] = b[i];
        }
        for (int col = 0; col < M; ++col) {
            int p = col;
            while (p < M && rows[p][col] % 2 == 0) ++p;
            if (p == M) return {};
            if (p != col) std::swap(rows[p], rows[col]);
            W inv = modinv_odd(rows[col][col]);
            for (int i = 0; i < M; ++i) {
                if (i == col || rows[i][col] == 0) continue;
                W factor = rows[i][col] * inv;
                for (int j = col; j <= M; ++j)
                    rows[i][j] -= factor * rows[col][j];
            }
        }
        std::array<W, M> x{};
        for (int i = 0; i < M; ++i)
            x[i] = rows[i][M] * modinv_odd(rows[i][i]);
        return x;
    }
};

} // namespace dyadic
