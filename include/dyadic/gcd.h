// dyadic/gcd.h — Polynomial GCD, Resultant, Discriminant (coefficient-wise ring)
// Requires: dyadic/core.h
//
// All operations are in monomial basis with `_cw` suffix (coefficient-wise ring).
// For the carry-chain ring, use operator* instead.

#pragma once

#include <dyadic/core.h>

namespace dyadic {

namespace detail {

template<int N, std::unsigned_integral W, typename Basis>
constexpr int poly_actual_degree(const Polynomial<N, W, Basis>& p) noexcept {
    return p.actual_degree();
}

template<int N, std::unsigned_integral W>
constexpr std::array<W, N> poly_diff_arr(const std::array<W, N>& a, int deg) noexcept {
    std::array<W, N> r{};
    for (int i = 0; i < deg; ++i) {
        r[i] = static_cast<W>(static_cast<dword_t<W>>(i + 1) *
                              static_cast<dword_t<W>>(a[i + 1]));
    }
    return r;
}

} // namespace detail

// ============================================================================
// Pseudo-remainder
// ============================================================================

template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<N, W, MonomialBasis>
pseudo_remainder_cw(const Polynomial<N, W, MonomialBasis>& A,
                    const Polynomial<M, W, MonomialBasis>& B) noexcept {
    int da = A.actual_degree();
    int db = B.actual_degree();
    if (db < 0) return Polynomial<N, W, MonomialBasis>{};
    if (da < db) return A;

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
        for (int i = 0; i <= nr; ++i) R[i] = R[i] * lcB;
        for (int j = 0; j <= db && k + j < N; ++j) {
            R[k + j] = R[k + j] - ck * B[j];
        }
    }

    Polynomial<N, W, MonomialBasis> result{};
    for (int i = 0; i < N; ++i) result[i] = R[i];
    return result;
}

// ============================================================================
// Polynomial division with remainder
// ============================================================================

template<int N, int M, std::unsigned_integral W>
constexpr std::pair<Polynomial<N, W, MonomialBasis>, Polynomial<M, W, MonomialBasis>>
poly_divmod_cw(const Polynomial<N, W, MonomialBasis>& A,
               const Polynomial<M, W, MonomialBasis>& B) noexcept {
    int da = A.actual_degree();
    int db = B.actual_degree();
    assert(db >= 0 && "poly_divmod_cw: division by zero polynomial");
    assert((B[db] & W(1)) != W(0) && "poly_divmod_cw: leading coefficient must be odd");

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

// ============================================================================
// Polynomial GCD (coefficient-wise ring)
// ============================================================================

template<int N, int M, std::unsigned_integral W>
constexpr Polynomial<(N > M ? N : M), W, MonomialBasis>
poly_gcd_cw(const Polynomial<N, W, MonomialBasis>& A_,
            const Polynomial<M, W, MonomialBasis>& B_) noexcept {
    constexpr int K = (N > M ? N : M);
    using Poly = Polynomial<K, W, MonomialBasis>;

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

// ============================================================================
// Determinant (Laplace expansion, dim <= 6)
// ============================================================================

namespace detail {

template<std::unsigned_integral W>
constexpr W det_laplace(const std::array<std::array<W, 6>, 6>& M, int n) noexcept {
    if (n == 0) return W(1);
    if (n == 1) return M[0][0];
    if (n == 2) return M[0][0] * M[1][1] - M[0][1] * M[1][0];

    W det = W(0);
    for (int j = 0; j < n; ++j) {
        if (M[0][j] == W(0)) continue;
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
        if (j % 2 == 0) {
            det += M[0][j] * cofactor;
        } else {
            det -= M[0][j] * cofactor;
        }
    }
    return det;
}

} // namespace detail

// ============================================================================
// Polynomial resultant
// ============================================================================

template<int N, int M, std::unsigned_integral W>
constexpr W polynomial_resultant_cw(const Polynomial<N, W, MonomialBasis>& A,
                                    const Polynomial<M, W, MonomialBasis>& B) noexcept {
    int m = A.actual_degree();
    int n = B.actual_degree();

    if (m < 0 || n < 0) return W(0);
    if (m == 0 && n == 0) return W(1);

    int dim = m + n;
    if (dim > 6) return W(0);

    std::array<std::array<W, 6>, 6> Mtx{};

    for (int k = 0; k < n; ++k) {
        for (int j = 0; j <= m; ++j) {
            Mtx[k][k + j] = A[m - j];
        }
    }
    for (int k = 0; k < m; ++k) {
        for (int j = 0; j <= n; ++j) {
            Mtx[n + k][k + j] = B[n - j];
        }
    }

    return detail::det_laplace<W>(Mtx, dim);
}

// ============================================================================
// Polynomial discriminant
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr W poly_discriminant_cw(const Polynomial<N, W, MonomialBasis>& P) noexcept {
    int d = P.actual_degree();
    if (d < 1) return W(0);
    if (d == 1) return W(1);

    W lc = P[d];
    if ((lc & W(1)) == W(0)) return W(0);

    auto P_prime = formal_derivative(P);
    W res = polynomial_resultant_cw(P, P_prime);
    int sign_exp = d * (d - 1) / 2;
    W sign = (sign_exp & 1) ? W(W(0) - W(1)) : W(1);
    return sign * res * modinv_odd(lc);
}

// ============================================================================
// Square-free check
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr bool poly_is_square_free_cw(const Polynomial<N, W, MonomialBasis>& P) noexcept {
    int d = P.actual_degree();
    if (d < 1) return false;

    auto P_prime = formal_derivative(P);
    auto g = poly_gcd_cw(P, P_prime);
    return g.actual_degree() == 0;
}

} // namespace dyadic
