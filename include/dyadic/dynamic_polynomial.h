// dynamic_polynomial.h — Runtime-degree polynomial (Dyadic extension)
// Requires: dyadic/core.h

#pragma once

#include <dyadic/core.h>

#include <vector>

namespace dyadic {

namespace detail {

template<std::unsigned_integral W>
struct DynamicStirlingCache {
    std::vector<std::vector<W>> S2;
    std::vector<std::vector<W>> s1;

    DynamicStirlingCache(int N) : S2(N, std::vector<W>(N)), s1(N, std::vector<W>(N)) {
        if (N > 0) {
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

} // namespace detail

template<std::unsigned_integral W, typename Basis = MonomialBasis>
struct DynamicPolynomial {
    using word_type = W;
    using basis_type = Basis;

    std::vector<W> coeff;

    DynamicPolynomial() = default;
    explicit DynamicPolynomial(int deg) : coeff(deg + 1) {}
    DynamicPolynomial(std::vector<W> c) : coeff(std::move(c)) {}
    DynamicPolynomial(std::initializer_list<W> il) : coeff(il) {}

    template<int N>
    DynamicPolynomial(const Polynomial<N, W, Basis>& p) : coeff(p.begin(), p.end()) {}

    W& operator[](int i) { return coeff[i]; }
    W operator[](int i) const { return coeff[i]; }
    int size() const { return static_cast<int>(coeff.size()); }
    int degree() const {
        if (coeff.empty()) return -1;
        int d = static_cast<int>(coeff.size()) - 1;
        while (d > 0 && coeff[d] == 0) --d;
        return d;
    }

    DynamicPolynomial& operator+=(const DynamicPolynomial& o) {
        if (o.coeff.size() > coeff.size()) coeff.resize(o.coeff.size());
        for (int i = 0; i < static_cast<int>(o.coeff.size()); ++i)
            coeff[i] += o.coeff[i];
        return *this;
    }

    DynamicPolynomial& operator-=(const DynamicPolynomial& o) {
        if (o.coeff.size() > coeff.size()) coeff.resize(o.coeff.size());
        for (int i = 0; i < static_cast<int>(o.coeff.size()); ++i)
            coeff[i] -= o.coeff[i];
        return *this;
    }

    W eval(W x) const noexcept {
        if (coeff.empty()) return 0;
        if constexpr (std::is_same_v<Basis, MonomialBasis>) {
            W r = 0;
            for (int i = static_cast<int>(coeff.size()) - 1; i >= 0; --i)
                r = r * x + coeff[i];
            return r;
        } else if constexpr (std::is_same_v<Basis, FallingFactorialBasis>) {
            W sum = 0;
            W falling = 1;
            for (int i = 0; i < static_cast<int>(coeff.size()); ++i) {
                if (i > 0) falling = falling * (x - static_cast<W>(i - 1));
                sum += coeff[i] * falling;
            }
            return sum;
        } else if constexpr (std::is_same_v<Basis, TaylorBasis>) {
            W sum = 0;
            W b = 1;
            for (int i = 0; i < static_cast<int>(coeff.size()); ++i) {
                sum += coeff[i] * b;
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
                "DynamicPolynomial::eval: unsupported basis type");
            return W{};
        }
    }

    friend DynamicPolynomial operator+(DynamicPolynomial a, const DynamicPolynomial& b) {
        a += b;
        return a;
    }

    friend DynamicPolynomial operator-(DynamicPolynomial a, const DynamicPolynomial& b) {
        a -= b;
        return a;
    }

    DynamicPolynomial operator*(const DynamicPolynomial& o) const {
        if (coeff.empty() || o.coeff.empty()) return DynamicPolynomial{};
        int na = static_cast<int>(coeff.size());
        int nb = static_cast<int>(o.coeff.size());
        DynamicPolynomial r(na + nb - 2);
        poly_mul(r.coeff.data(), coeff.data(), na, o.coeff.data(), nb);
        return r;
    }
};

// Conversion to/from compile-time Polynomial
template<int N, std::unsigned_integral W, typename Basis>
DynamicPolynomial<W, Basis> to_dynamic(const Polynomial<N, W, Basis>& p) {
    return DynamicPolynomial<W, Basis>(std::vector<W>(p.begin(), p.end()));
}

template<int N, std::unsigned_integral W, typename Basis>
Polynomial<N, W, Basis> to_static(const DynamicPolynomial<W, Basis>& p) {
    Polynomial<N, W, Basis> r{};
    int n = std::min(N, p.size());
    for (int i = 0; i < n; ++i) r[i] = p[i];
    return r;
}

// Basis conversion for DynamicPolynomial

template<std::unsigned_integral W>
DynamicPolynomial<W, FallingFactorialBasis>
monomial_to_falling(const DynamicPolynomial<W, MonomialBasis>& p) {
    int sz = p.size();
    if (sz == 0) return {};
    detail::DynamicStirlingCache<W> cache(sz);
    DynamicPolynomial<W, FallingFactorialBasis> r(sz - 1);
    for (int k = 0; k < sz; ++k) {
        W sum = 0;
        for (int n = k; n < sz; ++n)
            sum += p[n] * cache.S2[n][k];
        r[k] = sum;
    }
    return r;
}

template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
falling_to_monomial(const DynamicPolynomial<W, FallingFactorialBasis>& p) {
    int sz = p.size();
    if (sz == 0) return {};
    detail::DynamicStirlingCache<W> cache(sz);
    DynamicPolynomial<W, MonomialBasis> r(sz - 1);
    for (int k = 0; k < sz; ++k) {
        W sum = 0;
        for (int n = k; n < sz; ++n)
            sum += p[n] * cache.s1[n][k];
        r[k] = sum;
    }
    return r;
}

template<std::unsigned_integral W>
DynamicPolynomial<W, TaylorBasis>
monomial_to_taylor(const DynamicPolynomial<W, MonomialBasis>& p) {
    int sz = p.size();
    if (sz == 0) return {};
    detail::DynamicStirlingCache<W> cache(sz);
    DynamicPolynomial<W, TaylorBasis> r(sz - 1);
    W fact = 1;
    for (int k = 0; k < sz; ++k) {
        W sum = 0;
        if (k > 0) fact *= static_cast<W>(k);
        for (int n = k; n < sz; ++n)
            sum += p[n] * cache.S2[n][k] * fact;
        r[k] = sum;
    }
    return r;
}

template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
taylor_to_monomial(const DynamicPolynomial<W, TaylorBasis>& p) {
    int sz = p.size();
    if (sz == 0) return {};
    detail::DynamicStirlingCache<W> cache(sz);
    using accum_t = detail::widen_t<dword_t<W>>;
    std::vector<W> facts(sz);
    facts[0] = 1;
    for (int i = 1; i < sz; ++i) facts[i] = facts[i-1] * static_cast<W>(i);
    DynamicPolynomial<W, MonomialBasis> r(sz - 1);
    for (int k = 0; k < sz; ++k) {
        accum_t sum{};
        for (int j = k; j < sz; ++j) {
            W fact = facts[j];
            accum_t val = static_cast<accum_t>(p[j]) *
                          static_cast<accum_t>(cache.s1[j][k]);
            if (fact % 2 == 1) {
                sum += val * static_cast<accum_t>(modinv_odd(fact));
            } else {
                int shift = v2(fact);
                W odd = fact >> shift;
                accum_t shifted = div_2k_adic(val, shift);
                sum += shifted * static_cast<accum_t>(modinv_odd(odd));
            }
        }
        r[k] = static_cast<W>(sum);
    }
    return r;
}

template<typename ToBasis, std::unsigned_integral W, typename FromBasis>
DynamicPolynomial<W, ToBasis> change_basis(const DynamicPolynomial<W, FromBasis>& p) {
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

// Formal derivative and forward difference for DynamicPolynomial

template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
formal_derivative(const DynamicPolynomial<W, MonomialBasis>& p) {
    int sz = p.size();
    if (sz <= 1) return DynamicPolynomial<W, MonomialBasis>{};
    DynamicPolynomial<W, MonomialBasis> r(sz - 2);
    for (int i = 1; i < sz; ++i)
        r[i - 1] = static_cast<W>(i) * p[i];
    return r;
}

template<std::unsigned_integral W>
DynamicPolynomial<W, FallingFactorialBasis>
formal_derivative(const DynamicPolynomial<W, FallingFactorialBasis>& p) {
    return change_basis<FallingFactorialBasis>(
        formal_derivative(change_basis<MonomialBasis>(p)));
}

template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
forward_difference(const DynamicPolynomial<W, MonomialBasis>& p) {
    int sz = p.size();
    if (sz <= 1) return DynamicPolynomial<W, MonomialBasis>{};
    DynamicPolynomial<W, MonomialBasis> r(sz - 2);
    std::vector<W> C(sz);
    C[0] = 1;
    for (int n = 1; n < sz; ++n) {
        for (int k = n; k >= 1; --k)
            C[k] = C[k-1] + C[k];
        int lim = n - 1;
        if (lim > sz - 2) lim = sz - 2;
        for (int i = 0; i <= lim; ++i)
            r[i] += p[n] * C[i];
    }
    return r;
}

template<std::unsigned_integral W>
DynamicPolynomial<W, FallingFactorialBasis>
forward_difference(const DynamicPolynomial<W, FallingFactorialBasis>& p) {
    int sz = p.size();
    if (sz <= 1) return DynamicPolynomial<W, FallingFactorialBasis>{};
    DynamicPolynomial<W, FallingFactorialBasis> r(sz - 2);
    for (int i = 1; i < sz; ++i)
        r[i - 1] = static_cast<W>(i) * p[i];
    return r;
}

// Semantic alias for the six-axiom naming scheme.
template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
exact_formal_derivative(const DynamicPolynomial<W, MonomialBasis>& p) {
    return formal_derivative(p);
}

template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
exact_forward_difference(const DynamicPolynomial<W, MonomialBasis>& p) {
    return forward_difference(p);
}

} // namespace dyadic
