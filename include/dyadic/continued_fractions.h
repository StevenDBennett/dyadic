// continued_fractions.h — Continued Fraction Expansion (Dyadic extension)
// Requires: dyadic/dynamic_polynomial.h

#pragma once

#include <dyadic/dynamic_polynomial.h>

#include <utility>
#include <vector>

namespace dyadic {

template<std::unsigned_integral W>
std::vector<W> cf_expand(const DynamicPolynomial<W, MonomialBasis>& series, int max_terms) {
    int n = std::min(max_terms, series.size());
    std::vector<W> result;
    result.reserve(n);
    for (int i = 0; i < n; ++i)
        result.push_back(series[i]);
    return result;
}

template<std::unsigned_integral W>
std::pair<DynamicPolynomial<W, MonomialBasis>, DynamicPolynomial<W, MonomialBasis>>
cf_convergent(const std::vector<W>& c, int n) {
    using DP = DynamicPolynomial<W, MonomialBasis>;
    if (n < 0 || c.empty()) return {DP{}, DP{}};
    std::vector<W> Pm1 = {1};
    std::vector<W> P0  = {c[0]};
    std::vector<W> Qm1 = {0};
    std::vector<W> Q0  = {1};
    for (int k = 1; k <= n && k < (int)c.size(); ++k) {
        int deg = std::max((int)P0.size(), (int)Pm1.size() + 1);
        std::vector<W> Pk(deg, 0);
        for (int i = 0; i < (int)P0.size(); ++i)
            Pk[i] = c[k] * P0[i];
        for (int i = 0; i < (int)Pm1.size(); ++i)
            Pk[i + 1] += Pm1[i];
        std::vector<W> Qk(deg, 0);
        for (int i = 0; i < (int)Q0.size(); ++i)
            Qk[i] = c[k] * Q0[i];
        for (int i = 0; i < (int)Qm1.size(); ++i)
            Qk[i + 1] += Qm1[i];
        Pm1 = std::move(P0);
        P0  = std::move(Pk);
        Qm1 = std::move(Q0);
        Q0  = std::move(Qk);
    }
    return {DP(P0), DP(Q0)};
}

// Returns a DynamicPolynomial containing a single constant term (the evaluated
// convergent at x). This is not a full polynomial — the return type is an artifact
// of reusing cf_convergent's internal machinery.
template<std::unsigned_integral W>
DynamicPolynomial<W, MonomialBasis>
cf_eval(const std::vector<W>& c, int max_terms, W x) {
    using DP = DynamicPolynomial<W, MonomialBasis>;
    auto [P, Q] = cf_convergent(c, max_terms);
    if (Q.eval(x) % 2 == 1)
        return DP({P.eval(x) * modinv_odd(Q.eval(x))});
    return DP({});
}

} // namespace dyadic
