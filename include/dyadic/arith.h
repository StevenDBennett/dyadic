// dyadic/arith.h — Carry-save arithmetic and Newton division
// Requires: dyadic/core.h
//
// All operations treat Polynomial<N,W,MonomialBasis> as an unsigned big
// integer (little-endian limbs).  This is the "carry-chain ring" view where
// coefficient arrays behave like N-limb base-2^W integers.
//
// Newton division provides O(N log N) big-integer division via reciprocal
// iteration; the correction step uses Knuth D for the final adjustment.

#pragma once

#include <dyadic/core.h>
#include <utility>

namespace dyadic {

// ============================================================================
// Internal helpers for big-integer arithmetic on Polynomial
// ============================================================================

namespace detail {

template<int M, int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<M, W, Basis>
zero_extend(const Polynomial<N, W, Basis>& p) noexcept {
    static_assert(M >= N, "zero_extend: target size must be >= source size");
    Polynomial<M, W, Basis> r{};
    for (int i = 0; i < N; ++i) r[i] = p[i];
    return r;
}

template<int N, std::unsigned_integral W, typename Basis>
constexpr bool
unsigned_lt(const Polynomial<N, W, Basis>& a, const Polynomial<N, W, Basis>& b) noexcept {
    for (int i = N - 1; i >= 0; --i) {
        if (a[i] != b[i]) return a[i] < b[i];
    }
    return false;
}

template<int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N, W, Basis>
shift_left(const Polynomial<N, W, Basis>& p, int shift) noexcept {
    if (shift <= 0) return p;
    constexpr int BITS = 8 * sizeof(W);
    int wshift = shift / BITS;
    int bshift = shift % BITS;
    Polynomial<N, W, Basis> r{};
    for (int i = N - 1; i >= wshift; --i) r[i] = p[i - wshift];
    if (bshift > 0) {
        W carry = 0;
        for (int i = 0; i < N; ++i) {
            W old = r[i];
            r[i] = static_cast<W>((old << bshift) | carry);
            carry = old >> (BITS - bshift);
        }
    }
    return r;
}

// Full-width unsigned multiplication returning N+M limbs.
// Uses quad_width accumulators for unsaturated convolution, then a
// full carry-chain that can produce a carry into the (N+M-1)-th limb.
template<int N, int M, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N + M, W, Basis>
mul_unsigned(const Polynomial<N, W, Basis>& a,
             const Polynomial<M, W, Basis>& b) noexcept {
    constexpr int R = N + M;
    using accum_t = quad_width<W>;
    std::array<accum_t, R> buf{};
    for (int i = 0; i < N; ++i) {
        accum_t ai = static_cast<accum_t>(a[i]);
        if (ai == 0) continue;
        for (int j = 0; j < M; ++j) {
            buf[i + j] += ai * static_cast<accum_t>(b[j]);
        }
    }
    Polynomial<R, W, Basis> r{};
    accum_t carry = 0;
    for (int i = 0; i < R; ++i) {
        accum_t sum = buf[i] + carry;
        r[i] = static_cast<W>(sum);
        carry = sum >> (8 * sizeof(W));
    }
    return r;
}

template<std::unsigned_integral W>
constexpr int clz_normalize(W x) noexcept {
    return std::countl_zero(x);
}

// Single-limb division: multi-limb dividend / single-word divisor.
// Returns (quotient, remainder word).
template<int NL, std::unsigned_integral W, typename Basis>
constexpr std::pair<Polynomial<NL, W, Basis>, W>
divmod_single(const Polynomial<NL, W, Basis>& u, W v) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);
    Polynomial<NL, W, Basis> q{};
    dw_t carry = 0;
    dw_t vv = static_cast<dw_t>(v);
    for (int i = NL - 1; i >= 0; --i) {
        carry = (carry << BITS) | static_cast<dw_t>(u[i]);
        dw_t quot = carry / vv;
        q[i] = static_cast<W>(quot);
        carry -= quot * vv;
    }
    return {q, static_cast<W>(carry)};
}

// Knuth Algorithm D for unsigned division.
// Precondition: v[NR-1] != 0, NR >= 2.
template<int NL, int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<NL, W, MonomialBasis>, Polynomial<NR, W, MonomialBasis>>
div_unsigned_knuth(const Polynomial<NL, W, MonomialBasis>& u,
                   const Polynomial<NR, W, MonomialBasis>& v) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);
    constexpr int n = NR;
    constexpr int m = (NL >= NR) ? (NL - NR) : 0;

    if constexpr (NL < NR) {
        Polynomial<NL, W, MonomialBasis> q{};
        Polynomial<NR, W, MonomialBasis> r{};
        for (int i = 0; i < NL; ++i) r[i] = u[i];
        return {q, r};
    }

    int d = clz_normalize(v[n - 1]);

    constexpr int UN_SIZE = n + m + 1;
    std::array<W, UN_SIZE> un{};
    for (int i = 0; i < NL && i < UN_SIZE; ++i) un[i] = u[i];
    if (d > 0) {
        W carry = 0;
        for (int i = 0; i < UN_SIZE; ++i) {
            W old = un[i];
            un[i] = static_cast<W>((static_cast<dw_t>(old) << d) | static_cast<dw_t>(carry));
            carry = static_cast<W>(static_cast<dw_t>(old) >> (BITS - d));
        }
    }

    Polynomial<n, W, MonomialBasis> vn{};
    if (d > 0) {
        W carry = 0;
        for (int i = 0; i < n; ++i) {
            W old = v[i];
            vn[i] = static_cast<W>((static_cast<dw_t>(old) << d) | carry);
            carry = static_cast<W>(static_cast<dw_t>(old) >> (BITS - d));
        }
    } else {
        for (int i = 0; i < n; ++i) vn[i] = v[i];
    }

    std::array<W, m + 1> q{};

    for (int j = m; j >= 0; --j) {
        dw_t uhat = (static_cast<dw_t>(un[j + n]) << BITS) | static_cast<dw_t>(un[j + n - 1]);
        W qhat, rhat;

        if (static_cast<dw_t>(un[j + n]) >= static_cast<dw_t>(vn[n - 1])) {
            qhat = ~W(0);
            rhat = static_cast<W>(uhat - static_cast<dw_t>(qhat) * static_cast<dw_t>(vn[n - 1]));
            rhat += vn[n - 1];
        } else {
            qhat = static_cast<W>(uhat / static_cast<dw_t>(vn[n - 1]));
            rhat = static_cast<W>(uhat - static_cast<dw_t>(qhat) * static_cast<dw_t>(vn[n - 1]));
        }

        while (true) {
            dw_t lhs = static_cast<dw_t>(qhat) * static_cast<dw_t>(vn[n - 2]);
            dw_t rhs = (static_cast<dw_t>(rhat) << BITS) | static_cast<dw_t>(un[j + n - 2]);
            if (lhs <= rhs) break;
            qhat--;
            W old_rhat = rhat;
            rhat += vn[n - 1];
            if (rhat < old_rhat) break;
        }

        dw_t mul_carry = 0;
        dw_t sub_borrow = 0;
        for (int i = 0; i < n; ++i) {
            dw_t prod = static_cast<dw_t>(qhat) * static_cast<dw_t>(vn[i]) + mul_carry;
            mul_carry = prod >> BITS;
            dw_t diff = static_cast<dw_t>(un[j + i]) - static_cast<W>(prod) - sub_borrow;
            un[j + i] = static_cast<W>(diff);
            sub_borrow = (diff >> (2 * BITS - 1)) & 1;
        }
        {
            dw_t diff = static_cast<dw_t>(un[j + n]) - static_cast<W>(mul_carry) - sub_borrow;
            un[j + n] = static_cast<W>(diff);
            W borrow = static_cast<W>((diff >> (2 * BITS - 1)) & 1);

            q[j] = qhat;
            if (borrow) {
                q[j]--;
                dw_t add_carry = 0;
                for (int i = 0; i < n; ++i) {
                    dw_t sum = static_cast<dw_t>(un[j + i]) + static_cast<dw_t>(vn[i]) + add_carry;
                    un[j + i] = static_cast<W>(sum);
                    add_carry = sum >> BITS;
                }
                un[j + n] = static_cast<W>(static_cast<dw_t>(un[j + n]) + add_carry);
            }
        }
    }

    Polynomial<NR, W, MonomialBasis> r;
    if (d > 0) {
        W carry = 0;
        for (int i = n - 1; i >= 0; --i) {
            dw_t val = (static_cast<dw_t>(un[i]) >> d) | (static_cast<dw_t>(carry) << (BITS - d));
            r[i] = static_cast<W>(val);
            carry = un[i];
        }
    } else {
        for (int i = 0; i < n; ++i) r[i] = un[i];
    }

    Polynomial<NL, W, MonomialBasis> q_full{};
    for (int i = 0; i < NL && i < static_cast<int>(m + 1); ++i) q_full[i] = q[i];

    return {q_full, r};
}

} // namespace detail

// ============================================================================
// Carry-save arithmetic (Polynomial as big integer)
// ============================================================================

template<int N, std::unsigned_integral W, typename Basis>
struct CarrySavePair {
    Polynomial<N, W, Basis> sum;
    Polynomial<N + 1, W, Basis> carry;
};

// Decompose a + b into (sum, carry) such that a + b = sum + carry (as values).
// carry[i] ∈ {0, 1} for i > 0, carry[0] = 0.
template<int N, std::unsigned_integral W, typename Basis>
constexpr CarrySavePair<N, W, Basis>
carry_save_add(const Polynomial<N, W, Basis>& a,
               const Polynomial<N, W, Basis>& b) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);
    CarrySavePair<N, W, Basis> r;
    r.carry[0] = 0;
    for (int i = 0; i < N; ++i) {
        dw_t t = static_cast<dw_t>(a[i]) + static_cast<dw_t>(b[i]);
        r.sum[i] = static_cast<W>(t);
        r.carry[i + 1] = static_cast<W>(t >> BITS);
    }
    return r;
}

// 3→2 carry-save adder: reduce (a, b, c) to (sum, carry)
// such that a + b + c = sum + carry (as values).  carry[i] ∈ {0, 1, 2}.
template<int N, std::unsigned_integral W, typename Basis>
constexpr CarrySavePair<N, W, Basis>
csa3(const Polynomial<N, W, Basis>& a,
     const Polynomial<N, W, Basis>& b,
     const Polynomial<N, W, Basis>& c) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);
    CarrySavePair<N, W, Basis> r;
    r.carry[0] = 0;
    for (int i = 0; i < N; ++i) {
        dw_t t = static_cast<dw_t>(a[i]) + static_cast<dw_t>(b[i]) + static_cast<dw_t>(c[i]);
        r.sum[i] = static_cast<W>(t);
        r.carry[i + 1] = static_cast<W>(t >> BITS);
    }
    return r;
}

// Combine a carry-save pair back into a normal Polynomial,
// propagating all carries.  The final overflow beyond N limbs
// is discarded.
template<int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N, W, Basis>
combine_carry_save(const Polynomial<N, W, Basis>& sum,
                   const Polynomial<N + 1, W, Basis>& carry) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);
    Polynomial<N, W, Basis> r;
    dw_t acc = 0;
    for (int i = 0; i < N; ++i) {
        acc += static_cast<dw_t>(carry[i]) + static_cast<dw_t>(sum[i]);
        r[i] = static_cast<W>(acc);
        acc >>= BITS;
    }
    return r;
}

// Batched addition using carry-save accumulation.
// Sums all operands through a single final carry propagation.
template<int N, std::unsigned_integral W, typename Basis>
constexpr Polynomial<N, W, Basis>
batched_add_carry_save(const Polynomial<N, W, Basis>* operands,
                       int count) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);

    if (count <= 0) return Polynomial<N, W, Basis>{};
    if (count == 1) return operands[0];

    Polynomial<N, W, Basis> sum = operands[0];
    Polynomial<N + 1, W, Basis> carry{};
    carry[0] = 0;

    for (int j = 1; j < count; ++j) {
        for (int i = 0; i < N; ++i) {
            dw_t t = static_cast<dw_t>(sum[i]) + static_cast<dw_t>(operands[j][i]);
            sum[i] = static_cast<W>(t);
            carry[i + 1] += static_cast<W>(t >> BITS);
        }
    }

    return combine_carry_save(sum, carry);
}

// ============================================================================
// Unsigned division wrapper
// ============================================================================

template<int NL, int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<NL, W, MonomialBasis>, Polynomial<NR, W, MonomialBasis>>
div_unsigned(const Polynomial<NL, W, MonomialBasis>& u,
             const Polynomial<NR, W, MonomialBasis>& v) noexcept {
    int eff_nr = NR;
    while (eff_nr > 0 && v[eff_nr - 1] == 0) eff_nr--;

    if (eff_nr == 0) return {};

    if (eff_nr == 1) {
        auto [q, rem] = detail::divmod_single(u, v[eff_nr - 1]);
        Polynomial<NR, W, MonomialBasis> r{};  // zero-initialize all limbs
        r[0] = rem;
        return {q, r};
    }

    if constexpr (NR > 2) {
        if (eff_nr < NR) {
            Polynomial<NR - 1, W, MonomialBasis> v_eff;
            for (int i = 0; i < NR - 1; ++i) v_eff[i] = v[i];
            auto [q_tmp, r_tmp] = div_unsigned<NL, NR - 1>(u, v_eff);
            Polynomial<NR, W, MonomialBasis> r_full;
            for (int i = 0; i < NR - 1; ++i) r_full[i] = r_tmp[i];
            return {q_tmp, r_full};
        }
    }

    return detail::div_unsigned_knuth(u, v);
}

// ============================================================================
// Newton division — reciprocal via Newton iteration
// ============================================================================
// Newton's method for reciprocal: x_{n+1} = x_n * (2 - d * x_n).
// The scaled integer form (r ≈ B^(2n)/d) with n = NR limbs:
//   r_{k+1} = r_k * (2*B^(2n) - d * r_k) / B^(2n)

// Compute reciprocal r ≈ B^(2*NR) / d where d has NR limbs.
// d will be normalized internally if needed (top bit of d[NR-1] set).
// Returns (r, norm_shift) where r has 2*NR limbs.
template<int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<2 * NR, W, MonomialBasis>, int>
reciprocal_newton(const Polynomial<NR, W, MonomialBasis>& d) noexcept {
    constexpr int TWO_NR = 2 * NR;
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);

    Polynomial<NR, W, MonomialBasis> d_work;
    int norm_shift;
    {
        int top_idx = NR - 1;
        while (top_idx >= 0 && d[top_idx] == 0) top_idx--;

        int word_shift = 0;
        int bit_shift = 0;
        if (top_idx >= 0) {
            bit_shift = detail::clz_normalize(d[top_idx]);
            word_shift = NR - 1 - top_idx;
        }
        norm_shift = word_shift * BITS + bit_shift;

        for (int i = 0; i < NR; ++i) {
            d_work[i] = (i >= word_shift) ? d[i - word_shift] : W(0);
        }

        if (bit_shift > 0) {
            W carry = 0;
            for (int i = 0; i < NR; ++i) {
                W old = d_work[i];
                d_work[i] = static_cast<W>((static_cast<dw_t>(old) << bit_shift) | carry);
                carry = static_cast<W>(static_cast<dw_t>(old) >> (BITS - bit_shift));
            }
        }
    }

    Polynomial<TWO_NR, W, MonomialBasis> r{};

    W d_hi = d_work[NR - 1];
    dw_t init_full = (dw_t(1) << BITS) / static_cast<dw_t>(d_hi);
    W init = static_cast<W>(init_full);
    r[NR] = init;

    auto d_ext = detail::zero_extend<TWO_NR>(d_work);

    int bits_needed = (NR + 1) * BITS;
    int max_iter = 1;
    int bp = 2;
    while (bp < bits_needed) { bp *= 2; max_iter++; }

    for (int iter = 0; iter < max_iter; ++iter) {
        auto t = detail::mul_unsigned(d_ext, r);

        Polynomial<4 * NR, W, MonomialBasis> two_B{};
        two_B[TWO_NR] = W(2);

        Polynomial<4 * NR, W, MonomialBasis> error;
        dw_t carry = 1;
        for (int i = 0; i < 4 * NR; ++i) {
            carry = carry + static_cast<dw_t>(~t[i]) + static_cast<dw_t>(two_B[i]);
            error[i] = static_cast<W>(carry);
            carry >>= BITS;
        }

        auto prod = detail::mul_unsigned(r, error);

        for (int i = 0; i < TWO_NR; ++i)
            r[i] = prod[TWO_NR + i];
    }

    return {r, norm_shift};
}

// div_newton — unsigned division using Newton reciprocal.
// Computes (quotient, remainder) = u / v.
template<int NL, int NR, std::unsigned_integral W>
constexpr std::pair<Polynomial<NL, W, MonomialBasis>, Polynomial<NR, W, MonomialBasis>>
div_newton(const Polynomial<NL, W, MonomialBasis>& u,
           const Polynomial<NR, W, MonomialBasis>& v) noexcept {
    using dw_t = dword_t<W>;
    constexpr int BITS = 8 * sizeof(W);
    constexpr int TWO_NR = 2 * NR;

    // Handle leading zero limbs in divisor by falling through to div_unsigned
    {
        int eff_nr = NR;
        while (eff_nr > 0 && v[eff_nr - 1] == 0) eff_nr--;
        if (eff_nr < NR) {
            return div_unsigned(u, v);
        }
    }

    if constexpr (NR == 1) {
        auto [q, rem] = detail::divmod_single(u, v[0]);
        Polynomial<NR, W, MonomialBasis> r{};
        r[0] = rem;
        return {q, r};
    }

    {
        bool u_smaller = false;
        bool found_diff = false;
        int max_limbs = (NL > NR) ? NL : NR;
        for (int i = max_limbs - 1; i >= 0; --i) {
            W uv = (i < NL) ? u[i] : W(0);
            W vv = (i < NR) ? v[i] : W(0);
            if (uv != vv) {
                u_smaller = (uv < vv);
                found_diff = true;
                break;
            }
        }
        if (found_diff && u_smaller) {
            Polynomial<NL, W, MonomialBasis> q{};
            Polynomial<NR, W, MonomialBasis> r{};
            for (int i = 0; i < NL && i < NR; ++i) r[i] = u[i];
            return {q, r};
        }
    }

    auto [recip_norm, norm_shift] = reciprocal_newton(v);

    int wshift = norm_shift / BITS;
    int bshift = norm_shift % BITS;
    constexpr int U_EXT = NL + NR;
    Polynomial<U_EXT, W, MonomialBasis> u_shifted;
    for (int i = 0; i < U_EXT; ++i) {
        int src = (i >= wshift) ? i - wshift : -1;
        u_shifted[i] = (src >= 0 && src < NL) ? u[src] : W(0);
    }
    if (bshift > 0) {
        W carry = 0;
        for (int i = 0; i < U_EXT; ++i) {
            W old = u_shifted[i];
            u_shifted[i] = static_cast<W>((static_cast<dw_t>(old) << bshift) | carry);
            carry = static_cast<W>(static_cast<dw_t>(old) >> (BITS - bshift));
        }
    }

    auto prod = detail::mul_unsigned(u_shifted, recip_norm);

    Polynomial<NL, W, MonomialBasis> q;
    for (int i = 0; i < NL; ++i)
        q[i] = prod[TWO_NR + i];

    auto q_times_v = detail::mul_unsigned(q, v);
    auto u_ext = detail::zero_extend<NL + NR>(u);

    Polynomial<NL + NR, W, MonomialBasis> rem;
    dw_t carry = 1;
    for (int i = 0; i < NL + NR; ++i) {
        carry = carry + static_cast<dw_t>(u_ext[i]) + static_cast<dw_t>(~q_times_v[i]);
        rem[i] = static_cast<W>(carry);
        carry >>= BITS;
    }

    if (carry == 0) {
        for (int i = 0; i < NL; ++i) {
            if (q[i] > 0) { q[i]--; break; }
            q[i] = ~W(0);
        }
        dw_t add_carry = 0;
        for (int i = 0; i < NL + NR; ++i) {
            W v_limb = (i < NR) ? v[i] : W(0);
            add_carry += static_cast<dw_t>(rem[i]) + static_cast<dw_t>(v_limb);
            rem[i] = static_cast<W>(add_carry);
            add_carry >>= BITS;
        }
    }

    {
        bool rem_ge_v = true;
        for (int i = NL + NR - 1; i >= 0; --i) {
            W r_limb = rem[i];
            W v_limb = (i < NR) ? v[i] : W(0);
            if (r_limb != v_limb) {
                rem_ge_v = (r_limb > v_limb);
                break;
            }
        }
        if (rem_ge_v) {
            auto [k, rem_new] = div_unsigned<NL + NR, NR>(rem, v);

            dw_t addc = 0;
            for (int i = 0; i < NL; ++i) {
                W kw = (i < NL + NR) ? k[i] : W(0);
                addc += static_cast<dw_t>(q[i]) + static_cast<dw_t>(kw);
                q[i] = static_cast<W>(addc);
                addc >>= BITS;
            }

            for (int i = 0; i < NL + NR; ++i)
                rem[i] = (i < NR) ? rem_new[i] : W(0);
        }
    }

    Polynomial<NR, W, MonomialBasis> r_final;
    for (int i = 0; i < NR; ++i) r_final[i] = rem[i];

    return {q, r_final};
}

} // namespace dyadic
