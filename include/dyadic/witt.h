// dyadic/witt.h — Witt Vectors, Ghost Map, W_add/W_mul, W_log/W_exp/W_inverse
// Requires: dyadic/core.h
//
// Witt vector arithmetic with ghost-map and Newton recovery:
//   - Ghost map (Frobenius-module embedding)
//   - Witt addition and multiplication (ghost-op + recovery)
//   - Adams operation (Teichmüller-style ghost power)
//   - Witt logarithm, exponential, inverse (via p-adic series)
//   - Verification helpers

#pragma once

#include <dyadic/core.h>

namespace dyadic {

// ============================================================================
// 13. Witt Vectors
// ============================================================================

namespace detail {

// Precomputed powers a[i]^(2^k) for the ghost map.
// Reduces ghost sum from O(N³) to O(N²).
template<int N, std::unsigned_integral W>
struct GhostPowers {
    std::array<std::array<W, N>, N> pow{};

    constexpr void init(const W* a) noexcept {
        for (int i = 0; i < N; ++i) build_row(i, a[i]);
    }

    constexpr void set(int i, W val) noexcept {
        build_row(i, val);
    }

    constexpr W get(int i, int k) const noexcept { return pow[i][k]; }

  private:
    constexpr void build_row(int i, W val) noexcept {
        pow[i][0] = val;
        for (int k = 1; k < N - i; ++k) {
            pow[i][k] = static_cast<W>(1u * pow[i][k-1] * pow[i][k-1]);
        }
    }
};

// GhostPowersFull — gw_t-precision precomputed powers for ghost_map.
// Like GhostPowers but stores at gw_t (= widen_t<dword_t<W>>) precision
// Single O(N²) init.  Stores at gw_t (= quad_width<W>) precision to avoid
// truncation in log/exp/inverse recovery.
template<int N, std::unsigned_integral W>
struct GhostPowersFull {
    using gw_t = quad_width<W>;
    std::array<std::array<gw_t, N>, N> pow{};

    constexpr void init(const W* a) noexcept {
        for (int i = 0; i < N; ++i) {
            pow[i][0] = static_cast<gw_t>(a[i]);
            for (int k = 1; k < N - i; ++k) {
                pow[i][k] = pow[i][k-1] * pow[i][k-1];
            }
        }
    }

    constexpr gw_t ghost(int j) const noexcept {
        gw_t sum{};
        for (int i = 0; i <= j; ++i) {
            sum += (gw_t(uint64_t(1)) << i) * pow[i][j - i];
        }
        return sum;
    }
};

// Ghost sum using quad_width<W> (widen_t<dword_t<W>>) for extra headroom.
// For W=uint8_t this is uint32_t (4×); for larger types it's the next
// level in the 2× chain (dword_t widened once more).
template<int N, std::unsigned_integral W>
constexpr auto ghost_j_widen(const GhostPowers<N, W>& pows, int j) noexcept {
    using gw_t = quad_width<W>;
    gw_t sum{};
    for (int i = 0; i <= j; ++i) {
        sum += (gw_t(uint64_t(1)) << i) *
               static_cast<gw_t>(pows.get(i, j - i));
    }
    return sum;
}

// Wrapper from quad_width<W> → dword_t<W>.
// Ghost_j_widen accumulates at quad_width precision (4× headroom for uint8_t);
// this wrapper narrows to dword_t<W> for the WittVector::ghost() API,
// which further narrows to W on return. The dword_t intermediate prevents
// spurious overflow in ghost-compare operations.
template<int N, std::unsigned_integral W>
constexpr dword_t<W> ghost_j_dword(const GhostPowers<N, W>& pows, int j) noexcept {
    return static_cast<dword_t<W>>(ghost_j_widen(pows, j));
}

} // namespace detail

template<int N, std::unsigned_integral W>
struct WittVector {
    std::array<W, N> a{};

    constexpr WittVector() noexcept = default;
    constexpr WittVector(std::array<W, N> c) noexcept : a(c) {}

    constexpr W& operator[](int i) noexcept { return a[i]; }
    constexpr W operator[](int i) const noexcept { return a[i]; }

    constexpr W ghost(int j) const noexcept {
        return ghost_vector()[j];
    }

    constexpr std::array<W, N> ghost_vector() const noexcept {
        detail::GhostPowers<N, W> pows;
        pows.init(a.data());
        std::array<W, N> g{};
        for (int j = 0; j < N; ++j) {
            g[j] = static_cast<W>(detail::ghost_j_dword(pows, j));
        }
        return g;
    }

    constexpr WittVector frobenius() const noexcept {
        WittVector r;
        for (int i = 0; i < N; ++i) r.a[i] = static_cast<W>(1u * a[i] * a[i]);
        return r;
    }

    constexpr WittVector verschiebung() const noexcept {
        WittVector r;
        r.a[0] = 0;
        for (int i = 1; i < N; ++i) r.a[i] = a[i-1];
        return r;
    }

    constexpr WittVector FV() const noexcept {
        return frobenius().verschiebung();
    }

    constexpr WittVector VF() const noexcept {
        return verschiebung().frobenius();
    }

    constexpr friend WittVector operator+(const WittVector& a, const WittVector& b) noexcept {
        return witt_add(a, b);
    }

    constexpr friend WittVector operator*(const WittVector& a, const WittVector& b) noexcept {
        return witt_mul(a, b);
    }

    static constexpr WittVector zero() noexcept { return WittVector{}; }
    static constexpr WittVector one() noexcept {
        std::array<W, N> arr{}; arr[0] = W(1);
        return WittVector(arr);
    }
};

template<int N, std::unsigned_integral W>
constexpr bool check_ghost_ring(const WittVector<N,W>& a, const WittVector<N,W>& b) noexcept {
    auto ga = a.ghost_vector();
    auto gb = b.ghost_vector();
    auto gc = witt_add(a, b).ghost_vector();
    for (int j = 0; j < N; ++j) {
        if (gc[j] != ga[j] + gb[j]) return false;
    }
    return true;
}

// ============================================================================
// 18–19. Witt Addition & Multiplication (ghost-map + Newton recovery) — O(N²)
// ============================================================================

namespace detail {

// Compute ghost_j at full gw_t precision from raw Witt components.
// Unlike ghost_j_widen (which uses W-precision GhostPowers), this computes
// each a[i]^(2^{j-i}) at gw_t precision, avoiding truncation. Required when
// Witt components are large (e.g., after inverse/log/exp recovery).
template<int N, std::unsigned_integral W>
constexpr auto ghost_j_full(const W* a, int j) noexcept {
    using gw_t = quad_width<W>;
    gw_t sum{};
    for (int i = 0; i <= j; ++i) {
        gw_t pow = static_cast<gw_t>(a[i]);
        for (int k = 0; k < j - i; ++k) {
            pow = pow * pow;
        }
        sum += (gw_t(uint64_t(1)) << i) * pow;
    }
    return sum;
}

// Newton recovery from ghost values: given G[j] (widened precision), recover
// Witt components r[j] via r_j = (G_j - S_j) / 2^j where
// S_j = Σ_{i<j} 2^i · r_i^{2^{j-i}}.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> ghost_recover(const std::array<quad_width<W>, N>& G) noexcept {
    static_assert(N <= 32, "ghost_recover: N too large, stack overflow risk (N^2 * sizeof(gw_t))");
    using gw_t = quad_width<W>;
    WittVector<N, W> r;
    r[0] = static_cast<W>(G[0]);

    // Precompute r[i]^(2^k) at gw_t precision to avoid truncation
    std::array<std::array<gw_t, N>, N> pows{};
    pows[0][0] = static_cast<gw_t>(r[0]);
    for (int k = 1; k < N; ++k) {
        pows[0][k] = pows[0][k-1] * pows[0][k-1];
    }

    for (int j = 1; j < N; ++j) {
        gw_t S{};
        for (int i = 0; i < j; ++i) {
            S += (gw_t(uint64_t(1)) << i) * pows[i][j - i];
        }

        gw_t diff = G[j] - S;
        r[j] = static_cast<W>(diff >> j);

        pows[j][0] = static_cast<gw_t>(r[j]);
        for (int k = 1; k < N - j; ++k) {
            pows[j][k] = pows[j][k-1] * pows[j][k-1];
        }
    }

    return r;
}

template<int N, std::unsigned_integral W, typename Combine>
constexpr WittVector<N, W> ghost_op(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b,
                                    Combine combine) noexcept {
    static_assert(N <= 32, "ghost_op: N > 32 causes large stack allocation in ghost_recover");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> a_pows, b_pows;
    a_pows.init(a.a.data());
    b_pows.init(b.a.data());

    std::array<gw_t, N> G{};
    for (int j = 0; j < N; ++j) {
        G[j] = combine(a_pows.ghost(j), b_pows.ghost(j));
    }

    return ghost_recover<N, W>(G);
}

} // namespace detail

template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_add(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b) noexcept {
    return detail::ghost_op<N>(a, b, [](auto x, auto y) { return x + y; });
}

template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_mul(const WittVector<N, W>& a,
                                    const WittVector<N, W>& b) noexcept {
    return detail::ghost_op<N>(a, b, [](auto x, auto y) { return x * y; });
}

// Adams operation ψ^n: ghost_j(ψ^n(a)) = ghost_j(a)^n (componentwise power).
// A ring endomorphism on Witt vectors; ψ^p = Frobenius for prime p.
// ψ^{mn} = ψ^m ∘ ψ^n, ψ^1 = identity. Requires n ≥ 1.
// Ghost power uses quad_width<W> for extra headroom (4× for uint8_t,
// 2× of dword_t for others). For W=uint64, dword_t=uint128_t and no wider
// type exists, so overflow may occur for large ghost values with n ≥ 2.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> adams_operation(const WittVector<N, W>& a, int n) noexcept {
    static_assert(N <= 32, "adams_operation: N > 32 causes large stack allocation in ghost_recover");
    if (n <= 1) return a;
    assert(n <= 65536 && "adams_operation: n out of reasonable range");

    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> G{};
    for (int j = 0; j < N; ++j) {
        gw_t g = pows.ghost(j);
        gw_t gn = gw_t(uint64_t(1));
        for (int i = 0; i < n; ++i) gn = gn * g;
        G[j] = gn;
    }
    return detail::ghost_recover<N, W>(G);
}

// Teichmüller lift: τ(x) = (x, 0, ..., 0).
// Embeds the residue field into the Witt ring; τ(xy) = τ(x) × τ(y).
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> teichmueller_lift(W x) noexcept {
    WittVector<N, W> r{};
    r[0] = x;
    return r;
}

// ============================================================================
// 22. Witt Vector Ghost Recovery Precision Check
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr bool check_witt_recovery_precision(const WittVector<N, W>& w) noexcept {
    static_assert(N <= 32, "check_witt_recovery_precision: N > 32 causes large stack allocation");
    auto gv = w.ghost_vector();
    std::array<quad_width<W>, N> G{};
    for (int j = 0; j < N; ++j) G[j] = static_cast<quad_width<W>>(gv[j]);
    auto recovered = detail::ghost_recover<N, W>(G);
    auto rgv = recovered.ghost_vector();
    for (int j = 0; j < N; ++j) {
        if (gv[j] != rgv[j]) return false;
    }
    return true;
}

// ============================================================================
// 23. Witt Logarithm, Exponential, and Inverse
// ============================================================================

namespace detail {

// 2-adic primitives for uint128_t (not std::unsigned_integral)
// Named with _128 suffix to avoid hiding dyadic::v2, modinv_odd, div_2k_adic.
constexpr int v2_128(uint128_t x) noexcept {
    if (x.lo != 0) return std::countr_zero(x.lo);
    if (x.hi != 0) return 64 + std::countr_zero(x.hi);
    return 128;
}

constexpr uint128_t modinv_odd_128(uint128_t a) noexcept {
    uint128_t x = 1;
    constexpr int iterations = std::bit_width(128u) - 1;
    for (int i = 0; i < iterations; ++i) {
        x = x * (uint128_t(2) - a * x);
    }
    return x;
}

constexpr uint128_t div_2k_adic_128(uint128_t val, int k) noexcept {
    if (k <= 0) return val;
    bool negative = (val.hi >> 63) != 0;
    val = val >> k;
    if (negative) {
        uint128_t mask = ~(uint128_t::all_ones() >> k);
        val = val | mask;
    }
    return val;
}

// Divide a ℤ₂ value (generic type T) by integer n using 2-adic division.
template<typename T>
constexpr T div_by_int_generic(T val, int n) noexcept {
    T n_t = static_cast<T>(static_cast<uint64_t>(n));
    int shift;
    T odd_part;
    if constexpr (std::is_same_v<T, uint128_t>) {
        shift = v2_128(n_t);
        odd_part = n_t >> shift;
        return div_2k_adic_128(val, shift) * modinv_odd_128(odd_part);
    } else {
        shift = dyadic::v2(n_t);
        odd_part = n_t >> shift;
        return dyadic::div_2k_adic(val, shift) * dyadic::modinv_odd(odd_part);
    }
}

// 2-adic valuation for the generic types used in p-adic series.
template<typename T>
constexpr int p_adic_val(T x) noexcept {
    if constexpr (std::is_same_v<T, uint128_t>) {
        return v2_128(x);
    } else {
        return dyadic::v2(x);
    }
}

// ℤ₂ logarithm: log(1+y) = Σ_{n=1}^{∞} (-1)^{n+1} y^n / n
// Converges when v₂(y) ≥ 1 (y even). Terms have valuation n·v₂(y) − v₂(n),
// vanishing when this reaches the bit width of T.
// For v₂(y) = 1, needs ~bits + log₂(bits) terms (135 for 128-bit).
// Uses a 2× budget to cover the v₂(y) = 1 case.
// NOTE: This computes log(1+y), NOT log(y). For log(y), pass y-1.
template<typename T>
constexpr T p_adic_log_1plus_impl(T y) noexcept {
    if (y == T(0)) return T(0);
    assert(p_adic_val(y) >= 1 && "p_adic_log_1plus: v2(y) must be >= 1");
    int bits = 8 * int(sizeof(T));
    int max_terms = 2 * bits;  // generous budget for v₂=1 inputs
    T sum = y;
    T term = y;
    for (int n = 2; n <= max_terms; ++n) {
        term = term * y;
        term = term * T(static_cast<uint64_t>(n - 1));
        term = div_by_int_generic(term, n);
        if (term == T(0)) break;
        sum = (n & 1) ? (sum + term) : (sum - term);
    }
    return sum;
}

// ℤ₂ exponential: exp(x) = Σ_{n=0}^{∞} x^n / n!
// Converges when v₂(x) ≥ 2 (x ≡ 0 mod 4). Terms have valuation
// n·(v₂(x)−1) + s₂(n), vanishing when this reaches the bit width of T.
// For v₂(x) = 1, the series never converges at finite precision
// (s₂(n) ≤ log₂(n)+1). Uses a 2× budget and early-exits on vanishing terms.
template<typename T>
constexpr T p_adic_exp_impl(T x) noexcept {
    int bits = 8 * int(sizeof(T));
    int v = p_adic_val(x);
    assert(v >= 2 && "p_adic_exp: v2(x) must be >= 2");
    int max_terms = (v < 2) ? (2 * bits) : ((bits + v - 3) / (v - 1) + 1);
    if (max_terms > 2 * bits) max_terms = 2 * bits;
    T sum = T(1);
    T term = T(1);
    for (int n = 1; n <= max_terms; ++n) {
        term = term * x;
        term = div_by_int_generic(term, n);
        if (term == T(0)) break;
        sum = sum + term;
    }
    return sum;
}

} // namespace detail

// Witt vector logarithm: applies ℤ₂ log componentwise to ghost values.
// Defined when a[0] is odd (unit in the Witt ring).
// Returns b such that witt_exp(b) = a (when both are defined).
//
// p_adic_log_1plus_impl(y) computes log(1+y). We want log(ghost_j(a)),
// we pass ghost_j(a) - 1 so that log(1 + (ghost_j(a)-1)) = log(ghost_j(a)).
// Convergence requires ghost_j(a) ≡ 1 (mod 2), which holds when a[0] is odd
// (a[0]^{2^j} ≡ 1 mod 2 for odd a[0]).
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_log(const WittVector<N, W>& a) noexcept {
    static_assert(N <= 32, "witt_log: N > 32 causes large stack allocation in ghost_recover");
    assert((a[0] & W(1)) != W(0) && "witt_log: requires a[0] odd (unit in Witt ring)");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> H{};
    for (int j = 0; j < N; ++j) {
        H[j] = detail::p_adic_log_1plus_impl(pows.ghost(j) - gw_t(1));
    }
    return detail::ghost_recover<N, W>(H);
}

// Witt vector exponential: applies ℤ₂ exp componentwise to ghost values.
// Requires a[0] ≡ 0 (mod 4) for convergence in the 2-adic series.
// Returns b such that witt_log(b) = a (when both are defined).
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_exp(const WittVector<N, W>& a) noexcept {
    static_assert(N <= 32, "witt_exp: N > 32 causes large stack allocation in ghost_recover");
    assert((a[0] & W(3)) == W(0) && "witt_exp: requires a[0] ≡ 0 (mod 4), v₂(a[0]) ≥ 2");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> H{};
    for (int j = 0; j < N; ++j) {
        H[j] = detail::p_adic_exp_impl(pows.ghost(j));
    }
    return detail::ghost_recover<N, W>(H);
}

// Witt vector inverse: a^{-1} in the Witt ring.
// Defined when a[0] is odd (unit). Uses componentwise ghost inversion.
// Satisfies a * witt_inverse(a) = τ(1) = {1, 0, ..., 0}.
template<int N, std::unsigned_integral W>
constexpr WittVector<N, W> witt_inverse(const WittVector<N, W>& a) noexcept {
    static_assert(N <= 32, "witt_inverse: N > 32 causes large stack allocation in ghost_recover");
    assert((a[0] & W(1)) != W(0) && "witt_inverse: requires a[0] odd (unit in Witt ring)");
    using gw_t = quad_width<W>;

    detail::GhostPowersFull<N, W> pows;
    pows.init(a.a.data());

    std::array<gw_t, N> H{};
    for (int j = 0; j < N; ++j) {
        gw_t Gj = pows.ghost(j);
        if constexpr (std::is_same_v<gw_t, detail::uint128_t>) {
            H[j] = detail::modinv_odd_128(Gj);
        } else {
            H[j] = dyadic::modinv_odd(Gj);
        }
    }
    return detail::ghost_recover<N, W>(H);
}

// ============================================================================
// 24. Verification Helpers for Witt Log/Exp/Inverse
// ============================================================================

template<int N, std::unsigned_integral W>
constexpr bool check_witt_inverse(const WittVector<N, W>& a) noexcept {
    auto inv = witt_inverse(a);
    auto prod = a * inv;
    for (int i = 0; i < N; ++i) {
        if (i == 0 && prod[i] != 1) return false;
        if (i > 0 && prod[i] != 0) return false;
    }
    return true;
}

// Check exp∘log = identity on the multiplicative domain (odd a[0]).
// For the additive roundtrip (log∘exp = id on v₂(a₀) ≥ 2), see PROOF_WITT_LOG_EXP.
template<int N, std::unsigned_integral W>
constexpr bool check_witt_exp_log_roundtrip(const WittVector<N, W>& a) noexcept {
    auto b = witt_exp(witt_log(a));
    for (int i = 0; i < N; ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

} // namespace dyadic
