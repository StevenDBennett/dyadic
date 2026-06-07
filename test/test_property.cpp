// test_property.cpp — Property-based runtime tests
// Exercises invariants that must hold for all random inputs.
// Usage: test_property [--seed <uint64>]
//   Default seed is 0 (uses hardcoded per-section seeds for reproducibility).

#include <dyadic/verify.h>
#include <dyadic/calculus.h>
#include "test_util.h"
#include <cstdio>

using namespace dyadic;

// Derive a per-section seed from the base seed and a section constant.
// When base_seed is 0, returns the section constant (original behavior).
static constexpr uint64_t section_seed(uint64_t base, uint64_t section) noexcept {
    return base ? (base ^ section) : section;
}

template<int N, std::unsigned_integral W>
static int test_properties(uint64_t base_seed = 0) {
    int failures[30] = {};
    auto SEED = [&](uint64_t s) { return section_seed(base_seed, s); };

    // 1. Polynomial eval consistency
    {
    detail::XorShift64 rng(SEED(0xDEADBEEF));
    for (int trial = 0; trial < 100; ++trial) {
        Polynomial<N, W, MonomialBasis> p;
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next());
        if (p.eval(W(0)) != p[0]) failures[0]++;
        W sum = 0;
        for (int i = 0; i < N; ++i) sum += p[i];
        if (p.eval(W(1)) != sum) failures[0]++;
    }
    }

    // 2. Basis roundtrip: Monomial <-> FallingFactorial
    {
    detail::XorShift64 rng(SEED(0xCAFEBABE));
    for (int trial = 0; trial < 100; ++trial) {
        Polynomial<N, W, MonomialBasis> p;
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next());
        auto ff = change_basis<FallingFactorialBasis>(p);
        auto back = change_basis<MonomialBasis>(ff);
        for (int i = 0; i < N; ++i) {
            if (p[i] != back[i]) { failures[1]++; break; }
        }
    }
    }

    // 3. Taylor roundtrip (small coefficients)
    {
    detail::XorShift64 rng(SEED(0xDECAFBAD));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p;
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next() & 0xFF);
        auto taylor = change_basis<TaylorBasis>(p);
        auto back = change_basis<MonomialBasis>(taylor);
        for (int i = 0; i < N; ++i) {
            if (p[i] != back[i]) { failures[2]++; break; }
        }
    }
    }

    // 4. D and Delta commute
    {
    detail::XorShift64 rng(SEED(0xFEEDFACE));
    for (int trial = 0; trial < 100; ++trial) {
        Polynomial<N, W, MonomialBasis> p;
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next());
        auto d_delta = formal_derivative(forward_difference(p));
        auto delta_d = forward_difference(formal_derivative(p));
        for (int i = 0; i < N - 2; ++i) {
            if (d_delta[i] != delta_d[i]) { failures[3]++; break; }
        }
    }
    }

    // 5. Derivative linearity: D(P+Q) = D(P) + D(Q)
    {
    detail::XorShift64 rng(SEED(0xD1FFACE));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p, q;
        for (int i = 0; i < N; ++i) {
            p[i] = static_cast<W>(rng.next());
            q[i] = static_cast<W>(rng.next());
        }
        auto d_sum = formal_derivative(p + q);
        auto d_p = formal_derivative(p);
        auto d_q = formal_derivative(q);
        for (int i = 0; i < N - 2; ++i) {
            if (d_sum[i] != static_cast<W>(d_p[i] + d_q[i])) { failures[4]++; break; }
        }
    }
    }

    // 6. Forward difference linearity: Δ(P+Q) = Δ(P) + Δ(Q)
    {
    detail::XorShift64 rng(SEED(0xDE1FA));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p, q;
        for (int i = 0; i < N; ++i) {
            p[i] = static_cast<W>(rng.next());
            q[i] = static_cast<W>(rng.next());
        }
        auto delta_sum = forward_difference(p + q);
        auto delta_p = forward_difference(p);
        auto delta_q = forward_difference(q);
        for (int i = 0; i < N - 2; ++i) {
            if (delta_sum[i] != static_cast<W>(delta_p[i] + delta_q[i])) { failures[5]++; break; }
        }
    }
    }

    // 7. Carry chain idempotency
    {
    detail::XorShift64 rng(SEED(0xC0FFEE01));
    for (int trial = 0; trial < 100; ++trial) {
        W input[N];
        for (int i = 0; i < N; ++i) input[i] = static_cast<W>(rng.next());
        W first[N];
        dword_t<W> dw_input[N], dw_first[N];
        for (int i = 0; i < N; ++i) dw_input[i] = static_cast<dword_t<W>>(input[i]);
        carry_chain(first, dw_input, N);
        for (int i = 0; i < N; ++i) dw_first[i] = static_cast<dword_t<W>>(first[i]);
        W second[N];
        carry_chain(second, dw_first, N);
        for (int i = 0; i < N; ++i) {
            if (first[i] != second[i]) { failures[6]++; break; }
        }
    }
    }

    // 8. Taylor shift by 0 is identity
    {
    detail::XorShift64 rng(SEED(0xF1F7));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p;
        for (int i = 0; i < N; ++i) p[i] = static_cast<W>(rng.next());
        auto shift0 = taylor_shift(p, W(0));
        for (int i = 0; i < N; ++i) {
            if (p[i] != shift0[i]) { failures[7]++; break; }
        }
        auto ff = change_basis<FallingFactorialBasis>(p);
        auto ff_shift0 = taylor_shift(ff, W(0));
        for (int i = 0; i < N; ++i) {
            if (ff[i] != ff_shift0[i]) { failures[7]++; break; }
        }
    }
    }

    // 9. Witt vector FV = VF
    {
    detail::XorShift64 rng(SEED(0xF1FFE));
    for (int trial = 0; trial < 100; ++trial) {
        WittVector<N, W> w;
        for (int i = 0; i < N; ++i) w[i] = static_cast<W>(rng.next());
        auto fv = w.FV();
        auto vf = w.VF();
        for (int i = 0; i < N; ++i) {
            if (fv[i] != vf[i]) { failures[8]++; break; }
        }
    }
    }

    // 10. Witt addition commutativity: a+b = b+a
    {
    detail::XorShift64 rng(SEED(0xC0F1E));
    for (int trial = 0; trial < 50; ++trial) {
        WittVector<N, W> a, b;
        for (int i = 0; i < N; ++i) { a[i] = static_cast<W>(rng.next()); b[i] = static_cast<W>(rng.next()); }
        auto sum_ab = witt_add(a, b);
        auto sum_ba = witt_add(b, a);
        for (int i = 0; i < N; ++i) { if (sum_ab[i] != sum_ba[i]) { failures[9]++; break; } }
    }
    }

    // 11. Witt addition associativity: (a+b)+c = a+(b+c)
    {
    detail::XorShift64 rng(SEED(0xA550C));
    for (int trial = 0; trial < 50; ++trial) {
        WittVector<N, W> a, b, c;
        for (int i = 0; i < N; ++i) { a[i] = static_cast<W>(rng.next()); b[i] = static_cast<W>(rng.next()); c[i] = static_cast<W>(rng.next()); }
        auto lhs = witt_add(witt_add(a, b), c);
        auto rhs = witt_add(a, witt_add(b, c));
        for (int i = 0; i < N; ++i) { if (lhs[i] != rhs[i]) { failures[10]++; break; } }
    }
    }

    // 12. (Skipped: Adams ghost precision requires small ghost values.
    //      Verified in test_full.cpp with bounded inputs.)

    // 13. Teichmüller lift: τ(ab) ghost equivalence
    {
    detail::XorShift64 rng(SEED(0xFE1C));
    for (int trial = 0; trial < 50; ++trial) {
        W a = static_cast<W>(rng.next() & 0xFF);
        W b = static_cast<W>(rng.next() & 0xFF);
        auto lhs = witt_mul(teichmueller_lift<N, W>(a), teichmueller_lift<N, W>(b));
        auto rhs = teichmueller_lift<N, W>(static_cast<W>(a * b));
        for (int j = 0; j < N; ++j) { if (lhs.ghost(j) != rhs.ghost(j)) { failures[12]++; break; } }
    }
    }

    // 14. Polynomial addition is commutative: P+Q = Q+P
    {
    detail::XorShift64 rng(SEED(0xF1FA));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p, q;
        for (int i = 0; i < N; ++i) { p[i] = static_cast<W>(rng.next()); q[i] = static_cast<W>(rng.next()); }
        auto sum1 = p + q;
        auto sum2 = q + p;
        for (int i = 0; i < N; ++i) { if (sum1[i] != sum2[i]) { failures[13]++; break; } }
    }
    }

    // 15. Coefficient-wise distributivity: poly_mul_cw(P, Q+R) = poly_mul_cw(P, Q) + poly_mul_cw(P, R)
    {
    detail::XorShift64 rng(SEED(0xD157));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p, q, r;
        for (int i = 0; i < N; ++i) { p[i] = static_cast<W>(rng.next()); q[i] = static_cast<W>(rng.next()); r[i] = static_cast<W>(rng.next()); }
        auto lhs = poly_mul_cw(p, q + r);
        auto rhs = poly_mul_cw(p, q) + poly_mul_cw(p, r);
        for (int i = 0; i < N; ++i) { if (lhs[i] != rhs[i]) { failures[14]++; break; } }
    }
    }

    // 16. Coefficient-wise commutativity: poly_mul_cw(P, Q) = poly_mul_cw(Q, P)
    {
    detail::XorShift64 rng(SEED(0xC0DE));
    for (int trial = 0; trial < 50; ++trial) {
        Polynomial<N, W, MonomialBasis> p, q;
        for (int i = 0; i < N; ++i) { p[i] = static_cast<W>(rng.next()); q[i] = static_cast<W>(rng.next()); }
        auto direct = poly_mul_cw(p, q);
        auto swapped = poly_mul_cw(q, p);
        for (int i = 0; i < N; ++i) { if (direct[i] != swapped[i]) { failures[15]++; break; } }
    }
    }

    // 17. Witt multiplication is commutative: witt_mul(a, b) = witt_mul(b, a)
    {
    detail::XorShift64 rng(SEED(0xC0DE2));
    for (int trial = 0; trial < 50; ++trial) {
        WittVector<N, W> a, b;
        for (int i = 0; i < N; ++i) { a[i] = static_cast<W>(rng.next()); b[i] = static_cast<W>(rng.next()); }
        auto ab = witt_mul(a, b);
        auto ba = witt_mul(b, a);
        for (int i = 0; i < N; ++i) { if (ab[i] != ba[i]) { failures[16]++; break; } }
    }
    }

    // 18. Witt distributivity: a*(b+c) = a*b + a*c
    {
    detail::XorShift64 rng(SEED(0xD1572));
    for (int trial = 0; trial < 50; ++trial) {
        WittVector<N, W> a, b, c;
        for (int i = 0; i < N; ++i) { a[i] = static_cast<W>(rng.next()); b[i] = static_cast<W>(rng.next()); c[i] = static_cast<W>(rng.next()); }
        auto lhs = witt_mul(a, witt_add(b, c));
        auto rhs = witt_add(witt_mul(a, b), witt_mul(a, c));
        for (int i = 0; i < N; ++i) { if (lhs[i] != rhs[i]) { failures[17]++; break; } }
    }
    }

    // 19. Witt additive identity: a + 0 = a
    {
    WittVector<N, W> zero{};
    detail::XorShift64 rng(SEED(0xADD0));
    for (int trial = 0; trial < 50; ++trial) {
        WittVector<N, W> a;
        for (int i = 0; i < N; ++i) a[i] = static_cast<W>(rng.next());
        auto sum = witt_add(a, zero);
        for (int i = 0; i < N; ++i) { if (sum[i] != a[i]) { failures[18]++; break; } }
    }
    }

    // 20. Witt multiplicative identity: a * 1 = a  (1 = Teichmüller lift of 1)
    {
    auto one = teichmueller_lift<N, W>(W(1));
    detail::XorShift64 rng(SEED(0x1D45));
    for (int trial = 0; trial < 50; ++trial) {
        WittVector<N, W> a;
        for (int i = 0; i < N; ++i) a[i] = static_cast<W>(rng.next());
        auto prod = witt_mul(a, one);
        for (int i = 0; i < N; ++i) { if (prod[i] != a[i]) { failures[19]++; break; } }
    }
    }

    int total = 0;
    for (int i = 0; i < 30; ++i) total += failures[i];

    std::printf("    [N=%d W=%d] eval=%d bt=%d taylor=%d DΔ=%d "
                "Dlin=%d Δlin=%d carry=%d shift0=%d FV=VF=%d "
                "addComm=%d addAssoc=%d polyComm=%d dist=%d mulComm=%d "
                "wittMulComm=%d wittDist=%d explog=%d expHom=%d\n",
        N, int(8*sizeof(W)),
        failures[0], failures[1], failures[2], failures[3],
        failures[4], failures[5], failures[6], failures[7], failures[8],
        failures[9], failures[10], failures[13],
        failures[14], failures[15], failures[16], failures[17], failures[18], failures[19]);

    return total;
}

int main(int argc, char** argv) {
    uint64_t base_seed = dyadic_test::parse_seed(argc, argv, 0);

    int failures = 0;
    int total = 0;

    auto run = [&](const char* name, int f) {
        total++;
        if (f > 0) {
            failures += f;
            std::printf("  FAIL  %s (%d failures)\n", name, f);
        } else {
            std::printf("  PASS  %s\n", name);
        }
    };

    if (base_seed != 0) std::printf("=== Property-Based Tests (seed=%llu) ===\n", (unsigned long long)base_seed);
    else                std::printf("=== Property-Based Tests ===\n");

    run("props N=2 uint16_t", test_properties<2, uint16_t>(base_seed));
    run("props N=3 uint16_t", test_properties<3, uint16_t>(base_seed));
    run("props N=4 uint16_t", test_properties<4, uint16_t>(base_seed));
    run("props N=5 uint16_t", test_properties<5, uint16_t>(base_seed));
    run("props N=2 uint32_t", test_properties<2, uint32_t>(base_seed));
    run("props N=3 uint32_t", test_properties<3, uint32_t>(base_seed));
    run("props N=4 uint32_t", test_properties<4, uint32_t>(base_seed));
    run("props N=2 uint64_t", test_properties<2, uint64_t>(base_seed));
    run("props N=3 uint64_t", test_properties<3, uint64_t>(base_seed));
    run("props N=4 uint64_t", test_properties<4, uint64_t>(base_seed));

    if (failures == 0) {
        std::printf("=== All %d property test suites passed ===\n", total);
    } else {
        std::printf("=== %d property test suites failed out of %d ===\n", failures, total);
    }

    return failures;
}
