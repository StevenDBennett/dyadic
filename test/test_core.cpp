// Test suite for dyadic-core — covers all six axioms and core machinery.
// Each test prints PASS/FAIL and returns non-zero on failure.

#include <dyadic/core.h>
#include <cstdio>
#include <cstdint>

using namespace dyadic;

static int failures = 0;

#define TEST(name, expr) do { \
    bool ok = !!(expr); \
    if (!ok) { std::printf("  FAIL: %s\n", name); ++failures; } \
    else { std::printf("  PASS: %s\n", name); } \
} while(0)

int main() {
    std::printf("=== dyadic-core test ===\n\n");

    // -----------------------------------------------------------------------
    // 2-adic primitives
    // -----------------------------------------------------------------------
    std::printf("--- 2-adic primitives ---\n");
    TEST("v2(0) == 32", v2(uint32_t(0)) == 32);
    TEST("v2(1) == 0",  v2(uint32_t(1)) == 0);
    TEST("v2(8) == 3",  v2(uint32_t(8)) == 3);
    TEST("v2(6) == 1",  v2(uint32_t(6)) == 1);

    TEST("modinv_odd(1) == 1",  modinv_odd(uint32_t(1)) == 1);
    TEST("modinv_odd(3) * 3 == 1", modinv_odd(uint32_t(3)) * uint32_t(3) == 1);
    TEST("modinv_odd(7) * 7 == 1", modinv_odd(uint32_t(7)) * uint32_t(7) == 1);
    TEST("modinv_odd(9) * 9 == 1", modinv_odd(uint32_t(9)) * uint32_t(9) == 1);

    TEST("exact_divide(6, 3) == 2", exact_divide(uint32_t(6), uint32_t(3)) == 2);
    TEST("exact_divide(42, 7) == 6", exact_divide(uint32_t(42), uint32_t(7)) == 6);

    TEST("div_2k_adic(8, 2) == 2",  div_2k_adic(uint32_t(8), 2) == 2);
    TEST("div_2k_adic(0, 4) == 0",  div_2k_adic(uint32_t(0), 4) == 0);

    // Artin-Schreier: P(0) = 0, P(1) = 0, P(x) = P(1-x)
    TEST("artin_schreier(0) == 0",  artin_schreier(uint32_t(0)) == 0);
    TEST("artin_schreier(1) == 0",  artin_schreier(uint32_t(1)) == 0);
    TEST("artin_schreier(5) == artin_schreier(1-5)",
         artin_schreier(uint32_t(5)) == artin_schreier(uint32_t(uint32_t(1) - uint32_t(5))));

    // -----------------------------------------------------------------------
    // Polynomial construction, eval (monomial basis), actual_degree
    // -----------------------------------------------------------------------
    std::printf("\n--- Polynomial ---\n");
    using Poly4 = Polynomial<4, uint32_t, MonomialBasis>;

    Poly4 p{{1, 2, 3, 4}};  // 1 + 2x + 3x^2 + 4x^3
    TEST("degree of 1+2x+3x^2+4x^3 == 3", p.actual_degree() == 3);
    TEST("eval at 0 == 1",    p.eval(0) == 1);
    TEST("eval at 1 == 10",   p.eval(1) == 10);
    TEST("eval at 2 == 49",   p.eval(2) == 49);

    Poly4 zero{};
    TEST("zero poly degree == -1", zero.actual_degree() == -1);
    TEST("eval zero poly == 0", zero.eval(5) == 0);

    Poly4 q{{0, 0, 0, 5}};  // 5x^3
    TEST("degree of 5x^3 == 3", q.actual_degree() == 3);

    Poly4 r = p + q;
    TEST("p+q: eval at 2 == 89", r.eval(2) == 89);
    Poly4 s = p - q;
    TEST("p-q: eval at 2 == 9",  s.eval(2) == 9);

    // -----------------------------------------------------------------------
    // synthetic_divide (Axiom 1)
    // -----------------------------------------------------------------------
    std::printf("\n--- synthetic_divide (Axiom 1) ---\n");
    {
        using P64 = Polynomial<4, uint64_t, MonomialBasis>;
        P64 a{{UINT64_MAX, 1, 0, 0}};
        auto sd = synthetic_divide(a);
        // After synthetic divide, each word is split in half
        // The high halves cascade as carries into the next word.
        // For uint64_t, shift=32, mask=0xFFFFFFFF.
        // a[0] = 0xFFFF... -> low=0xFFFFFFFF, carry=0xFFFFFFFF00000000>>32=0xFFFFFFFF
        // a[1] = 1 + carry = 0x100000000 -> low=0, carry=1
        // a[2] = 0 + carry = 1 -> low=1, carry=0
        TEST("synthetic_divide splits words", sd[0] == 0xFFFFFFFF);
    }

    // -----------------------------------------------------------------------
    // carry_normalize
    // -----------------------------------------------------------------------
    std::printf("\n--- carry_normalize ---\n");
    {
        auto rn = carry_normalize(Ring<uint32_t, Unnormalized>{42});
        TEST("normalize unnormalized", rn.value == 42);
        auto rm = carry_normalize(Ring<uint32_t, Modular<uint32_t>>{7});
        TEST("normalize modular", rm.value == 7);
        auto re = carry_normalize(Ring<uint32_t, Exact>{99});
        TEST("normalize exact", re.value == 99);
    }

    // -----------------------------------------------------------------------
    // Formal Derivative (Axiom 2)
    // -----------------------------------------------------------------------
    std::printf("\n--- formal_derivative (Axiom 2) ---\n");
    {
        Poly4 p{{3, 1, 4, 1}};  // 3 + x + 4x^2 + x^3
        auto d = formal_derivative(p);  // 1 + 8x + 3x^2
        using Poly3 = Polynomial<3, uint32_t, MonomialBasis>;
        Poly3 expected{{1, 8, 3}};
        bool match = true;
        for (int i = 0; i < 3; ++i) if (d[i] != expected[i]) match = false;
        TEST("D(3 + x + 4x^2 + x^3) == 1 + 8x + 3x^2", match);

        // D^3 of cubic = 6 (constant, 3! * leading coeff)
        auto d2d = formal_derivative(d);
        auto d3d = formal_derivative(d2d);
        TEST("D^3 of cubic is constant 6", d3d[0] == 6 && d3d.actual_degree() == 0);

        // D^4 of cubic is zero
        auto d4d = formal_derivative(d3d);
        TEST("D^4 of cubic is zero", d4d.actual_degree() == -1);
    }

    // -----------------------------------------------------------------------
    // Forward Difference (Axiom 3)
    // -----------------------------------------------------------------------
    std::printf("\n--- forward_difference (Axiom 3) ---\n");
    {
        // Δx^3 = 3x^2 + 3x + 1  (in monomial basis: 1 + 3x + 3x^2)
        Polynomial<4, uint32_t, MonomialBasis> p{{0, 0, 0, 1}};
        auto d = forward_difference(p);
        Polynomial<3, uint32_t, MonomialBasis> expected{{1, 3, 3}};
        bool match = true;
        for (int i = 0; i < 3; ++i) if (d[i] != expected[i]) match = false;
        TEST("Δ(x^3) == 1 + 3x + 3x^2", match);

        // Δ^3 of degree-3 polynomial is 6 (constant)
        auto d2 = forward_difference(d);
        auto d3 = forward_difference(d2);
        TEST("Δ^3(x^3) == 6 (constant)", d3[0] == 6 && d3.actual_degree() == 0);
    }

    // -----------------------------------------------------------------------
    // Exact variants (Axioms 4 & 5) — identity with D/Δ
    // -----------------------------------------------------------------------
    std::printf("\n--- exact variants (Axioms 4 & 5) ---\n");
    {
        Poly4 p{{3, 1, 4, 1}};
        auto d1 = formal_derivative(p);
        auto d2 = exact_formal_derivative(p);
        bool match_d = true;
        for (int i = 0; i < 3; ++i) if (d1[i] != d2[i]) match_d = false;
        TEST("exact_D == D", match_d);

        auto f1 = forward_difference(p);
        auto f2 = exact_forward_difference(p);
        bool match_f = true;
        for (int i = 0; i < 3; ++i) if (f1[i] != f2[i]) match_f = false;
        TEST("exact_Δ == Δ", match_f);
    }

    // -----------------------------------------------------------------------
    // Carry chain — Accum is dword_t<W> (wider type)
    // -----------------------------------------------------------------------
    std::printf("\n--- carry chain ---\n");
    {
        using dw = dword_t<uint32_t>;  // uint64_t
        uint32_t r[4];
        dw v[4] = {0xFFFFFFFF, 0xFFFFFFFF, 0x00000001, 0x00000000};
        auto carry = carry_chain<uint32_t>(r, v, 4);
        TEST("carry_chain carry=0", carry == 0);
        TEST("carry_chain[2]=1",    r[2] == 1);

        dw v2[3] = {0xFFFFFFFF, 0x00000001, 0x00000000};
        auto c2 = carry_chain<uint32_t>(r, v2, 3);
        TEST("carry_chain overflow[0]=0xFFFFFFFF", r[0] == 0xFFFFFFFF);
        TEST("carry_chain overflow[1]=1",         r[1] == 1);
    }

    // -----------------------------------------------------------------------
    // Polynomial multiplication (carry-chain)
    // -----------------------------------------------------------------------
    std::printf("\n--- polynomial multiplication ---\n");
    {
        Polynomial<3, uint32_t, MonomialBasis> a{{1, 1, 0}};  // 1 + x
        Polynomial<3, uint32_t, MonomialBasis> b{{1, 2, 0}};  // 1 + 2x
        auto c = a * b;  // (1+x)(1+2x) = 1 + 3x + 2x^2  (carry-chain ring)
        Polynomial<5, uint32_t, MonomialBasis> expected{{1, 3, 2, 0, 0}};
        bool match = true;
        for (int i = 0; i < 5; ++i) if (c[i] != expected[i]) match = false;
        TEST("(1+x)*(1+2x) carry-chain == 1+3x+2x^2", match);

        // Coefficient-wise version matches (no carries needed here)
        auto cw = poly_mul_cw(a, b);
        bool match_cw = true;
        for (int i = 0; i < 5; ++i) if (c[i] != cw[i]) match_cw = false;
        TEST("poly_mul_cw matches operator* for small polys", match_cw);

        // Verify identity: (1+x)^2 = 1 + 2x + x^2 (carry-chain)
        auto sq = a * a;
        Polynomial<5, uint32_t, MonomialBasis> sq_exp{{1, 2, 1, 0, 0}};
        bool sq_match = true;
        for (int i = 0; i < 5; ++i) if (sq[i] != sq_exp[i]) sq_match = false;
        TEST("(1+x)^2 == 1+2x+x^2 (carry-chain)", sq_match);
    }

    // -----------------------------------------------------------------------
    // verify_divmod_cw
    // -----------------------------------------------------------------------
    std::printf("\n--- verify_divmod_cw ---\n");
    {
        // A = (x+1)(x+2) = x^2 + 3x + 2, stored as [2,3,1,0] (pad to N=4)
        Polynomial<4, uint32_t, MonomialBasis> A{{2, 3, 1, 0}};
        // B = x + 1, stored as [1,1]
        Polynomial<2, uint32_t, MonomialBasis> B{{1, 1}};
        // Q = x + 2, pad to N=4
        Polynomial<4, uint32_t, MonomialBasis> Q{{2, 1, 0, 0}};
        Polynomial<2, uint32_t, MonomialBasis> R{};

        TEST("verify A = (x+2)*(x+1) + 0", verify_divmod_cw(A, Q, B, R));
    }

    // -----------------------------------------------------------------------
    // double-width type basics
    // -----------------------------------------------------------------------
    std::printf("\n--- double-width types ---\n");
    {
        // clang-format off
        TEST("dword_u8_is_u16", (std::is_same_v<dword_t<uint8_t>, uint16_t>));
        TEST("dword_u32_is_u64", (std::is_same_v<dword_t<uint32_t>, uint64_t>));
        { using dw64 = dword_t<uint64_t>;
          TEST("dword_u64_is_uint128", (std::is_same_v<dw64, detail::uint128_t>)); }
        // clang-format on

        detail::uint128_t x{uint64_t(0xAAAAAAAAAAAAAAAA), uint64_t(0xBBBBBBBBBBBBBBBB)};
        detail::uint128_t y{0, 1};
        auto sum = x + y;
        TEST("uint128_t add lo", sum.lo == uint64_t(0xBBBBBBBBBBBBBBBC));
        TEST("uint128_t add hi", sum.hi == uint64_t(0xAAAAAAAAAAAAAAAA));

        auto z = detail::uint128_t(3);
        auto m = detail::uint128_t(7);
        TEST("uint128_t mul", (z * m).lo == 21);
    }

    std::printf("\n=== %s ===\n", failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
    return failures;
}
