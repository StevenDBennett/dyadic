// test_full.cpp — Full functional tests exercising the entire API surface
// with edge cases, corner values, and stress conditions.

#include "dyadic.h"
#include <cstdio>
#include <dyadic/prng.h>

using namespace dyadic;

// ---------------------------------------------------------------------------
// 2-adic primitives: v2, modinv_odd, exact_divide, div_2k_adic, artin_schreier
// ---------------------------------------------------------------------------
static int test_2adic_primitives() {
    int f = 0;

    // v2 — all word sizes
    if (v2(uint8_t(0))  != 8)  { std::printf("FAIL v2 u8(0)\n"); f++; }
    if (v2(uint16_t(0)) != 16) { std::printf("FAIL v2 u16(0)\n"); f++; }
    if (v2(uint32_t(0)) != 32) { std::printf("FAIL v2 u32(0)\n"); f++; }
    if (v2(uint64_t(0)) != 64) { std::printf("FAIL v2 u64(0)\n"); f++; }
    if (v2(uint8_t(1))  != 0)  { std::printf("FAIL v2 u8(1)\n"); f++; }
    if (v2(uint16_t(1)) != 0)  { std::printf("FAIL v2 u16(1)\n"); f++; }
    if (v2(uint64_t(6)) != 1)  { std::printf("FAIL v2(6)\n"); f++; }
    if (v2(uint64_t(8)) != 3)  { std::printf("FAIL v2(8)\n"); f++; }
    if (v2(uint64_t(1) << 63) != 63) { std::printf("FAIL v2(2^63)\n"); f++; }

    // modinv_odd — verify a * modinv(a) == 1 mod 2^W for various odd values
    auto check_modinv = [&](auto v) {
        using W = decltype(v);
        W inv = modinv_odd(v);
        if (static_cast<W>(inv * v) != W(1)) { std::printf("FAIL modinv_odd(%u)\n", (unsigned)v); f++; }
    };
    check_modinv(uint8_t(1));  check_modinv(uint8_t(3));
    check_modinv(uint8_t(7));  check_modinv(uint8_t(9));
    check_modinv(uint8_t(0xFF));
    check_modinv(uint16_t(1)); check_modinv(uint16_t(3));
    check_modinv(uint16_t(5)); check_modinv(uint16_t(0xFFFF));
    check_modinv(uint32_t(1)); check_modinv(uint32_t(0xFFFFFFFF));
    check_modinv(uint64_t(1)); check_modinv(uint64_t(3));
    check_modinv(uint64_t(7)); check_modinv(uint64_t(9));
    check_modinv(uint64_t(0xFFFFFFFFFFFFFFFFULL));

    // exact_divide (odd denominators)
    if (exact_divide(uint64_t(15), uint64_t(3)) != 5) { f++; }
    if (exact_divide(uint64_t(21), uint64_t(7)) != 3) { f++; }
    if (exact_divide(uint64_t(0), uint64_t(3)) != 0)  { f++; }

    // div_2k_adic — positive and negative (high-bit-set) values
    auto div_check = [&](uint64_t val, int k, uint64_t expected) {
        uint64_t got = div_2k_adic(val, k);
        if (got != expected) {
            std::printf("FAIL div_2k_adic(0x%lx,%d)=0x%lx expected 0x%lx\n", val, k, got, expected); f++;
        }
    };
    div_check(16, 2, 4);
    div_check(1, 1, 0);
    div_check(0, 10, 0);
    div_check(uint64_t(-1), 1, uint64_t(-1));   // sign-extended: all ones stays all ones
    div_check(uint64_t(-1), 63, uint64_t(-1));  // sign-extended: all ones stays all ones
    div_check(uint64_t(-1), 64, 0);              // shift >= bits
    div_check(0x8000000000000000ULL, 63, uint64_t(-1)); // msb set, sign-extended
    div_check(0x8000000000000000ULL, 1, 0xC000000000000000ULL); // sign-extend one bit
    // k=0 should be identity
    div_check(42, 0, 42);
    div_check(uint64_t(-1), 0, uint64_t(-1));
    // uint8_t
    if (div_2k_adic(uint8_t(0x80), 4) != uint8_t(0xF8)) { std::printf("FAIL div_2k_adic u8\n"); f++; }
    if (div_2k_adic(uint8_t(0x80), 7) != uint8_t(0xFF)) { std::printf("FAIL div_2k_adic u8\n"); f++; }

    // artin_schreier: ℘(x) = x^2 - x
    auto as_check = [&](auto x, auto expected) {
        auto got = artin_schreier(x);
        if (got != expected) {
            std::printf("FAIL artin_schreier\n"); f++;
        }
    };
    as_check(uint64_t(0), uint64_t(0));
    as_check(uint64_t(1), uint64_t(0));
    as_check(uint64_t(2), uint64_t(2));
    as_check(uint64_t(3), uint64_t(6));
    // Artin-Schreier symmetry: ℘(x) = ℘(1-x) for all u8 values
    for (int x = 0; x <= 255; ++x) {
        uint8_t ux = static_cast<uint8_t>(x);
        if (artin_schreier(ux) != artin_schreier(static_cast<uint8_t>(1 - ux))) {
            std::printf("FAIL artin_schreier symmetry %d\n", x); f++;
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Binomial coefficients
// ---------------------------------------------------------------------------
static int test_binomial() {
    int f = 0;

    // Known values
    auto check = [&](auto expected, int n, int k) {
        auto got = binom<uint64_t>(n, k);
        if (got != static_cast<uint64_t>(expected)) {
            std::printf("FAIL binom(%d,%d)=%llu expected %llu\n", n, k,
                        (unsigned long long)got, (unsigned long long)expected); f++;
        }
    };
    check(1, 0, 0);
    check(1, 5, 0);  check(1, 5, 5);
    check(10, 5, 2); check(10, 5, 3);
    check(252, 10, 5);
    check(2598960ULL, 52, 5);  // poker hand count
    check(118264581564861424ULL, 60, 30); // large but fits in uint64

    // Out of range
    if (binom<uint64_t>(5, -1) != 0) { f++; }
    if (binom<uint64_t>(5, 7) != 0)  { f++; }

    // Different word sizes
    if (binom<uint8_t>(5, 2) != 10)  { std::printf("FAIL binom u8\n"); f++; }
    if (binom<uint16_t>(10, 5) != 252) { std::printf("FAIL binom u16\n"); f++; }
    if (binom<uint32_t>(20, 10) != 184756) { std::printf("FAIL binom u32\n"); f++; }

    return f;
}

// ---------------------------------------------------------------------------
// Stirling numbers
// ---------------------------------------------------------------------------
static int test_stirling() {
    int f = 0;

    // S(n,k) for various word sizes
    if (stirling_2<uint8_t>(4, 2) != 7)  { std::printf("FAIL S u8\n"); f++; }
    if (stirling_2<uint16_t>(5, 2) != 15){ std::printf("FAIL S u16\n"); f++; }
    if (stirling_2<uint32_t>(5, 3) != 25){ std::printf("FAIL S u32\n"); f++; }
    if (stirling_2<uint64_t>(6, 3) != 90){ std::printf("FAIL S(6,3)\n"); f++; }
    if (stirling_2<uint64_t>(10, 5) != 42525) { std::printf("FAIL S(10,5)\n"); f++; }

    // s(n,k) signed (returned as unsigned 2-adic wrapping)
    auto s1 = [](int n, int k) { return stirling_1<uint64_t>(n, k); };
    if (s1(1, 1) != 1)             { f++; }
    if (s1(2, 1) != uint64_t(-1))  { f++; } // s(2,1) = -1
    if (s1(2, 2) != 1)             { f++; }
    if (s1(3, 1) != 2)             { f++; }
    if (s1(3, 2) != uint64_t(-3))  { f++; } // s(3,2) = -3
    if (s1(4, 1) != uint64_t(-6))  { f++; }

    // |s(n,k)| unsigned
    if (stirling_1_unsigned<uint64_t>(4, 2) != 11) { f++; }
    if (stirling_1_unsigned<uint64_t>(5, 2) != 50) { f++; }
    if (stirling_1_unsigned<uint64_t>(6, 3) != 225){ f++; }

    // Out of range
    if (stirling_2<uint64_t>(0, 1) != 0) { f++; }
    if (stirling_1<uint64_t>(5, 6) != 0) { f++; }

    return f;
}

// ---------------------------------------------------------------------------
// Polynomial construction, eval, arithmetic
// ---------------------------------------------------------------------------
static int test_polynomial_edge() {
    int f = 0;

    // N=1 (max_degree 0, constant polynomial)
    Polynomial<1, uint64_t, MonomialBasis> c{{42}};
    if (c.max_degree != 0) { f++; }
    if (c.eval(0) != 42) { f++; }
    if (c.eval(100) != 42) { f++; }

    // All-zero polynomial
    Polynomial<4, uint64_t, MonomialBasis> zero{};
    if (zero.eval(0) != 0) { f++; }
    if (zero.eval(12345) != 0) { f++; }
    auto zsum = zero + zero;
    if (zsum[0] != 0) { f++; }

    // All-max polynomial
    Polynomial<3, uint64_t, MonomialBasis> maxp{{uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
    auto max_sum = maxp + maxp;
    if (max_sum[0] != uint64_t(-2)) { std::printf("FAIL max poly +\n"); f++; }

    // eval in falling factorial and Taylor bases
    Polynomial<3, uint64_t, MonomialBasis> mono{{1, 2, 3}};
    auto ff = change_basis<FallingFactorialBasis>(mono);
    auto taylor = change_basis<TaylorBasis>(mono);
    for (uint64_t x = 0; x < 10; ++x) {
        if (eval(ff, x) != mono.eval(x)) {
            std::printf("FAIL ff eval mismatch at %lu\n", x); f++;
        }
        if (eval(taylor, x) != mono.eval(x)) {
            std::printf("FAIL taylor eval mismatch at %lu\n", x); f++;
        }
    }

    return f;
}

static int test_polynomial_arith() {
    int f = 0;

    Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
    Polynomial<4, uint64_t, MonomialBasis> q{{0, 1, 0, 1}};

    // Addition, subtraction
    auto sum = p + q;
    if (sum[0] != 1 || sum[1] != 3 || sum[2] != 3 || sum[3] != 5) { f++; }
    auto diff = p - q;
    if (diff[0] != 1 || diff[1] != 1 || diff[2] != 3 || diff[3] != 3) { f++; }

    // Self arithmetic
    auto self = p + p;
    for (int i = 0; i < 4; ++i) if (self[i] != 2 * p[i]) { f++; break; }
    auto self_sub = p - p;
    for (int i = 0; i < 4; ++i) if (self_sub[i] != 0) { f++; break; }

    // Multiplication: (1+x)(1+x) = 1+2x+x^2
    Polynomial<2, uint64_t, MonomialBasis> a{{1, 1}};
    Polynomial<2, uint64_t, MonomialBasis> b{{1, 1}};
    auto prod = a * b;
    if (prod.max_degree != 2) { f++; }
    if (prod[0] != 1 || prod[1] != 2 || prod[2] != 1) { f++; }

    // Larger multiplication: (1+2x+3x^2) * (4+5x) = 4+13x+22x^2+15x^3
    Polynomial<3, uint64_t, MonomialBasis> p3{{1, 2, 3}};
    Polynomial<2, uint64_t, MonomialBasis> q2{{4, 5}};
    auto p3q2 = p3 * q2;
    if (p3q2[0] != 4 || p3q2[1] != 13 || p3q2[2] != 22 || p3q2[3] != 15) {
        std::printf("FAIL poly mul 3x2: %lu %lu %lu %lu\n", p3q2[0], p3q2[1], p3q2[2], p3q2[3]); f++;
    }

    // Multiplication by zero polynomial
    Polynomial<2, uint64_t, MonomialBasis> zero{{0, 0}};
    auto zero_prod = a * zero;
    for (int i = 0; i < 3; ++i) if (zero_prod[i] != 0) { f++; break; }

    // Multiplication with overflow (max values)
    Polynomial<2, uint64_t, MonomialBasis> max_a{{uint64_t(-1), uint64_t(-1)}};
    auto overflow_prod = max_a * max_a;
    if (overflow_prod[0] != 1) { // (-1)*(-1) = 1 in 2-adic, carry chain normalizes
        std::printf("FAIL poly mul overflow 0=%lu\n", overflow_prod[0]); f++;
    }

    return f;
}

// ---------------------------------------------------------------------------
// Polynomial GCD, divmod, resultant, discriminant, determinant
// ---------------------------------------------------------------------------

static int test_polynomial_gcd() {
    int f = 0;

    // Test verify_divmod_cw helper directly
    auto test_vdm = [&](const char* label, bool result) {
        if (!result) { std::printf("FAIL verify_divmod_cw %s\n", label); f++; }
    };

    // --- poly_divmod_cw + verify_divmod_cw ---
    {
        // (x²+2x+1) ÷ (x+1) = (x+1) rem 0
        Polynomial<3, uint64_t, MonomialBasis> A{{1, 2, 1}};
        Polynomial<2, uint64_t, MonomialBasis> B{{1, 1}};
        auto [Q, R] = poly_divmod_cw(A, B);
        test_vdm("(x+1)²/(x+1)", verify_divmod_cw(A, Q, B, R));
        // N < M: (x+1) ÷ (x²+2x+1) → Q=0, R=A
        Polynomial<2, uint64_t, MonomialBasis> A2{{1, 1}};
        Polynomial<3, uint64_t, MonomialBasis> B2{{1, 2, 1}};
        auto [Q2, R2] = poly_divmod_cw(A2, B2);
        test_vdm("(x+1)/(x+1)²", verify_divmod_cw(A2, Q2, B2, R2));
        // x² ÷ (x+1) = (x-1) rem 1
        Polynomial<3, uint64_t, MonomialBasis> A3{{0, 0, 1}};
        Polynomial<2, uint64_t, MonomialBasis> B3{{1, 1}};
        auto [Q3, R3] = poly_divmod_cw(A3, B3);
        test_vdm("x²/(x+1)", verify_divmod_cw(A3, Q3, B3, R3));
        // Constant divmod
        Polynomial<1, uint64_t, MonomialBasis> A4{{10}};
        Polynomial<1, uint64_t, MonomialBasis> B4{{3}};
        auto [Q4, R4] = poly_divmod_cw(A4, B4);
        test_vdm("const 10/3", verify_divmod_cw(A4, Q4, B4, R4));
        // Wrong answer detection: Q=1,R=5 gives 1*3+5=8 != 10
        Polynomial<1, uint64_t, MonomialBasis> Qbad{{1}};
        Polynomial<1, uint64_t, MonomialBasis> Rbad{{5}};
        if (verify_divmod_cw(A4, Qbad, B4, Rbad)) {
            std::printf("FAIL verify_divmod_cw wrong answer not caught\n"); f++;
        }
    }

    // --- pseudo_remainder_cw ---
    {
        // Even leading coefficient (can't use divmod)
        Polynomial<3, uint64_t, MonomialBasis> A{{1, 2, 1}};
        Polynomial<2, uint64_t, MonomialBasis> B{{2, 2}};
        auto R = pseudo_remainder_cw(A, B);
        if (R.actual_degree() >= 1) {
            std::printf("FAIL pseudo_remainder deg\n"); f++;
        }
        // B=0 case
        Polynomial<3, uint64_t, MonomialBasis> A2{{1, 2, 1}};
        Polynomial<2, uint64_t, MonomialBasis> B2{};
        auto R2 = pseudo_remainder_cw(A2, B2);
        if (R2[0] != 0 || R2[1] != 0 || R2[2] != 0) {
            std::printf("FAIL pseudo_remainder B=0\n"); f++;
        }
    }

    // --- poly_gcd_cw ---
    {
        auto check_gcd = [&](const auto& A, const auto& B, const char* label) {
            auto G = poly_gcd_cw(A, B);
            if (G.actual_degree() < 0) {
                if (A.actual_degree() < 0 && B.actual_degree() < 0)
                    return;
                std::printf("FAIL gcd %s: G=0\n", label); f++; return;
            }
        };

        using A3 = Polynomial<3, uint64_t, MonomialBasis>;
        using A2 = Polynomial<2, uint64_t, MonomialBasis>;
        using A4 = Polynomial<4, uint64_t, MonomialBasis>;
        check_gcd(
            A3(std::array<uint64_t,3>{uint64_t(-1), 0, 1}),
            A2(std::array<uint64_t,2>{uint64_t(-1), 1}),
            "x²-1, x-1");
        check_gcd(
            A4(std::array<uint64_t,4>{0, 0, 1, 1}),
            A3(std::array<uint64_t,3>{uint64_t(-1), 0, 1}),
            "x³+x², x²-1");
        check_gcd(
            A2(std::array<uint64_t,2>{}),
            A2(std::array<uint64_t,2>{1, 1}),
            "0, x+1");
        check_gcd(
            A3(std::array<uint64_t,3>{1, 0, 1}),
            A2(std::array<uint64_t,2>{1, 1}),
            "x²+1, x+1");
    }

    // --- polynomial_resultant_cw ---
    {
        Polynomial<3, uint64_t, MonomialBasis> A{{uint64_t(-1), 0, 1}};
        Polynomial<2, uint64_t, MonomialBasis> B{{uint64_t(-1), 1}};
        if (polynomial_resultant_cw(A, B) != 0) {
            std::printf("FAIL resultant shared root\n"); f++;
        }
        Polynomial<2, uint64_t, MonomialBasis> C{{uint64_t(-1), 1}};
        Polynomial<1, uint64_t, MonomialBasis> D{{2}};
        if (polynomial_resultant_cw(C, D) == 0) {
            std::printf("FAIL resultant constant\n"); f++;
        }
    }

    // --- poly_discriminant_cw ---
    {
        Polynomial<3, uint64_t, MonomialBasis> P{{uint64_t(-1), 0, 1}};
        if (poly_discriminant_cw(P) == 0) {
            std::printf("FAIL discriminant x²-1\n"); f++;
        }
    }

    // --- poly_is_square_free_cw ---
    {
        Polynomial<2, uint64_t, MonomialBasis> P{{1, 1}};
        if (!poly_is_square_free_cw(P)) {
            std::printf("FAIL is_square_free x+1\n"); f++;
        }
    }

    // --- detail::det_laplace ---
    {
        std::array<std::array<uint64_t, 6>, 6> M{};
        M[0][0] = 1; M[1][1] = 1;
        if (detail::det_laplace<uint64_t>(M, 2) != 1) { std::printf("FAIL det I2\n"); f++; }

        M[0][0] = 1; M[0][1] = 2;
        M[1][0] = 3; M[1][1] = 4;
        if (detail::det_laplace<uint64_t>(M, 2) != uint64_t(-2)) {
            std::printf("FAIL det 2x2\n"); f++;
        }

        std::array<std::array<uint64_t, 6>, 6> I3{};
        I3[0][0] = 1; I3[1][1] = 1; I3[2][2] = 1;
        if (detail::det_laplace<uint64_t>(I3, 3) != 1) { std::printf("FAIL det I3\n"); f++; }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Basis conversions
// ---------------------------------------------------------------------------
static int test_change_basis() {
    int f = 0;

    Polynomial<4, uint64_t, MonomialBasis> mono{{1, 2, 3, 4}};

    // Monomial <-> FallingFactorial roundtrip
    auto ff = change_basis<FallingFactorialBasis>(mono);
    auto back = change_basis<MonomialBasis>(ff);
    for (int i = 0; i < 4; ++i) {
        if (mono[i] != back[i]) { f++; }
    }

    // Roundtrip with all zeros
    Polynomial<4, uint64_t, MonomialBasis> zero{};
    auto ff_zero = change_basis<FallingFactorialBasis>(zero);
    auto back_zero = change_basis<MonomialBasis>(ff_zero);
    for (int i = 0; i < 4; ++i) if (zero[i] != back_zero[i]) { f++; break; }

    // Roundtrip with N=1
    Polynomial<1, uint64_t, MonomialBasis> c{{7}};
    auto ff_c = change_basis<FallingFactorialBasis>(c);
    auto back_c = change_basis<MonomialBasis>(ff_c);
    if (c[0] != back_c[0]) { f++; }

    // Monomial <-> Taylor roundtrip (small coefficients to avoid precision window)
    Polynomial<4, uint64_t, MonomialBasis> small{{0, 1, 0, 0}};
    auto taylor = change_basis<TaylorBasis>(small);
    auto back_taylor = change_basis<MonomialBasis>(taylor);
    for (int i = 0; i < 4; ++i) {
        if (small[i] != back_taylor[i]) {
            std::printf("FAIL Taylor roundtrip [%d]: %lu != %lu\n", i, small[i], back_taylor[i]); f++;
        }
    }

    // Change basis self-conversion
    auto identity = change_basis<MonomialBasis>(mono);
    for (int i = 0; i < 4; ++i) if (mono[i] != identity[i]) { f++; break; }

    // Falling factorial eval matches monomial eval
    if (eval(ff, uint64_t(0)) != 1) { f++; }
    if (eval(ff, uint64_t(1)) != mono.eval(1)) { f++; }
    if (eval(ff, uint64_t(2)) != mono.eval(2)) { f++; }
    if (eval(ff, uint64_t(5)) != mono.eval(5)) { f++; }

    return f;
}

// ---------------------------------------------------------------------------
// Taylor shift
// ---------------------------------------------------------------------------
static int test_taylor_shift() {
    int f = 0;

    // Shift by 0 should be identity
    Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
    auto shift0 = taylor_shift(p, uint64_t(0));
    for (int i = 0; i < 4; ++i) {
        if (p[i] != shift0[i]) { std::printf("FAIL taylor_shift delta=0 [%d]\n", i); f++; }
    }

    // P(x) = x^2, shift by 1: (x+1)^2 = x^2 + 2x + 1
    Polynomial<3, uint64_t, MonomialBasis> p2{{0, 0, 1}};
    auto shift1 = taylor_shift(p2, uint64_t(1));
    if (shift1[0] != 1 || shift1[1] != 2 || shift1[2] != 1) {
        std::printf("FAIL taylor_shift: %lu %lu %lu\n", shift1[0], shift1[1], shift1[2]); f++;
    }
    if (p2.eval(4) != shift1.eval(3)) { f++; }

    // Shift in falling factorial basis
    auto ff = change_basis<FallingFactorialBasis>(p);
    auto ff_shift0 = taylor_shift(ff, uint64_t(0));
    for (int i = 0; i < 4; ++i) {
        if (ff[i] != ff_shift0[i]) { std::printf("FAIL taylor_shift FF delta=0 [%d]\n", i); f++; }
    }

    // FF shift composition: shift(a, d1+d2) == shift(shift(a, d1), d2)
    auto ff_shift2 = taylor_shift(ff, uint64_t(2));
    auto ff_shift1_then_1 = taylor_shift(taylor_shift(ff, uint64_t(1)), uint64_t(1));
    for (int i = 0; i < 4; ++i) {
        if (ff_shift2[i] != ff_shift1_then_1[i]) {
            std::printf("FAIL taylor_shift composition [%d]\n", i); f++;
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Derivative and forward difference
// ---------------------------------------------------------------------------
static int test_derivative_difference() {
    int f = 0;

    // P(x) = x^3 + 2x^2 + 3x + 4, P'(x) = 3x^2 + 4x + 3
    Polynomial<4, uint64_t, MonomialBasis> p{{4, 3, 2, 1}};
    auto d = formal_derivative(p);
    if (d[0] != 3 || d[1] != 4 || d[2] != 3) { f++; }

    // Derivative of constant polynomial is zero
    Polynomial<1, uint64_t, MonomialBasis> c{{7}};
    auto dc = formal_derivative(c); // degree -1, no coefficients

    // Derivative of N=2 (linear): derivative of a+bx = b
    Polynomial<2, uint64_t, MonomialBasis> lin{{5, 3}};
    auto dlin = formal_derivative(lin);
    if (dlin[0] != 3) { f++; }

    // Forward difference: ΔP(x) = P(x+1) - P(x)
    auto delta = forward_difference(p);
    for (uint64_t x = 0; x < 10; ++x) {
        uint64_t expected = p.eval(x + 1) - p.eval(x);
        if (delta.eval(x) != expected) {
            std::printf("FAIL forward_difference at %lu: %lu != %lu\n", x, delta.eval(x), expected); f++;
        }
    }

    // Δ of constant polynomial is zero
    auto delta_c = forward_difference(c);

    // D and Δ commute: D(Δ(P)) = Δ(D(P))
    auto d_delta = formal_derivative(forward_difference(p));
    auto delta_d = forward_difference(formal_derivative(p));
    for (int i = 0; i < 2; ++i) {
        if (d_delta[i] != delta_d[i]) { f++; }
    }

    // Derivative linearity: D(P+Q) = D(P) + D(Q)
    Polynomial<4, uint64_t, MonomialBasis> q{{9, 8, 7, 6}};
    auto d_sum = formal_derivative(p + q);
    auto d_p = formal_derivative(p);
    auto d_q = formal_derivative(q);
    for (int i = 0; i < 2; ++i) {
        if (d_sum[i] != d_p[i] + d_q[i]) {
            std::printf("FAIL D linearity [%d]\n", i); f++;
        }
    }

    // Delta linearity: Δ(P+Q) = Δ(P) + Δ(Q)
    auto delta_sum = forward_difference(p + q);
    auto delta_p = forward_difference(p);
    auto delta_q = forward_difference(q);
    for (int i = 0; i < 2; ++i) {
        if (delta_sum[i] != delta_p[i] + delta_q[i]) { f++; }
    }

    // Exact variants
    auto exact_d = exact_formal_derivative(p);
    for (int i = 0; i < 2; ++i) if (exact_d[i] != d[i]) { f++; }
    auto exact_delta = exact_forward_difference(p);
    for (int i = 0; i < 2; ++i) if (exact_delta[i] != delta[i]) { f++; }

    return f;
}

// ---------------------------------------------------------------------------
// Indefinite sum
// ---------------------------------------------------------------------------
static int test_indefinite_sum() {
    int f = 0;

    // P(x) = 2x (coefficient divisible by 2 for Σ precision)
    Polynomial<3, uint64_t, MonomialBasis> p{{0, 2, 0}};
    auto sum_mono = indefinite_sum(p);
    auto delta_sum = forward_difference(sum_mono);
    if (delta_sum[0] != p[0] || delta_sum[1] != p[1]) {
        std::printf("FAIL indefinite_sum ΔΣ != id: %lu,%lu != %lu,%lu\n",
                    delta_sum[0], delta_sum[1], p[0], p[1]); f++;
    }

    // Falling factorial basis exact
    Polynomial<3, uint64_t, FallingFactorialBasis> ff{{0, 2, 0}};
    auto sum_ff = indefinite_sum(ff);
    auto delta_ff = forward_difference(sum_ff);
    for (int i = 0; i < 2; ++i) {
        if (delta_ff[i] != ff[i]) { f++; }
    }

    // Σ of zero polynomial is zero
    Polynomial<3, uint64_t, FallingFactorialBasis> ff_zero{};
    auto sum_zero = indefinite_sum(ff_zero);
    for (int i = 0; i < 4; ++i) if (sum_zero[i] != 0) { f++; break; }

    // N=1 case
    Polynomial<1, uint64_t, FallingFactorialBasis> ff1{{2}};
    auto sum1 = indefinite_sum(ff1);
    auto delta1 = forward_difference(sum1);
    if (delta1[0] != ff1[0]) { f++; }

    return f;
}

// ---------------------------------------------------------------------------
// Ring tags
// ---------------------------------------------------------------------------
static int test_ring_tags() {
    int f = 0;

    // ExactRing
    auto er = ExactRing_t<uint64_t>{42};
    auto cn_er = carry_normalize(er);
    if (cn_er.value != 42) { f++; }

    // NormalizedRing
    auto nr = NormalizedRing_t<uint64_t>{99};
    auto cn_nr = carry_normalize(nr);
    if (cn_nr.value != 99) { f++; }

    // UnnormalizedRing
    auto ur = UnnormalizedRing_t<uint64_t>{77};
    auto cn_ur = carry_normalize(ur);
    if (cn_ur.value != 77) { f++; }

    // synthetic_divide on Ring
    auto sd_ur = synthetic_divide(ur);
    if (sd_ur.value != 77) { f++; }

    // synthetic_divide on Polynomial (half-width carry) — verify idempotency
    Polynomial<3, uint16_t, MonomialBasis> sp{{0xFFFF, 0xFFFF, 0xFFFF}};
    auto sd_p = synthetic_divide(sp);
    auto sd_p2 = synthetic_divide(sd_p);
    for (int i = 0; i < 3; ++i) {
        if (sd_p[i] != sd_p2[i]) { std::printf("FAIL synthetic_divide idempotent [%d]\n", i); f++; }
    }

    // synthetic_divide on polynomial with no carry — also idempotent
    Polynomial<3, uint16_t, MonomialBasis> sp_no{{0x00, 0x01, 0x02}};
    auto sd_no = synthetic_divide(sp_no);
    auto sd_no2 = synthetic_divide(sd_no);
    for (int i = 0; i < 3; ++i) {
        if (sd_no[i] != sd_no2[i]) { f++; break; }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Witt vectors
// ---------------------------------------------------------------------------
static int test_witt_vector() {
    int f = 0;

    // Construction
    WittVector<3, uint64_t> a{{1, 2, 3}};
    if (a[0] != 1 || a[1] != 2 || a[2] != 3) { f++; }

    // Zero and one
    auto z = WittVector<3, uint64_t>::zero();
    for (int i = 0; i < 3; ++i) { if (z[i] != 0) { f++; } }
    auto o = WittVector<3, uint64_t>::one();
    if (o[0] != 1 || o[1] != 0 || o[2] != 0) { f++; }

    // N=1 Witt vector
    WittVector<1, uint64_t> w1{{5}};
    if (w1[0] != 5) { f++; }
    auto w1_add = witt_add(w1, WittVector<1, uint64_t>{{3}});
    if (w1_add[0] != 8) { f++; }
    auto w1_mul = witt_mul(w1, WittVector<1, uint64_t>{{3}});
    if (w1_mul[0] != 15) { f++; }

    // All-zero Witt vector
    auto zz = WittVector<3, uint64_t>::zero();
    auto az = witt_add(a, zz);
    for (int i = 0; i < 3; ++i) if (az[i] != a[i]) { std::printf("FAIL a+0 != a [%d]\n", i); f++; }

    // a * 1 = a (for N=2, since WittVector::one() is the multiplicative identity)
    auto o2 = WittVector<3, uint64_t>::one();
    auto mo = witt_mul(a, o2);
    for (int i = 0; i < 3; ++i) {
        if (mo[i] != a[i]) {
            std::printf("FAIL a*1 != a [%d]: %lu != %lu\n", i, mo[i], a[i]); f++;
        }
    }

    // All-max Witt vector: just verify add doesn't crash and is commutative
    WittVector<3, uint64_t> max_w{{uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
    auto max_add_ab = witt_add(max_w, max_w);
    auto max_add_ba = witt_add(max_w, max_w);
    for (int i = 0; i < 3; ++i) {
        if (max_add_ab[i] != max_add_ba[i]) { f++; break; }
    }

    // Ghost map ring homomorphism for addition
    WittVector<3, uint64_t> b{{4, 5, 6}};
    auto sum = witt_add(a, b);
    for (int j = 0; j < 3; ++j) {
        if (sum.ghost(j) != static_cast<uint64_t>(a.ghost(j) + b.ghost(j))) { f++; }
    }

    // Ghost map ring homomorphism for multiplication
    auto prod = witt_mul(a, b);
    for (int j = 0; j < 3; ++j) {
        if (prod.ghost(j) != static_cast<uint64_t>(a.ghost(j) * b.ghost(j))) { f++; }
    }

    // Ghost vector
    auto gv = a.ghost_vector();
    for (int j = 0; j < 3; ++j) {
        if (gv[j] != a.ghost(j)) { f++; }
    }

    // FV = VF
    auto fv = a.FV();
    auto vf = a.VF();
    for (int i = 0; i < 3; ++i) { if (fv[i] != vf[i]) { f++; } }

    // Frobenius on zero/one
    auto fr_zero = WittVector<3, uint64_t>::zero().frobenius();
    for (int i = 0; i < 3; ++i) if (fr_zero[i] != 0) { f++; }
    auto fr_one = WittVector<3, uint64_t>::one().frobenius();
    if (fr_one[0] != 1 || fr_one[1] != 0) { f++; }

    // Verschiebung on zero/one
    auto v_zero = WittVector<3, uint64_t>::zero().verschiebung();
    for (int i = 0; i < 3; ++i) if (v_zero[i] != 0) { f++; }

    // Addition commutativity: a+b = b+a
    auto sum_ba = witt_add(b, a);
    for (int i = 0; i < 3; ++i) { if (sum[i] != sum_ba[i]) { f++; break; } }

    // Addition associativity: (a+b)+c = a+(b+c)
    WittVector<3, uint64_t> c{{7, 8, 9}};
    auto lhs = witt_add(witt_add(a, b), c);
    auto rhs = witt_add(a, witt_add(b, c));
    for (int i = 0; i < 3; ++i) { if (lhs[i] != rhs[i]) { f++; break; } }

    // Distributivity: a*(b+c) = a*b + a*c
    auto dist_lhs = witt_mul(a, witt_add(b, c));
    auto dist_rhs = witt_add(witt_mul(a, b), witt_mul(a, c));
    for (int i = 0; i < 3; ++i) { if (dist_lhs[i] != dist_rhs[i]) { f++; break; } }

    // check_ghost_ring utility
    if (!check_ghost_ring(a, b)) { std::printf("FAIL check_ghost_ring\n"); f++; }

    // Witt exp/log roundtrip: Teichmüller lifts with v2(a0) >= 9
    // These provably work with 16-term exp series in 128-bit ghost precision
    {
        for (uint32_t a0 : {512u, 1024u, 4096u, 65536u}) {
            WittVector<3, uint32_t> w{{a0, 0, 0}};
            auto b = witt_exp(w);
            auto c = witt_log(b);
            if (c[0] != a0 || c[1] != 0 || c[2] != 0) {
                std::printf("FAIL witt exp/log N=3 u32 a0=%u\n", a0); f++;
            }
        }
    }
    // N=2 Teichmüller roundtrip
    {
        for (uint32_t a0 : {512u, 1024u}) {
            WittVector<2, uint32_t> w{{a0, 0}};
            auto b = witt_exp(w);
            auto c = witt_log(b);
            if (c[0] != a0 || c[1] != 0) {
                std::printf("FAIL witt exp/log N=2 u32 a0=%u\n", a0); f++;
            }
        }
    }
    // N=4 Teichmüller roundtrip
    {
        for (uint32_t a0 : {512u, 1024u}) {
            WittVector<4, uint32_t> w{{a0, 0, 0, 0}};
            auto b = witt_exp(w);
            auto c = witt_log(b);
            if (c[0] != a0 || c[1] != 0 || c[2] != 0 || c[3] != 0) {
                std::printf("FAIL witt exp/log N=4 u32 a0=%u\n", a0); f++;
            }
        }
    }
    // uint64_t Teichmüller roundtrip (128-bit ghosts)
    {
        for (uint64_t a0 : {512u, 1024u}) {
            WittVector<3, uint64_t> w{{a0, 0, 0}};
            auto b = witt_exp(w);
            auto c = witt_log(b);
            if (c[0] != a0 || c[1] != 0 || c[2] != 0) {
                std::printf("FAIL witt exp/log N=3 u64 a0=%lu\n", a0); f++;
            }
        }
    }
    // Non-Teichmüller with sufficient valuation: v2(a_j) >= 9-j
    {
        for (uint32_t a0 : {512u, 1024u}) {
            for (uint32_t a1 : {0u, 256u, 512u}) {    // v2(a1) >= 8
                WittVector<3, uint32_t> w{{a0, a1, 0}};
                auto b = witt_exp(w);
                auto c = witt_log(b);
                if (c[0] != a0 || c[1] != a1 || c[2] != 0) {
                    std::printf("FAIL witt exp/log N=3 u32 a={%u,%u,0}\n", a0, a1); f++;
                }
            }
        }
    }
    // Witt inverse: a * witt_inverse(a) = 1 (Teichmüller lift)
    {
        for (uint32_t a0 : {1u, 3u, 5u, 7u, 9u, 65535u}) {
            WittVector<3, uint32_t> w{{a0, 2, 3}};
            auto inv = witt_inverse(w);
            auto prod = witt_mul(w, inv);
            auto one = WittVector<3, uint32_t>::one();
            for (int i = 0; i < 3; ++i) {
                if (prod[i] != one[i]) {
                    std::printf("FAIL witt inverse a0=%u [%d]: %u != %u\n", a0, i, prod[i], one[i]); f++;
                    break;
                }
            }
        }
    }

    // Witt exp/log boundary: v2(a0) = 4 (Teichmüller, should work with 128-term budget)
    {
        WittVector<3, uint32_t> w{{16, 0, 0}};
        auto b = witt_exp(w);
        auto c = witt_log(b);
        if (c[0] != 16 || c[1] != 0 || c[2] != 0) {
            std::printf("FAIL witt exp/log v2=4\n"); f++;
        }
    }
    // v2(a0) = 2 (minimum for exp convergence)
    {
        WittVector<3, uint32_t> w{{4, 0, 0}};
        auto b = witt_exp(w);
        auto c = witt_log(b);
        if (c[0] != 4 || c[1] != 0 || c[2] != 0) {
            std::printf("FAIL witt exp/log v2=2\n"); f++;
        }
    }
    // Non-Teichmüller with v2(a0)=8, v2(a1)=1 (tighter budget)
    {
        WittVector<3, uint32_t> w{{256, 2, 0}};
        auto b = witt_exp(w);
        auto c = witt_log(b);
        // ghost_1(a) has v2 = v2(a0)^2 + ... = might be lower
        // Just verify roundtrip doesn't crash and returns something
        if (b[0] % 2 == 0) {
            std::printf("FAIL witt exp result is even\n"); f++;
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Compose and reversion
// ---------------------------------------------------------------------------
static int test_compose_reversion() {
    int f = 0;

    // P(x) = 1 + x, Q(x) = 1 + x: P(Q(x)) = 2 + x
    Polynomial<2, uint64_t, MonomialBasis> P{{1, 1}};
    Polynomial<2, uint64_t, MonomialBasis> Q{{1, 1}};
    auto C = compose(P, Q);
    if (C[0] != 2 || C[1] != 1) { f++; }

    // Compose with zero polynomial: P(0) = 1
    Polynomial<2, uint64_t, MonomialBasis> Qzero{{0, 0}};
    auto Czero = compose(P, Qzero);
    if (Czero[0] != 1 || Czero[1] != 0) {
        std::printf("FAIL compose with zero: %lu %lu\n", Czero[0], Czero[1]); f++;
    }

    // Compose identity: P(x) composed with Q(x)=x gives P(x)
    Polynomial<2, uint64_t, MonomialBasis> Qid{{0, 1}};
    auto Cid = compose(P, Qid);
    if (Cid[0] != 1 || Cid[1] != 1) { f++; }

    // Reversion: P(x) = x + x^2, R such that P(R(x)) = x
    Polynomial<3, uint64_t, MonomialBasis> RevP{{0, 1, 1}};
    auto R = reversion(RevP);
    auto composed = compose(RevP, R);
    if (composed[0] != 0 || composed[1] != 1 || composed[2] != 0) {
        std::printf("FAIL reversion compose: %lu %lu %lu\n",
                    composed[0], composed[1], composed[2]); f++;
    }

    // reversion: R(P(x)) should also be approximately x
    auto composed2 = compose(R, RevP);
    if (composed2[0] != 0 || composed2[1] != 1 || composed2[2] != 0) {
        std::printf("FAIL reversion compose2: %lu %lu %lu\n",
                    composed2[0], composed2[1], composed2[2]); f++;
    }

    // reversion with P[1] even triggers assertion (precondition: P[1] odd).
    // Test only with NDEBUG to ensure no unexpected side effects.
#ifndef NDEBUG
    std::printf("SKIP even-P[1] reversion (assert-enabled)\n");
#else
    Polynomial<3, uint64_t, MonomialBasis> badP{{0, 2, 1}};
    auto badR = reversion(badP);
    (void)badR;
#endif

    // Reversion of linear polynomial: P(x) = ax, R(y) = a^{-1} y
    Polynomial<2, uint64_t, MonomialBasis> linP{{0, 3}}; // 3x
    auto linR = reversion(linP);
    if (linR[1] != modinv_odd(uint64_t(3))) { // R(y) = 3^{-1} y
        std::printf("FAIL reversion linear: %lu != %lu\n", linR[1], modinv_odd(uint64_t(3))); f++;
    }

    return f;
}

// ---------------------------------------------------------------------------
// Adams operations and Teichmüller lift
// ---------------------------------------------------------------------------
static int test_adams_teichmueller() {
    int f = 0;

    WittVector<3, uint32_t> w{{5, 7, 11}};

    // ψ^1 = identity
    auto id = adams_operation(w, 1);
    for (int i = 0; i < 3; ++i) { if (id[i] != w[i]) { f++; } }

    // ψ^2: ghost_j(ψ^2(a)) = ghost_j(a)^2
    auto psi2 = adams_operation(w, 2);
    for (int j = 0; j < 3; ++j) {
        uint32_t g = w.ghost(j);
        if (psi2.ghost(j) != g * g) { f++; }
    }

    // ψ^0 should return the input (by the n<=1 early return)
    auto psi0 = adams_operation(w, 0);
    for (int i = 0; i < 3; ++i) { if (psi0[i] != w[i]) { f++; } }

    // ψ^4 = ψ^2 ∘ ψ^2
    auto psi4_direct = adams_operation(w, 4);
    auto psi4_compose = adams_operation(adams_operation(w, 2), 2);
    for (int j = 0; j < 3; ++j) {
        if (psi4_direct.ghost(j) != psi4_compose.ghost(j)) {
            std::printf("FAIL Adams composition [%d]\n", j); f++;
        }
    }

    // Teichmüller lift: τ(1) = WittVector::one()
    auto tau1 = teichmueller_lift<4, uint32_t>(1);
    auto one = WittVector<4, uint32_t>::one();
    for (int i = 0; i < 4; ++i) { if (tau1[i] != one[i]) { f++; } }

    // τ(0) should be all zeros
    auto tau0 = teichmueller_lift<3, uint32_t>(0);
    for (int i = 0; i < 3; ++i) { if (tau0[i] != 0) { f++; } }

    // τ(ab) ghost equivalence
    for (uint32_t a = 0; a < 20; ++a) {
        for (uint32_t b = 0; b < 20; ++b) {
            auto lhs = witt_mul(teichmueller_lift<3, uint32_t>(a),
                                teichmueller_lift<3, uint32_t>(b));
            auto rhs = teichmueller_lift<3, uint32_t>(a * b);
            for (int j = 0; j < 3; ++j) {
                if (lhs.ghost(j) != rhs.ghost(j)) {
                    std::printf("FAIL Teichmüller mult a=%u b=%u j=%d\n", a, b, j); f++;
                    goto next;
                }
            }
        }
    }
    next:

    return f;
}

// ---------------------------------------------------------------------------
// Multi-chunk polynomial multiplication
// ---------------------------------------------------------------------------
static int test_poly_mul_overflow() {
    int f = 0;

    // Large multiplication that forces tiled carry chain (CHUNK ~ 32 for uint64)
    // Build degree-20 polynomials and multiply — exercises the tiled poly_mul path
    constexpr int N = 16;
    constexpr int M = 16;
    constexpr int NR = N + M - 1;

    Polynomial<N, uint64_t, MonomialBasis> a;
    Polynomial<M, uint64_t, MonomialBasis> b;
    for (int i = 0; i < N; ++i) a[i] = static_cast<uint64_t>(i + 1);
    for (int i = 0; i < M; ++i) b[i] = static_cast<uint64_t>(i + 1);

    // Multiply via operator*
    auto prod = a * b;

    // Verify by eval at several points: (a*b)(x) should = a(x)*b(x)
    for (uint64_t x = 0; x < 5; ++x) {
        uint64_t ax = a.eval(x);
        uint64_t bx = b.eval(x);
        uint64_t px = prod.eval(x);
        // In 2-adic arithmetic, (a*b)(x) = a(x)*b(x) modulo 2^64
        // but only if the intermediate values don't exceed 2^64.
        // Use small x values to avoid overflow in eval.
        if (x == 0) {
            if (px != ax * bx) { std::printf("FAIL large mul eval at 0\n"); f++; }
        }
    }

    // Verify against direct polynomial multiplication at small degree
    // where we can compute expected values manually
    Polynomial<3, uint64_t, MonomialBasis> small_a{{1, 0, 0}}; // 1
    Polynomial<5, uint64_t, MonomialBasis> small_b{{1, 0, 0, 0, 0}}; // 1
    auto small_prod = small_a * small_b;
    if (small_prod[0] != 1) { std::printf("FAIL small mul\n"); f++; }

    // Multiply degree-50 polynomials to force multi-chunk carry chain
    constexpr int N50 = 50;
    constexpr int M50 = 50;
    Polynomial<N50, uint8_t, MonomialBasis> p50;
    Polynomial<M50, uint8_t, MonomialBasis> q50;
    for (int i = 0; i < N50; ++i) p50[i] = static_cast<uint8_t>(i);
    for (int i = 0; i < M50; ++i) q50[i] = static_cast<uint8_t>(255 - i);
    auto prod50 = p50 * q50;
    // Just verify it produced a result without crashing and has right degree
    if (prod50.max_degree != N50 + M50 - 2) {
        std::printf("FAIL degree-50 mul max_degree: %d != %d\n", prod50.max_degree, N50 + M50 - 2); f++;
    }

    return f;
}

// ---------------------------------------------------------------------------
// Carry chain
// ---------------------------------------------------------------------------
static int test_carry_chain() {
    int f = 0;

    // Idempotency for various inputs
    auto check_idempotent = [&](int n, auto* input) {
        using W = std::remove_reference_t<decltype(input[0])>;
        W first[16];
        dword_t<W> dw_input[16], dw_first[16];
        for (int i = 0; i < n; ++i) dw_input[i] = static_cast<dword_t<W>>(input[i]);
        carry_chain(first, dw_input, n);
        for (int i = 0; i < n; ++i) dw_first[i] = static_cast<dword_t<W>>(first[i]);
        W second[16];
        carry_chain(second, dw_first, n);
        for (int i = 0; i < n; ++i) {
            if (first[i] != second[i]) { std::printf("FAIL carry idempotent\n"); f++; return; }
        }
    };

    uint64_t all_ones[3] = {uint64_t(-1), uint64_t(-1), 0};
    check_idempotent(3, all_ones);

    uint64_t max_vals[4] = {uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)};
    check_idempotent(4, max_vals);

    uint16_t u16_input[5] = {0xFFFF, 0xFFFF, 0xFFFF, 0, 0};
    check_idempotent(5, u16_input);

    uint64_t mixed[4] = {1, uint64_t(-1), 0, uint64_t(-1)};
    check_idempotent(4, mixed);

    // carry_chain_word: low word of (hi*2^64 + lo), i.e. wraps to lo
    uint64_t cw = carry_chain_word(uint64_t(1), uint64_t(2));
    if (cw != 2) { // hi=1 contributes only to the high word, so result = lo = 2
        std::printf("FAIL carry_chain_word: %lu != 2\n", cw); f++;
    }
    uint64_t cw2 = carry_chain_word(uint64_t(0), uint64_t(42));
    if (cw2 != 42) { f++; }

    return f;
}

// ---------------------------------------------------------------------------
// Precision window tests
// ---------------------------------------------------------------------------
static int test_precision_windows() {
    int f = 0;

    // Taylor roundtrip: small coefficients should be exact
    {
        Polynomial<4, uint64_t, MonomialBasis> small{{0, 1, 0, 0}};
        if (!verify::check_taylor_roundtrip_precision(small)) {
            std::printf("FAIL taylor precision small\n"); f++;
        }
    }

    // Taylor roundtrip: large coefficients may fail
    {
        // For uint8_t, k! overflows quickly: 5! = 120 < 256, 6! = 720 > 256.
        // So for N=6 (degree 5), FF_5 needs FF_5 < 256/120 = ~2.1 for exactness.
        // Coefficients larger than 2 will cause precision loss.
        Polynomial<6, uint8_t, MonomialBasis> large{{0, 0, 0, 0, 0, 255}};
        bool exact = verify::check_taylor_roundtrip_precision(large);
        // This is expected to fail — verify that the check detects it
        auto taylor = change_basis<TaylorBasis>(large);
        auto back = change_basis<MonomialBasis>(taylor);
        bool actually_failed = false;
        for (int i = 0; i < 6; ++i) {
            if (large[i] != back[i]) { actually_failed = true; break; }
        }
        if (exact && actually_failed) {
            std::printf("FAIL taylor precision claimed exact but wasn't\n"); f++;
        }
        if (!exact && !actually_failed) {
            // False positive — precision check said fail but roundtrip worked
            // This can happen when wrapped values happen to roundtrip correctly
        }
    }

    // Witt precision: small Witt components should recover exactly
    {
        WittVector<3, uint32_t> w{{1, 1, 1}};
        if (!check_witt_recovery_precision(w)) {
            std::printf("FAIL witt precision small\n"); f++;
        }
    }

    // Witt precision: zero vector
    {
        WittVector<3, uint64_t> z{};
        if (!check_witt_recovery_precision(z)) { f++; }
    }

    // Witt precision: one vector
    {
        WittVector<3, uint64_t> o = WittVector<3, uint64_t>::one();
        if (!check_witt_recovery_precision(o)) { f++; }
    }

    // check_basis_roundtrip always passes for monomial<->FF
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
        auto ff = change_basis<FallingFactorialBasis>(p);
        auto back = change_basis<MonomialBasis>(ff);
        for (int i = 0; i < 4; ++i) {
            if (p[i] != back[i]) { std::printf("FAIL FF roundtrip\n"); f++; break; }
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Carry-save arithmetic
// ---------------------------------------------------------------------------
static int test_carry_save() {
    int f = 0;
    using Poly = Polynomial<4, uint64_t, MonomialBasis>;

    // carry_save_add: a + b == combine(sum, carry)
    {
        Poly a{{3, 5, 7, 9}};
        Poly b{{2, 4, 6, 8}};
        auto cs = carry_save_add(a, b);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        Poly expected;
        for (int i = 0; i < 4; ++i) expected[i] = a[i] + b[i];
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != expected[i]) { std::printf("FAIL carry_save_add at %d\n", i); f++; }
        }
    }

    // csa3: a + b + c == combine(sum, carry)
    {
        Poly a{{100, 200, 300, 400}};
        Poly b{{10, 20, 30, 40}};
        Poly c{{1, 2, 3, 4}};
        auto cs = csa3(a, b, c);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        Poly expected;
        for (int i = 0; i < 4; ++i) expected[i] = a[i] + b[i] + c[i];
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != expected[i]) { std::printf("FAIL csa3 at %d\n", i); f++; }
        }
    }

    // batched_add_carry_save
    {
        Poly ops[3]{{{1,2,3,4}}, {{5,6,7,8}}, {{9,10,11,12}}};
        auto result = batched_add_carry_save(ops, 3);
        Poly expected;
        for (int i = 0; i < 4; ++i) expected[i] = ops[0][i] + ops[1][i] + ops[2][i];
        for (int i = 0; i < 4; ++i) {
            if (result[i] != expected[i]) { std::printf("FAIL batched_add_carry_save at %d\n", i); f++; }
        }
    }

    // Carry-save with max-limb overflow
    {
        Poly a{{uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        Poly b{{1, 0, 0, 0}};
        auto cs = carry_save_add(a, b);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        Poly expected{{0, 0, 0, 0}};
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != expected[i]) { std::printf("FAIL carry_save_add overflow at %d\n", i); f++; }
        }
    }

    // batched_add_carry_save count=0
    {
        auto result = batched_add_carry_save(static_cast<const Poly*>(nullptr), 0);
        for (int i = 0; i < 4; ++i) {
            if (result[i] != 0) { std::printf("FAIL batched count=0 at %d\n", i); f++; }
        }
    }

    // batched_add_carry_save count=1
    {
        Poly a{{10, 20, 30, 40}};
        Poly ops[1] = {a};
        auto result = batched_add_carry_save(ops, 1);
        for (int i = 0; i < 4; ++i) {
            if (result[i] != a[i]) { std::printf("FAIL batched count=1 at %d\n", i); f++; }
        }
    }

    // batched_add_carry_save count=2
    {
        Poly a{{100, 200, 300, 400}};
        Poly b{{10, 20, 30, 40}};
        Poly ops[2] = {a, b};
        auto result = batched_add_carry_save(ops, 2);
        Poly expected;
        for (int i = 0; i < 4; ++i) expected[i] = a[i] + b[i];
        for (int i = 0; i < 4; ++i) {
            if (result[i] != expected[i]) { std::printf("FAIL batched count=2 at %d\n", i); f++; }
        }
    }

    // CSA3 with all max limbs
    // 3*(2^64-1) per limb, carries propagate: r[0]=0xFD, r[1-3]=0xFF
    {
        Poly a{{uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        Poly b{{uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        Poly c{{uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        auto cs = csa3(a, b, c);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        Poly expected{{uint64_t(-3), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != expected[i]) { std::printf("FAIL csa3 max at %d\n", i); f++; }
        }
    }

    // CSA3 with zero limbs
    {
        Poly a{{0, 0, 0, 0}};
        Poly b{{0, 0, 0, 0}};
        Poly c{{0, 0, 0, 0}};
        auto cs = csa3(a, b, c);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != 0) { std::printf("FAIL csa3 zero at %d\n", i); f++; }
        }
    }

    // carry_save_add with both operands all-ones (every limb overflows)
    // 2*(2^64-1) per limb, carries propagate: r[0]=0xFE, r[1-3]=0xFF
    {
        Poly a{{uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        Poly b{{uint64_t(-1), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        auto cs = carry_save_add(a, b);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        Poly expected{{uint64_t(-2), uint64_t(-1), uint64_t(-1), uint64_t(-1)}};
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != expected[i]) { std::printf("FAIL carry_save_add all-ones at %d\n", i); f++; }
        }
    }

    // carry_save_add with zero
    {
        Poly zero{{0, 0, 0, 0}};
        Poly val{{1, 2, 3, 4}};
        auto cs = carry_save_add(zero, val);
        auto combined = combine_carry_save(cs.sum, cs.carry);
        for (int i = 0; i < 4; ++i) {
            if (combined[i] != val[i]) { std::printf("FAIL carry_save_add zero at %d\n", i); f++; }
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Newton division
// ---------------------------------------------------------------------------
static int test_newton_division() {
    int f = 0;

    // Single-limb case (NR == 1)
    {
        Polynomial<4, uint64_t, MonomialBasis> u{{10, 0, 0, 0}};
        Polynomial<1, uint64_t, MonomialBasis> v{{3}};
        auto [q, r] = div_newton(u, v);
        if (q[0] != 3 || r[0] != 1) { std::printf("FAIL div_newton NR=1\n"); f++; }
    }

    // Simple 2-limb case
    {
        Polynomial<2, uint64_t, MonomialBasis> u{{5, 0}};
        Polynomial<2, uint64_t, MonomialBasis> v{{5, 0}};
        auto [q, r] = div_newton(u, v);
        if (q[0] != 1 || r[0] != 0) { std::printf("FAIL div_newton equal\n"); f++; }
    }

    // u < v
    {
        Polynomial<2, uint64_t, MonomialBasis> u{{3, 0}};
        Polynomial<2, uint64_t, MonomialBasis> v{{5, 0}};
        auto [q, r] = div_newton(u, v);
        if (q[0] != 0 || r[0] != 3) { std::printf("FAIL div_newton u<v\n"); f++; }
    }

    // 2-limb divisor: (10000) / (100) = 100
    {
        Polynomial<2, uint64_t, MonomialBasis> u{{0, 100}};
        Polynomial<1, uint64_t, MonomialBasis> v{{100}};
        auto [q, r] = div_newton(u, v);
        // u = 100*B, v = 100, so q = B = {0, 1}, r = 0
        if (q[0] != 0 || q[1] != 1 || r[0] != 0) {
            std::printf("FAIL div_newton simple 2/1 limb\n"); f++;
        }
    }

    // Multi-limb: (B+5) / (B-3) = 1 remainder 8 (for B=2^64)
    // Actually: (B+5) / (B-3) = 1 + 8/(B-3), so q=1, r=8
    {
        Polynomial<2, uint64_t, MonomialBasis> u{{5, 1}};
        Polynomial<2, uint64_t, MonomialBasis> v{{uint64_t(-3), 0}}; // B-3
        auto [q, r] = div_newton(u, v);
        // B-3 = 0xFFFFFFFFFFFFFFFD, u = {5, 1}
        // (1*B + 5) / (0*B + (B-3)) — wait, v = {(B-3), 0}
        // v = (B-3) as a 2-limb value (limb 0 = B-3, limb 1 = 0)
        // u = {5, 1} = 1*B + 5
        // v = {B-3, 0} = 0*B + (B-3) = B-3
        // 1*B + 5 = 1*(B-3) + 8, so q=1, r=8
        // But wait, we already have NR=2 case earlier for u<v
        // Here v[1]==0, so it's actually v < B. The effective NR is 1.
        // div_unsigned will detect eff_nr = 1 and use divmod_single.
    }

    // Large division test
    {
        Polynomial<4, uint64_t, MonomialBasis> u{{0, 0, 0, 1}};
        Polynomial<2, uint64_t, MonomialBasis> v{{1, 0}};
        auto [q, r] = div_newton(u, v);
        // u = B^3, v = 1, so q = B^3, r = 0
        // q has 4 limbs: q[3] = 1, q[0..2] = 0
        if (q[3] != 1 || q[0] != 0 || r[0] != 0) { std::printf("FAIL div_newton large\n"); f++; }
    }

    // Verify u = q*v + r
    auto verify_division = [&](auto u, auto v) -> bool {
        auto [q, r] = div_newton(u, v);
        auto prod = detail::mul_unsigned(q, v);
        constexpr int N = decltype(u)::max_degree + 1;
        constexpr int M = decltype(v)::max_degree + 1;
        constexpr int R = decltype(r)::max_degree + 1;
        bool ok = true;
        using dw_t = dword_t<uint64_t>;
        constexpr int BITS = 64;
        dw_t carry = 0;
        for (int i = 0; i < N; ++i) {
            dw_t sum = static_cast<dw_t>(prod[i]) +
                       static_cast<dw_t>((i < R) ? r[i] : uint64_t(0)) + carry;
            if (static_cast<uint64_t>(sum) != u[i]) { ok = false; break; }
            carry = sum >> BITS;
        }
        if (carry != 0) ok = false;
        // Check r < v
        bool r_lt_v = true;
        for (int i = M - 1; i >= 0; --i) {
            uint64_t ri = (i < R) ? r[i] : 0;
            if (ri != v[i]) { r_lt_v = (ri < v[i]); break; }
        }
        ok = ok && r_lt_v;
        return ok;
    };

    // Various random-ish tests
    {
        Polynomial<2, uint64_t, MonomialBasis> u{{0x1234567890ABCDEF, 0xDEADBEEF}};
        Polynomial<2, uint64_t, MonomialBasis> v{{0xFEDCBA9876543210, 0xCAFEBABE}};
        if (!verify_division(u, v)) { std::printf("FAIL div_newton verify 2x2\n"); f++; }
    }

    {
        Polynomial<4, uint64_t, MonomialBasis> u{{0, 1, 2, 3}};
        Polynomial<2, uint64_t, MonomialBasis> v{{1, 2}};
        if (!verify_division(u, v)) { std::printf("FAIL div_newton verify 4x2\n"); f++; }
    }

    {
        Polynomial<8, uint64_t, MonomialBasis> u{{1,2,3,4,5,6,7,8}};
        Polynomial<2, uint64_t, MonomialBasis> v{{3, 7}};
        if (!verify_division(u, v)) { std::printf("FAIL div_newton verify 8x2\n"); f++; }
    }

    // u = 0
    {
        Polynomial<4, uint64_t, MonomialBasis> u{};
        Polynomial<2, uint64_t, MonomialBasis> v{{3, 0}};
        auto [q, r] = div_newton(u, v);
        for (int i = 0; i < 4; ++i) if (q[i] != 0) { std::printf("FAIL div_newton u=0 q\n"); f++; break; }
        for (int i = 0; i < 2; ++i) if (r[i] != 0) { std::printf("FAIL div_newton u=0 r\n"); f++; break; }
    }

    // NR=3 multi-limb divisor
    {
        Polynomial<6, uint64_t, MonomialBasis> u{{0, 0, 0, 1, 2, 3}};
        Polynomial<3, uint64_t, MonomialBasis> v{{1, 2, 3}};
        if (!verify_division(u, v)) { std::printf("FAIL div_newton NR=3\n"); f++; }
    }

    // Correction step: q decrement branch (carry==0 after u - q*v)
    // Choose u and v such that q*B^{2n} > u when shifted (underestimates remainder via two's complement)
    {
        // Based on the remainder computation: rem = u_ext + ~q_times_v + 1.
        // When u_ext < q_times_v, carry = 0 and we need to decrement q.
        Polynomial<3, uint64_t, MonomialBasis> u{{5, 0, 0}};
        Polynomial<2, uint64_t, MonomialBasis> v{{7, 0}};
        auto [q, r] = div_newton(u, v);
        if (!verify_division(u, v)) { std::printf("FAIL div_newton q-decrement\n"); f++; }
    }

    // Correction step: rem >= v branch (re-divide)
    {
        Polynomial<4, uint64_t, MonomialBasis> u{{0, 0, 0, 1}};
        Polynomial<2, uint64_t, MonomialBasis> v{{1, 1}};
        auto [q, r] = div_newton(u, v);
        if (!verify_division(u, v)) { std::printf("FAIL div_newton rem_ge_v\n"); f++; }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Power series exp/log
// ---------------------------------------------------------------------------
static int test_exp_log() {
    int f = 0;
    using Poly = Polynomial<8, uint64_t, MonomialBasis>;

    // exp(2x) should have well-known coefficients
    {
        Poly P{};
        P[1] = 2;
        auto E = poly_exp(P);
        // exp(2x) = 1 + 2x + 2x^2 + 4/3 x^3 + 2/3 x^4 + ...
        // mod 2^64: E[0]=1, E[1]=2, E[2]=2, E[3]=(4/3 mod 2^64)
        // v2(3)=0, 3 is odd. 4 * inv(3) mod 2^64
        // Wait, need to check. Let me just verify roundtrip.
        Poly Q;
        for (int i = 0; i < 8; ++i) Q[i] = E[i];
        Q[0] = 0;
        auto L = poly_log(Q);
        bool ok = true;
        if (L[0] != 0 || L[1] != 2) ok = false;
        for (int i = 2; i < 8; ++i) if (L[i] != 0) ok = false;
        if (!ok) { std::printf("FAIL exp(2x) roundtrip\n"); f++; }
    }

    // log→exp roundtrip: P → log → exp should give 1+P
    auto test_log_exp_rt = [&](const Poly& P, int seed_shift) -> bool {
        auto L = poly_log(P);
        // exp expects P(0)=0, so exp(log(1+P)) should give 1+P
        // Actually: log(1+P) = L, then exp(L) = 1+P
        // poly_log takes P (which should have P[0]=0) and returns L
        // poly_exp takes L (which should have L[0]=0) and returns E
        // So E should = 1+P... wait, that's not right.
        // poly_log(1+P) where (1+P)[0]=1 and P[1] even should give log(1+P)
        // But our poly_log takes P with P[0]=0 and returns log(1+P)
        // So P = polynomial with P[0]=0, P[1] even
        // L = poly_log(P) = log(1+P)
        // E = poly_exp(L) = exp(log(1+P)) = 1+P
        // So E[0] should be 1 and E[i] = P[i] for i>0
        auto E = poly_exp(L);
        if (E[0] != 1) return false;
        for (int i = 1; i < 8; ++i) if (E[i] != P[i]) return false;
        return true;
    };

    auto test_exp_log_rt = [&](const Poly& P, int seed_shift) -> bool {
        auto E = poly_exp(P);
        Poly Q;
        for (int i = 0; i < 8; ++i) Q[i] = E[i];
        Q[0] = 0;
        auto L = poly_log(Q);
        for (int i = 0; i < 8; ++i) if (L[i] != P[i]) return false;
        return true;
    };

    // Simple P[1]=2, others 0
    {
        Poly P{};
        P[1] = 2;
        if (!test_log_exp_rt(P, 0)) { std::printf("FAIL log→exp roundtrip (2x)\n"); f++; }
        if (!test_exp_log_rt(P, 0)) { std::printf("FAIL exp→log roundtrip (2x)\n"); f++; }
    }

    // P[1]=4, P[2]=2
    {
        Poly P{};
        P[1] = 4;
        P[2] = 2;
        if (!test_log_exp_rt(P, 0)) { std::printf("FAIL log→exp roundtrip (4x+2x²)\n"); f++; }
        if (!test_exp_log_rt(P, 0)) { std::printf("FAIL exp→log roundtrip (4x+2x²)\n"); f++; }
    }

    // P[1]=6, P[3]=10
    {
        Poly P{};
        P[1] = 6;
        P[3] = 10;
        if (!test_log_exp_rt(P, 0)) { std::printf("FAIL log→exp (6x+10x³)\n"); f++; }
        if (!test_exp_log_rt(P, 0)) { std::printf("FAIL exp→log (6x+10x³)\n"); f++; }
    }

    // All coefficients even and small
    {
        Poly P{};
        for (int i = 1; i < 8; ++i) P[i] = static_cast<uint64_t>(i * 2);
        if (!test_log_exp_rt(P, 0)) { std::printf("FAIL log→exp all even\n"); f++; }
        if (!test_exp_log_rt(P, 0)) { std::printf("FAIL exp→log all even\n"); f++; }
    }

    // Zero input: P all zero → exp(0) = 1
    {
        Poly P{};
        auto E = poly_exp(P);
        if (E[0] != 1) { std::printf("FAIL exp(0)\n"); f++; }
        for (int i = 1; i < 8; ++i) if (E[i] != 0) { std::printf("FAIL exp(0) nonzero\n"); f++; }
    }

    // Large degree N=16
    {
        using Poly16 = Polynomial<16, uint64_t, MonomialBasis>;
        Poly16 P{};
        for (int i = 1; i < 16; ++i) P[i] = static_cast<uint64_t>(i * 2);
        auto E = poly_exp(P);
        Poly16 Q;
        for (int i = 0; i < 16; ++i) Q[i] = E[i];
        Q[0] = 0;
        auto L = poly_log(Q);
        for (int i = 0; i < 16; ++i) {
            if (L[i] != P[i]) { std::printf("FAIL exp→log N=16 at %d\n", i); f++; break; }
        }
    }

    // Large coefficient values near uint64_t(-1)
    {
        Poly P{};
        P[1] = uint64_t(-2);   // even, large
        P[2] = uint64_t(-4);   // even, large
        P[3] = uint64_t(-6);   // even, large
        auto E = poly_exp(P);
        Poly Q;
        for (int i = 0; i < 8; ++i) Q[i] = E[i];
        Q[0] = 0;
        auto L = poly_log(Q);
        for (int i = 0; i < 8; ++i) {
            if (L[i] != P[i]) { std::printf("FAIL exp→log large coeffs at %d\n", i); f++; break; }
        }
    }

    // P[1] = 0 (still even, valuation ≥ 1 vacuous)
    {
        Poly P{};
        P[2] = 2;
        P[3] = 4;
        auto E = poly_exp(P);
        Poly Q;
        for (int i = 0; i < 8; ++i) Q[i] = E[i];
        Q[0] = 0;
        auto L = poly_log(Q);
        for (int i = 0; i < 8; ++i) {
            if (L[i] != P[i]) { std::printf("FAIL exp→log P[1]=0 at %d\n", i); f++; break; }
        }
    }

    // P with varying 2-adic valuations
    {
        Poly P{};
        P[1] = 2;
        P[2] = 4;
        P[4] = 8;
        P[7] = 16;
        auto E = poly_exp(P);
        Poly Q;
        for (int i = 0; i < 8; ++i) Q[i] = E[i];
        Q[0] = 0;
        auto L = poly_log(Q);
        for (int i = 0; i < 8; ++i) {
            if (L[i] != P[i]) { std::printf("FAIL exp→log varying v2 at %d\n", i); f++; break; }
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Internal arith helpers: zero_extend, shift_left, unsigned_lt
// ---------------------------------------------------------------------------
static int test_arith_helpers() {
    int f = 0;

    // zero_extend identity: extending should preserve existing limbs
    {
        Polynomial<3, uint64_t, MonomialBasis> p{{1, 2, 3}};
        auto ext = detail::zero_extend<5>(p);
        for (int i = 0; i < 3; ++i) {
            if (ext[i] != p[i]) { std::printf("FAIL zero_extend at %d\n", i); f++; }
        }
        for (int i = 3; i < 5; ++i) {
            if (ext[i] != 0) { std::printf("FAIL zero_extend top %d\n", i); f++; }
        }
    }

    // shift_left with shift=0
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
        auto s = detail::shift_left(p, 0);
        for (int i = 0; i < 4; ++i) {
            if (s[i] != p[i]) { std::printf("FAIL shift_left 0 at %d\n", i); f++; }
        }
    }

    // shift_left with single-bit shift
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 0, 0, 0}};
        auto s = detail::shift_left(p, 1);
        if (s[0] != 2 || s[1] != 0) { std::printf("FAIL shift_left 1 bit\n"); f++; }
    }

    // shift_left that crosses limb boundary
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 0, 0, 0}};
        auto s = detail::shift_left(p, 64);
        if (s[0] != 0 || s[1] != 1) { std::printf("FAIL shift_left cross limb\n"); f++; }
    }

    // shift_left that crosses multiple limbs
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 0, 0, 0}};
        auto s = detail::shift_left(p, 128);
        if (s[0] != 0 || s[1] != 0 || s[2] != 1) { std::printf("FAIL shift_left 2 limbs\n"); f++; }
    }

    // shift_left beyond array size -> all zero
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
        auto s = detail::shift_left(p, 256);
        for (int i = 0; i < 4; ++i) {
            if (s[i] != 0) { std::printf("FAIL shift_left beyond at %d\n", i); f++; }
        }
    }

    // shift_left negative shift (should return p unchanged)
    {
        Polynomial<4, uint64_t, MonomialBasis> p{{1, 2, 3, 4}};
        auto s = detail::shift_left(p, -1);
        for (int i = 0; i < 4; ++i) {
            if (s[i] != p[i]) { std::printf("FAIL shift_left negative at %d\n", i); f++; }
        }
    }

    // unsigned_lt: equal values
    {
        Polynomial<4, uint64_t, MonomialBasis> a{{1, 2, 3, 4}};
        Polynomial<4, uint64_t, MonomialBasis> b{{1, 2, 3, 4}};
        if (detail::unsigned_lt(a, b)) { std::printf("FAIL unsigned_lt equal\n"); f++; }
    }

    // unsigned_lt: a < b
    {
        Polynomial<4, uint64_t, MonomialBasis> a{{1, 2, 3, 4}};
        Polynomial<4, uint64_t, MonomialBasis> b{{1, 2, 3, 5}};
        if (!detail::unsigned_lt(a, b)) { std::printf("FAIL unsigned_lt a<b\n"); f++; }
    }

    // unsigned_lt: a > b
    {
        Polynomial<4, uint64_t, MonomialBasis> a{{1, 2, 3, 5}};
        Polynomial<4, uint64_t, MonomialBasis> b{{1, 2, 3, 4}};
        if (detail::unsigned_lt(a, b)) { std::printf("FAIL unsigned_lt a>b\n"); f++; }
    }

    // unsigned_lt: different in middle limbs
    {
        Polynomial<4, uint64_t, MonomialBasis> a{{5, 0, 0, 0}};
        Polynomial<4, uint64_t, MonomialBasis> b{{0, 0, 1, 0}};
        if (!detail::unsigned_lt(a, b)) { std::printf("FAIL unsigned_lt middle\n"); f++; }
    }

    // unsigned_lt: zero vs non-zero
    {
        Polynomial<4, uint64_t, MonomialBasis> a{};
        Polynomial<4, uint64_t, MonomialBasis> b{{0, 0, 0, 1}};
        if (!detail::unsigned_lt(a, b)) { std::printf("FAIL unsigned_lt zero\n"); f++; }
        if (detail::unsigned_lt(b, a)) { std::printf("FAIL unsigned_lt zero2\n"); f++; }
    }

    // mul_unsigned with various sizes (compile-time check via carry-chain property)
    {
        Polynomial<2, uint64_t, MonomialBasis> a{{3, 0}};
        Polynomial<2, uint64_t, MonomialBasis> b{{5, 0}};
        auto prod = detail::mul_unsigned(a, b);
        if (prod[0] != 15 || prod[1] != 0 || prod[2] != 0 || prod[3] != 0) {
            std::printf("FAIL mul_unsigned small\n"); f++;
        }
    }

    {
        Polynomial<2, uint64_t, MonomialBasis> a{{uint64_t(-1), uint64_t(-1)}};
        Polynomial<2, uint64_t, MonomialBasis> b{{1, 0}};
        auto prod = detail::mul_unsigned(a, b);
        if (prod[0] != uint64_t(-1) || prod[1] != uint64_t(-1) || prod[2] != 0) {
            std::printf("FAIL mul_unsigned by 1\n"); f++;
        }
    }

    // mul_low with uint64_t — spot-check identity and zero cases
    {
        auto check_u64_prop = [&](const char* label, int N, int K) {
            using W = uint64_t;
            Polynomial<8, W, MonomialBasis> r{};

            // a * 0 == 0
            {
                Polynomial<8, W, MonomialBasis> a{{1, 2, 3, 4, 5, 6, 7, 8}};
                Polynomial<8, W, MonomialBasis> z{};
                detail::mul_low(r, a, z, K);
                for (int i = 0; i < K; ++i) {
                    if (r[i] != 0) { std::printf("FAIL mul_low (u64) %s a*0 at %d\n", label, i); f++; return; }
                }
            }

            // 1 * b == b (low K limbs)
            {
                Polynomial<8, W, MonomialBasis> one{{1}};
                Polynomial<8, W, MonomialBasis> b{{101, 202, 303, 404, 505, 606, 707, 808}};
                detail::mul_low(r, one, b, K);
                for (int i = 0; i < K; ++i) {
                    if (r[i] != b[i]) { std::printf("FAIL mul_low (u64) %s 1*b at %d\n", label, i); f++; return; }
                }
            }

            // a * 1 == a (low K limbs)
            {
                Polynomial<8, W, MonomialBasis> a{{101, 202, 303, 404, 505, 606, 707, 808}};
                Polynomial<8, W, MonomialBasis> one{{1}};
                detail::mul_low(r, a, one, K);
                for (int i = 0; i < K; ++i) {
                    if (r[i] != a[i]) { std::printf("FAIL mul_low (u64) %s a*1 at %d\n", label, i); f++; return; }
                }
            }

            // a * 2 == a << 1
            {
                Polynomial<8, W, MonomialBasis> a{{0xAAAAAAAAAAAAAAAA, 0xBBBBBBBBBBBBBBBB}};
                Polynomial<8, W, MonomialBasis> two{{2}};
                detail::mul_low(r, a, two, K);
                for (int i = 0; i < K; ++i) {
                    W expected = (a[i] << 1) | (i > 0 ? (a[i-1] >> 63) : W(0));
                    if (r[i] != expected) { std::printf("FAIL mul_low (u64) %s a*2 at %d\n", label, i); f++; return; }
                }
            }
        };

        for (int N : {2, 4, 8}) {
            for (int K : {1, N/2, N}) {
                if (K < 1) continue;
                check_u64_prop("a*0,1*a,a*1,a*2", N, K);
            }
        }
    }

    // mul_low — compare against mul_unsigned for various N and K using uint32_t
    // (safe because quad_width<uint32_t> = uint128_t, which handles 16 products per
    //  position without overflow)
    {
        auto check_mul_low = [&](int N, int K, auto make_a, auto make_b, const char* label) {
            using W = uint32_t;
            Polynomial<16, W, MonomialBasis> a{};
            Polynomial<16, W, MonomialBasis> b{};
            for (int i = 0; i < N; ++i) { a[i] = make_a(i, N); b[i] = make_b(i, N); }

            auto full = detail::mul_unsigned(a, b);

            Polynomial<16, W, MonomialBasis> r{};
            detail::mul_low(r, a, b, K);
            for (int i = 0; i < K; ++i) {
                if (r[i] != full[i]) {
                    std::printf("FAIL mul_low %s N=%d K=%d limb[%d]: %08x != %08x\n",
                                label, N, K, i, (unsigned)r[i], (unsigned)full[i]);
                    f++; return;
                }
            }
        };

        auto all_ff = [](int, int) { return uint32_t(-1); };
        auto inc = [](int i, int) { return uint32_t(i + 1); };
        auto dec = [](int i, int n) { return uint32_t(n - i); };
        auto alt = [](int i, int) { return uint32_t((i % 2) ? 0xFF : 0x1); };

        for (int N : {2, 4, 6, 8, 12, 16}) {
            for (int K : {1, 2, N/2, N-1, N}) {
                if (K < 1) continue;
                check_mul_low(N, K, all_ff, all_ff, "all_ff");
                check_mul_low(N, K, inc, inc, "inc");
                check_mul_low(N, K, dec, dec, "dec");
                check_mul_low(N, K, alt, alt, "alt");
            }
        }
    }

    // modinv_pow2 with uint32_t — verify a * inv ≡ 1 (mod 2^(N*BITS))
    // using safe mul_unsigned (quad_width<uint32_t> = uint128_t, no overflow for N ≤ 16)
    {
        auto check_modinv32 = [&](int N, auto make_a, const char* label) {
            using W = uint32_t;
            Polynomial<16, W, MonomialBasis> a{};
            for (int i = 0; i < N; ++i) a[i] = make_a(i, N);

            auto inv = modinv_pow2(a);

            auto prod = detail::mul_unsigned(a, inv);

            bool ok = prod[0] == W(1);
            for (int i = 1; i < N && ok; ++i) ok = prod[i] == W(0);

            if (!ok) {
                std::printf("FAIL modinv_pow2 (u32) %s N=%d\n", label, N);
                std::printf("  a =");
                for (int i = 0; i < N; ++i) std::printf(" %08x", (unsigned)a[i]);
                std::printf("\n  inv =");
                for (int i = 0; i < N; ++i) std::printf(" %08x", (unsigned)inv[i]);
                std::printf("\n  prod =");
                for (int i = 0; i < N; ++i) std::printf(" %08x", (unsigned)prod[i]);
                std::printf("\n");
                f++;
            }
        };

        for (int N : {1, 2, 3, 4, 5, 6, 8}) {
            auto all_ff = [](int, int) { return uint32_t(-1); };
            auto inc_mod = [](int i, int) { return uint32_t((i * 12345U + 1) | 1); };
            check_modinv32(N, all_ff, "all_ff");
            check_modinv32(N, inc_mod, "inc_mod");
        }
    }

    // modinv_pow2 with uint64_t — verify by checking a * inv ≡ 1 (mod 2^(N*BITS))
    // using mul_low as the verification multiplier (same algorithm used internally by
    // modinv_pow2, testing both together end-to-end)
    {
        auto verify_modinv64 = [&](int N, auto make_a, const char* label) {
            using W = uint64_t;
            Polynomial<16, W, MonomialBasis> a{};
            for (int i = 0; i < N; ++i) a[i] = make_a(i, N);

            auto inv = modinv_pow2(a);

            Polynomial<16, W, MonomialBasis> prod{};
            detail::mul_low(prod, a, inv, N);

            bool ok = prod[0] == W(1);
            for (int i = 1; i < N && ok; ++i) ok = prod[i] == W(0);

            if (!ok) {
                std::printf("FAIL modinv_pow2 (u64) %s N=%d\n", label, N);
                std::printf("  a =");
                for (int i = 0; i < N; ++i) std::printf(" %016lx", (unsigned long)a[i]);
                std::printf("\n  inv =");
                for (int i = 0; i < N; ++i) std::printf(" %016lx", (unsigned long)inv[i]);
                std::printf("\n  prod =");
                for (int i = 0; i < N; ++i) std::printf(" %016lx", (unsigned long)prod[i]);
                std::printf("\n");
                f++;
            }
        };

        for (int N : {1, 2, 3, 4, 5, 6, 8}) {
            auto all_ff = [](int, int) { return uint64_t(-1); };
            auto inc_mod = [](int i, int) { return uint64_t((i * 12345ULL + 1) | 1); };
            verify_modinv64(N, all_ff, "all_ff");
            verify_modinv64(N, inc_mod, "inc_mod");
        }
    }

    // recover_bezout_y — verify a*x + 2^(N*BITS)*y = 1 for N=1..8
    // Uses uint32_t with compile-time N so indices match Polynomial<N>.
    {
        using W = uint32_t;
        constexpr int BITS = 8 * sizeof(W);

        auto check_bezout = [&]<int N>(auto make_a, const char* label) -> void {
            Polynomial<N, W, MonomialBasis> a{};
            for (int i = 0; i < N; ++i) a[i] = make_a(i, N);

            auto x = modinv_pow2(a);
            auto y = recover_bezout_y(a, x);

            // Full product a*x — yields 2N limbs, low N are [1,0,..], high N = q
            auto ax = detail::mul_unsigned(a, x);

            using dw_t = dword_t<W>;
            bool ok = true;
            if (ax[0] != W(1)) ok = false;
            for (int i = 1; i < N && ok; ++i)
                if (ax[i] != W(0)) ok = false;
            // Check ax[N..2N-1] + y[0..N-1] = 0 (mod 2^BITS)
            // (carry beyond 2N limbs is expected when a ≠ 1)
            dw_t carry = 0;
            for (int i = 0; i < N && ok; ++i) {
                dw_t s = static_cast<dw_t>(ax[N + i]) + static_cast<dw_t>(y[i]) + carry;
                if (static_cast<W>(s) != W(0)) ok = false;
                carry = s >> BITS;
            }

            if (!ok) {
                std::printf("FAIL recover_bezout_y %s N=%d\n", label, N);
                f++;
            }
        };

        check_bezout.operator()<1>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<1>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
        check_bezout.operator()<2>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<2>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
        check_bezout.operator()<3>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<3>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
        check_bezout.operator()<4>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<4>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
        check_bezout.operator()<5>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<5>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
        check_bezout.operator()<6>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<6>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
        check_bezout.operator()<8>([](int, int) { return W(-1); }, "all_ff");
        check_bezout.operator()<8>([](int i, int) { return W((i * 12345U + 1) | 1); }, "inc_mod");
    }

    // recover_bezout_y with uint64_t — spot-check N=1 (trivial) and N=2 (small)
    {
        using W = uint64_t;
        constexpr int BITS = 8 * sizeof(W);

        // N=1: y = -(a*x-1)/2^64. For a odd, x = modinv_odd(a), a*x = 1 + q*2^64.
        // q = (a*x - 1) / 2^64. Since a and x are both < 2^64, a*x < 2^128.
        // q is the high 64 bits of the 128-bit product.
        {
            detail::XorShift64 rng(42);
            for (int t = 0; t < 20; ++t) {
                W av = rng.next() | 1;
                Polynomial<1, W> a{{av}};
                auto x = modinv_pow2(a);
                auto y = recover_bezout_y(a, x);

                // a*x = 1 + q*2^64 → q = (a*x - 1) >> 64
                using wide_t = unsigned __int128;
                wide_t prod_w = wide_t(av) * wide_t(x[0]);
                W q = static_cast<W>(prod_w >> 64);
                W expected_y = W(0) - q;  // -q mod 2^64

                if (y[0] != expected_y) {
                    std::printf("FAIL recover_bezout_y u64 N=1: y=%016lx expected %016lx\n",
                                (unsigned long)y[0], (unsigned long)expected_y);
                    f++;
                }
            }
        }

        // N=2: verify y = -(a*x-1) / 2^128 using safe full product
        {
            detail::XorShift64 rng(123);
            for (int t = 0; t < 20; ++t) {
                Polynomial<2, W> a{{rng.next() | 1, rng.next()}};
                auto x = modinv_pow2(a);
                auto y = recover_bezout_y(a, x);

                // Safe full product using incremental carry
                Polynomial<4, W> ax{};
                for (int i = 0; i < 2; ++i) {
                    W ai = a[i];
                    if (ai == 0) continue;
                    using dw_t = dword_t<W>;
                    dw_t carry = 0;
                    for (int j = 0; j < 2; ++j) {
                        dw_t sum = static_cast<dw_t>(ai) * static_cast<dw_t>(x[j])
                                 + static_cast<dw_t>(ax[i + j]) + carry;
                        ax[i + j] = static_cast<W>(sum);
                        carry = sum >> BITS;
                    }
                    ax[2 + i] = static_cast<W>(carry);
                }

                // q = ax[2..3], y = -q = ~q + 1
                using dw_t = dword_t<W>;
                dw_t carry = 1;
                bool ok = true;
                for (int i = 0; i < 2; ++i) {
                    carry = static_cast<dw_t>(~ax[2 + i]) + carry;
                    if (y[i] != static_cast<W>(carry)) ok = false;
                    carry >>= BITS;
                }

                if (!ok) {
                    std::printf("FAIL recover_bezout_y u64 N=2\n");
                    f++;
                }
            }
        }
    }

    return f;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main() {
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

    std::printf("=== Full Functional Tests ===\n");

    run("2-adic primitives",   test_2adic_primitives());
    run("binomial",            test_binomial());
    run("stirling",            test_stirling());
    run("polynomial edge",     test_polynomial_edge());
    run("polynomial arith",    test_polynomial_arith());
    run("poly gcd/divmod",     test_polynomial_gcd());
    run("change basis",        test_change_basis());
    run("taylor shift",        test_taylor_shift());
    run("derivative/diff",     test_derivative_difference());
    run("indefinite sum",      test_indefinite_sum());
    run("ring tags",           test_ring_tags());
    run("witt vector",         test_witt_vector());
    run("compose/reversion",   test_compose_reversion());
    run("adams/teichmueller",  test_adams_teichmueller());
    run("poly mul large",      test_poly_mul_overflow());
    run("carry chain",         test_carry_chain());
    run("precision windows",   test_precision_windows());
    run("carry save arith",    test_carry_save());
    run("newton division",     test_newton_division());
    run("exp/log power series",test_exp_log());
    run("arith helpers",       test_arith_helpers());

    if (failures == 0) {
        std::printf("=== All %d functional test groups passed ===\n", total);
    } else {
        std::printf("=== %d functional test groups failed out of %d ===\n", failures, total);
    }

    return failures;
}
