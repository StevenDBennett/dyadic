// test_negatives.cpp — Verify that assert preconditions fire on invalid input.
// Each invalid case is run in a forked child so the assertion SIGABRT
// doesn't kill the test harness.
// POSIX-only (fork/waitpid); skipped on Windows.

#include "dyadic.h"
#include "dyadic_verify.h"
#include <cstdio>

#if defined(_WIN32) || defined(_WIN64)
int main() {
    std::printf("=== Negative Tests (assert preconditions) ===\n");
    std::printf("SKIP (POSIX-only: fork not available on Windows)\n");
    return 0;
}
#else

#include <cstdlib>
#include <sys/wait.h>
#include <unistd.h>

using namespace dyadic;

static int failures = 0;

static void check_assert(const char* name, void (*fn)()) {
    pid_t pid = fork();
    if (pid == 0) {
        fn();
        _exit(0);
    }
    int status;
    waitpid(pid, &status, 0);
    if (WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) {
        std::printf("  PASS  %s (assert fired)\n", name);
    } else if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
        failures++;
        std::printf("  FAIL  %s (no assert)\n", name);
    } else {
        failures++;
        std::printf("  FAIL  %s (unexpected exit)\n", name);
    }
}

// --- 1. poly_divmod_cw: B == 0 ---
static void poly_divmod_cw_Bzero() {
    Polynomial<3, uint64_t, MonomialBasis> A{{1, 2, 3}};
    Polynomial<3, uint64_t, MonomialBasis> B{}; // B == 0
    auto result = poly_divmod_cw(A, B);
    (void)result;
}
// poly_divmod_cw: leading coefficient of B is even ---
static void poly_divmod_cw_evenLC() {
    Polynomial<3, uint64_t, MonomialBasis> A{{1, 2, 3}};
    Polynomial<3, uint64_t, MonomialBasis> B{{0, 2, 4}}; // B[2]=4 (even)
    auto result = poly_divmod_cw(A, B);
    (void)result;
}
// witt_exp: v2(a0) < 2 ---
static void witt_exp_low_valuation() {
    WittVector<3, uint64_t> w;
    w[0] = 2; // v2 = 1, < 2
    witt_exp(w);
}
// witt_log: a0 even ---
static void witt_log_even() {
    WittVector<3, uint64_t> w;
    w[0] = 0; // even
    witt_log(w);
}
// witt_inverse: a0 even ---
static void witt_inverse_even() {
    WittVector<3, uint64_t> w;
    w[0] = 0; // even
    witt_inverse(w);
}
// reversion: P[1] even ---
static void reversion_even_P1() {
    Polynomial<3, uint64_t, MonomialBasis> P{{0, 2, 1}};
    reversion(P);
}
// uint128_t division by zero ---
static void uint128_div_by_zero() {
    detail::uint128_t a(1), b(0);
    auto r = a / b;
    (void)r;
}
// uint128_t modulo by zero ---
static void uint128_mod_by_zero() {
    detail::uint128_t a(1), b(0);
    auto r = a % b;
    (void)r;
}

int main() {
    std::printf("=== Negative Tests (assert preconditions) ===\n");

    check_assert("poly_divmod_cw B==0",         poly_divmod_cw_Bzero);
    check_assert("poly_divmod_cw even LC",      poly_divmod_cw_evenLC);
    check_assert("witt_exp v2(a0)<2",           witt_exp_low_valuation);
    check_assert("witt_log a0 even",            witt_log_even);
    check_assert("witt_inverse a0 even",        witt_inverse_even);
    check_assert("reversion P[1] even",         reversion_even_P1);
    check_assert("uint128_t div by zero",        uint128_div_by_zero);
    check_assert("uint128_t mod by zero",        uint128_mod_by_zero);

    std::printf("=== %s ===\n", failures ? "FAILED" : "All negative tests passed");
    return failures;
}

#endif // _WIN32 || _WIN64
