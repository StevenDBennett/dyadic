// benchmark.cpp — Performance benchmarks for key operations.
// Compile: g++ -std=c++20 -O2 -I.. test/benchmark.cpp -o bench && ./bench

#include "dyadic.h"
#include <cstdio>
#include <chrono>

using namespace dyadic;

// Timer helper
struct Timer {
    using Clock = std::chrono::steady_clock;
    Clock::time_point start;

    Timer() : start(Clock::now()) {}

    double elapsed() const {
        auto end = Clock::now();
        return std::chrono::duration<double>(end - start).count();
    }

    void report(const char* name, int ops) {
        double sec = elapsed();
        if (sec > 0)
            std::printf("  %-30s %8d ops  %8.3f s  %8.0f ops/s\n", name, ops, sec, ops / sec);
        else
            std::printf("  %-30s %8d ops  %8.3f s\n", name, ops, sec);
    }
};

// Accumulate a checksum to prevent dead-code elimination
static volatile uint64_t global_checksum = 0;

template<typename T>
static void sink(const T& val) {
    if constexpr (requires { val[0]; }) {
        global_checksum ^= static_cast<uint64_t>(val[0]);
    } else {
        global_checksum ^= static_cast<uint64_t>(val);
    }
}

// ============================================================================
// Benchmarks
// ============================================================================

static void bench_polynomial_eval(int reps) {
    Timer t;
    for (int r = 0; r < reps; ++r) {
        Polynomial<10, uint64_t, MonomialBasis> p{{1,2,3,4,5,6,7,8,9,10}};
        for (uint64_t x = 0; x < 100; ++x) {
            sink(p.eval(x));
        }
    }
    t.report("poly eval N=10", reps * 100);
}

static void bench_basis_roundtrip(int reps) {
    Timer t;
    Polynomial<10, uint64_t, MonomialBasis> p{{1,2,3,4,5,6,7,8,9,10}};
    for (int r = 0; r < reps; ++r) {
        auto ff = change_basis<FallingFactorialBasis>(p);
        auto back = change_basis<MonomialBasis>(ff);
        sink(back[0]);
    }
    t.report("basis roundtrip N=10", reps);
}

static void bench_taylor_shift(int reps) {
    Timer t;
    Polynomial<10, uint64_t, MonomialBasis> p{{1,2,3,4,5,6,7,8,9,10}};
    for (int r = 0; r < reps; ++r) {
        auto shifted = taylor_shift(p, uint64_t(r));
        sink(shifted[0]);
    }
    t.report("taylor_shift N=10", reps);
}

static void bench_formal_derivative(int reps) {
    Timer t;
    Polynomial<10, uint64_t, MonomialBasis> p{{1,2,3,4,5,6,7,8,9,10}};
    for (int r = 0; r < reps; ++r) {
        auto d = formal_derivative(p);
        sink(d[0]);
    }
    t.report("formal_derivative N=10", reps);
}

static void bench_forward_difference(int reps) {
    Timer t;
    Polynomial<10, uint64_t, MonomialBasis> p{{1,2,3,4,5,6,7,8,9,10}};
    for (int r = 0; r < reps; ++r) {
        auto d = forward_difference(p);
        sink(d[0]);
    }
    t.report("forward_difference N=10", reps);
}

static void bench_poly_mul(int reps) {
    Timer t;
    Polynomial<10, uint64_t, MonomialBasis> a{{1,2,3,4,5,6,7,8,9,10}};
    Polynomial<10, uint64_t, MonomialBasis> b{{10,9,8,7,6,5,4,3,2,1}};
    for (int r = 0; r < reps; ++r) {
        auto prod = a * b;
        sink(prod[0]);
    }
    t.report("poly_mul N=10x10", reps);
}

static void bench_witt_add(int reps) {
    Timer t;
    WittVector<5, uint64_t> a{{1,2,3,4,5}};
    WittVector<5, uint64_t> b{{5,4,3,2,1}};
    for (int r = 0; r < reps; ++r) {
        auto sum = witt_add(a, b);
        sink(sum[0]);
    }
    t.report("witt_add N=5", reps);
}

static void bench_witt_mul(int reps) {
    Timer t;
    WittVector<5, uint64_t> a{{1,2,3,4,5}};
    WittVector<5, uint64_t> b{{5,4,3,2,1}};
    for (int r = 0; r < reps; ++r) {
        auto prod = witt_mul(a, b);
        sink(prod[0]);
    }
    t.report("witt_mul N=5", reps);
}

static void bench_compose(int reps) {
    Timer t;
    Polynomial<5, uint64_t, MonomialBasis> P{{0,1,1,0,0}};
    Polynomial<5, uint64_t, MonomialBasis> Q{{0,1,1,0,0}};
    for (int r = 0; r < reps; ++r) {
        auto C = compose(P, Q);
        sink(C[0]);
    }
    t.report("compose N=5x5", reps);
}

static void bench_reversion(int reps) {
    Timer t;
    Polynomial<5, uint64_t, MonomialBasis> P{{0,1,1,0,0}};
    for (int r = 0; r < reps; ++r) {
        auto R = reversion(P);
        sink(R[0]);
    }
    t.report("reversion N=5", reps);
}

static void bench_adams(int reps) {
    Timer t;
    WittVector<4, uint32_t> w{{1,2,3,4}};
    for (int r = 0; r < reps; ++r) {
        auto psi = adams_operation(w, 3);
        sink(psi[0]);
    }
    t.report("adams_operation N=4 n=3", reps);
}

template<int N>
static void bench_modinv_pow2(int reps) {
    Timer t;
    Polynomial<N, uint64_t, MonomialBasis> a;
    a[0] = 0xDEADBEEFCAFEBABEULL | 1;
    for (int i = 1; i < N; ++i) a[i] = 0xCAFEBABEDEADBEEFULL + i;
    for (int r = 0; r < reps; ++r) {
        auto inv = modinv_pow2(a);
        sink(inv[0]);
    }
    char buf[64];
    std::sprintf(buf, "modinv_pow2 N=%-2d (%d bits)", N, N * 64);
    t.report(buf, reps);
}

static void bench_indefinite_sum(int reps) {
    Timer t;
    Polynomial<10, uint64_t, MonomialBasis> p{{1,2,3,4,5,6,7,8,9,10}};
    for (int r = 0; r < reps; ++r) {
        auto sum = indefinite_sum(p);
        sink(sum[0]);
    }
    t.report("indefinite_sum N=10", reps);
}

int main() {
    std::printf("=== dyadic Benchmarks ===\n\n");

    const int REPS = 50000;

    bench_polynomial_eval(REPS);
    bench_basis_roundtrip(REPS);
    bench_taylor_shift(REPS);
    bench_formal_derivative(REPS);
    bench_forward_difference(REPS);
    bench_poly_mul(REPS);
    bench_witt_add(REPS);
    bench_witt_mul(REPS);
    bench_compose(REPS);
    bench_reversion(REPS);
    bench_adams(REPS);
    bench_indefinite_sum(REPS);

    std::printf("\n--- modinv_pow2 ---\n");
    bench_modinv_pow2<1>(200000);
    bench_modinv_pow2<2>(100000);
    bench_modinv_pow2<4>(50000);
    bench_modinv_pow2<8>(20000);
    bench_modinv_pow2<16>(5000);
    bench_modinv_pow2<32>(2000);
    bench_modinv_pow2<64>(500);
    bench_modinv_pow2<128>(200);

    std::printf("\n=== Done ===\n");
    return 0;
}
