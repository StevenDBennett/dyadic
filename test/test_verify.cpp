// test_verify.cpp — Compile-time proofs + runtime verification
// Including the header triggers all static_assert compile-time proofs.
// The runtime verifications extend coverage beyond what constexpr can handle.
// Usage: test_verify [--seed <uint64>]

#include <dyadic/verify.h>
#include <cstdio>
#include <cstdlib>

using namespace dyadic;
using namespace dyadic::verify;

int main(int argc, char** argv) {
    uint64_t base_seed = 0xDEADBEEF;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-' && argv[i][1] == '-' && argv[i][2] == 's' &&
            argv[i][3] == 'e' && argv[i][4] == 'e' && argv[i][5] == 'd' &&
            argv[i][6] == '=') {
            base_seed = static_cast<uint64_t>(std::atoll(argv[i] + 7));
        } else if (argv[i][0] == '-' && argv[i][1] == '-' && argv[i][2] == 's' &&
                   argv[i][3] == 'e' && argv[i][4] == 'e' && argv[i][5] == 'd') {
            if (i + 1 < argc) base_seed = static_cast<uint64_t>(std::atoll(argv[++i]));
        }
    }

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

    std::printf("=== Runtime Verification Suite (seed=0x%llx) ===\n", (unsigned long long)base_seed);

    run("runtime N=2 uint16_t", run_all_verifications<2, uint16_t>(base_seed));
    run("runtime N=3 uint16_t", run_all_verifications<3, uint16_t>(base_seed));
    run("runtime N=4 uint16_t", run_all_verifications<4, uint16_t>(base_seed));
    run("runtime N=2 uint32_t", run_all_verifications<2, uint32_t>(base_seed));
    run("runtime N=3 uint32_t", run_all_verifications<3, uint32_t>(base_seed));
    run("runtime N=4 uint32_t", run_all_verifications<4, uint32_t>(base_seed));
    run("runtime N=2 uint64_t", run_all_verifications<2, uint64_t>(base_seed));
    run("runtime N=3 uint64_t", run_all_verifications<3, uint64_t>(base_seed));
    run("runtime N=4 uint64_t", run_all_verifications<4, uint64_t>(base_seed));

    if (failures == 0) {
        std::printf("=== All %d runtime checks passed ===\n", total);
    } else {
        std::printf("=== %d runtime checks failed out of %d ===\n", failures, total);
    }

    return failures;
}
