// test_util.h — Shared test utilities
#pragma once
#include <cstdlib>
#include <cstdint>
#include <string_view>

namespace dyadic_test {

inline uint64_t parse_seed(int argc, char** argv, uint64_t default_seed) {
    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);
        if (arg.starts_with("--seed=")) {
            return static_cast<uint64_t>(std::atoll(arg.data() + 7));
        } else if (arg == "--seed") {
            if (i + 1 < argc) return static_cast<uint64_t>(std::atoll(argv[++i]));
        }
    }
    return default_seed;
}

} // namespace dyadic_test
