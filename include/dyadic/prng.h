// dyadic/prng.h — XorShift64 PRNG for testing and verification
// Requires: (none — standalone, only needs <cstdint>)

#pragma once

#include <cstdint>

namespace dyadic {
namespace detail {

struct XorShift64 {
    uint64_t state;
    constexpr explicit XorShift64(uint64_t seed) noexcept : state(seed) {}
    constexpr uint64_t next() noexcept {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }
    static constexpr bool is_valid_seed(uint64_t s) noexcept { return s != 0; }
};

} // namespace detail
} // namespace dyadic
