#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
FAIL=0

# Detect available compilers
GCC=""
for c in g++ g++-14 g++-13 g++-12; do
    if command -v "$c" &>/dev/null; then GCC="$c"; break; fi
done

CLANG=""
for c in clang++-20 clang++-19 clang++-18 clang++-17 clang++; do
    if command -v "$c" &>/dev/null; then CLANG="$c"; break; fi
done

red() { printf '\033[31m%s\033[0m\n' "$1"; }
green() { printf '\033[32m%s\033[0m\n' "$1"; }

check() {
    local label="$1" cmd="$2"
    local out
    out=$(bash -c "$cmd" 2>&1) || {
        echo "$out"
        red "  FAIL  $label"
        FAIL=1
        return
    }
    # Only show filtered output (hide known warnings)
    echo "$out" | grep -v 'pragma once in main file' | grep -v 'ISO C++ does not support' || true
    green "  PASS  $label"
}

echo "=== Compile-time Verification (GCC: $GCC, Clang: ${CLANG:-none}) ==="

# --- GCC verify header ---
check "GCC: verify.h standalone (light proofs)" \
    "$GCC -std=c++20 -O2 -I '$DIR' -c -x c++ '$DIR/dyadic_verify.h' -o /dev/null -Wall -Wextra -Wpedantic -fconstexpr-ops-limit=200000000"

check "GCC: verify.h standalone (heavy proofs)" \
    "$GCC -std=c++20 -O2 -I '$DIR' -c -x c++ '$DIR/dyadic_verify.h' -o /dev/null -Wall -Wextra -Wpedantic -DDYADIC_HEAVY_PROOFS -fconstexpr-ops-limit=200000000"

# --- GCC v2.h header ---
check "GCC: v2.h standalone" \
    "$GCC -std=c++20 -O2 -I '$DIR' -c -x c++ '$DIR/dyadic.h' -o /dev/null -Wall -Wextra -Wpedantic"

# --- Clang verify header (exhaustive proofs need extra constexpr steps) ---
if [ -n "$CLANG" ]; then
    check "Clang: verify.h standalone" \
        "$CLANG -std=c++20 -O2 -I '$DIR' -c -x c++ '$DIR/dyadic_verify.h' -o /dev/null -Wall -Wpedantic -fconstexpr-steps=50000000"

    check "Clang: v2.h standalone" \
        "$CLANG -std=c++20 -O2 -I '$DIR' -c -x c++ '$DIR/dyadic.h' -o /dev/null -Wall -Wpedantic"
else
    red "  SKIP  Clang: not installed"
fi

echo ""
echo "=== Full Test Suite ==="

# --- GCC test suite ---
check "GCC: verify + runtime proofs" \
    "$GCC -std=c++20 -O2 -I '$DIR' '$DIR/test/test_verify.cpp' -o /tmp/_ci_verify 2>&1 && /tmp/_ci_verify"

check "GCC: property-based tests" \
    "$GCC -std=c++20 -O2 -I '$DIR' '$DIR/test/test_property.cpp' -o /tmp/_ci_property 2>&1 && /tmp/_ci_property"

check "GCC: full functional tests" \
    "$GCC -std=c++20 -O2 -I '$DIR' '$DIR/test/test_full.cpp' -o /tmp/_ci_full 2>&1 && /tmp/_ci_full"

# --- Clang test suite ---
if [ -n "$CLANG" ]; then
    check "Clang: verify + runtime proofs" \
        "$CLANG -std=c++20 -O2 -I '$DIR' -fconstexpr-steps=50000000 '$DIR/test/test_verify.cpp' -o /tmp/_ci_verify_c 2>&1 && /tmp/_ci_verify_c"

check "Clang: property-based tests" \
    "$CLANG -std=c++20 -O2 -I '$DIR' -fconstexpr-steps=50000000 '$DIR/test/test_property.cpp' -o /tmp/_ci_property_c 2>&1 && /tmp/_ci_property_c"

    check "Clang: full functional tests" \
        "$CLANG -std=c++20 -O2 -I '$DIR' -fconstexpr-steps=50000000 '$DIR/test/test_full.cpp' -o /tmp/_ci_full_c 2>&1 && /tmp/_ci_full_c"
fi

echo ""
echo "=== Sanitized Build (GCC) ==="
check "GCC: ASan + UBSan" \
    "$GCC -std=c++20 -O1 -fsanitize=address,undefined -I '$DIR' '$DIR/test/test_verify.cpp' -o /tmp/_ci_san 2>&1 && /tmp/_ci_san"

echo ""
if [ "$FAIL" -eq 0 ]; then
    green "=== All CI checks passed ==="
else
    red "=== Some CI checks failed ==="
fi
exit $FAIL
