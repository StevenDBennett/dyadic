#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
FAIL=0

# Include paths: root (for umbrella dyadic.h) and include/ (for sub-headers)
INC="-I '$DIR' -I '$DIR/include'"

# Detect available compilers — discover by pattern instead of hardcoded version ranges
find_compiler() {
    local base="$1"
    local found=""
    # Check base name first (might be a versioned wrapper already)
    if command -v "$base" &>/dev/null && "$base" -std=c++20 -x c++ - -o /dev/null </dev/null &>/dev/null; then
        echo "$base"
        return
    fi
    # Search PATH for version-suffixed names, pick the highest version
    for dir in $(echo "$PATH" | tr ':' ' '); do
        [ -d "$dir" ] || continue
        for f in "$dir/$base"-[0-9]*; do
            [ -x "$f" ] || continue
            local name="${f##*/}"
            # Take the first one that works (they're found in PATH order)
            if [ -z "$found" ] && "$name" -std=c++20 -x c++ - -o /dev/null </dev/null &>/dev/null; then
                found="$name"
            fi
        done
    done
    if [ -n "$found" ]; then
        echo "$found"
    else
        # Final fallback: try the base name
        command -v "$base" &>/dev/null && echo "$base" || echo ""
    fi
}

GCC=$(find_compiler g++)
CLANG=$(find_compiler clang++)

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

echo "=== Compile-Time Verification (GCC: $GCC, Clang: ${CLANG:-none}) ==="

# --- GCC verify header ---
check "GCC: verify.h standalone (light proofs)" \
    "$GCC -std=c++20 -O2 $INC -c -x c++ '$DIR/include/dyadic/verify.h' -o /dev/null -Wall -Wextra -Wpedantic -fconstexpr-ops-limit=200000000"

check "GCC: verify.h standalone (heavy proofs)" \
    "$GCC -std=c++20 -O2 $INC -c -x c++ '$DIR/include/dyadic/verify.h' -o /dev/null -Wall -Wextra -Wpedantic -DDYADIC_HEAVY_PROOFS -fconstexpr-ops-limit=200000000"

# --- GCC dyadic.h header ---
check "GCC: dyadic.h standalone" \
    "$GCC -std=c++20 -O2 $INC -c -x c++ '$DIR/dyadic.h' -o /dev/null -Wall -Wextra -Wpedantic"

# --- Clang verify header (exhaustive proofs need extra constexpr steps) ---
if [ -n "$CLANG" ]; then
    check "Clang: verify.h standalone" \
        "$CLANG -std=c++20 -O2 $INC -c -x c++ '$DIR/include/dyadic/verify.h' -o /dev/null -Wall -Wpedantic -fconstexpr-steps=50000000"

    check "Clang: dyadic.h standalone" \
        "$CLANG -std=c++20 -O2 $INC -c -x c++ '$DIR/dyadic.h' -o /dev/null -Wall -Wpedantic"
else
    red "  SKIP  Clang: not installed"
fi

echo ""
echo "=== Full Test Suite ==="

# --- GCC test suite ---
check "GCC: core unit tests" \
    "$GCC -std=c++20 -O2 $INC '$DIR/test/test_core.cpp' -o /tmp/_ci_core 2>&1 && /tmp/_ci_core"

check "GCC: verify + runtime proofs" \
    "$GCC -std=c++20 -O2 $INC '$DIR/test/test_verify.cpp' -o /tmp/_ci_verify 2>&1 && /tmp/_ci_verify"

check "GCC: property-based tests" \
    "$GCC -std=c++20 -O2 $INC '$DIR/test/test_property.cpp' -o /tmp/_ci_property 2>&1 && /tmp/_ci_property"

check "GCC: full functional tests" \
    "$GCC -std=c++20 -O2 $INC '$DIR/test/test_full.cpp' -o /tmp/_ci_full 2>&1 && /tmp/_ci_full"

# --- Clang test suite ---
if [ -n "$CLANG" ]; then
    check "Clang: core unit tests" \
        "$CLANG -std=c++20 -O2 $INC -fconstexpr-steps=50000000 '$DIR/test/test_core.cpp' -o /tmp/_ci_core_c 2>&1 && /tmp/_ci_core_c"

    check "Clang: verify + runtime proofs" \
        "$CLANG -std=c++20 -O2 $INC -fconstexpr-steps=50000000 '$DIR/test/test_verify.cpp' -o /tmp/_ci_verify_c 2>&1 && /tmp/_ci_verify_c"

    check "Clang: property-based tests" \
        "$CLANG -std=c++20 -O2 $INC -fconstexpr-steps=50000000 '$DIR/test/test_property.cpp' -o /tmp/_ci_property_c 2>&1 && /tmp/_ci_property_c"

    check "Clang: full functional tests" \
        "$CLANG -std=c++20 -O2 $INC -fconstexpr-steps=50000000 '$DIR/test/test_full.cpp' -o /tmp/_ci_full_c 2>&1 && /tmp/_ci_full_c"
fi

echo ""
echo "=== Sanitized Build (GCC) ==="
check "GCC: ASan + UBSan" \
    "$GCC -std=c++20 -O1 -fsanitize=address,undefined $INC '$DIR/test/test_verify.cpp' -o /tmp/_ci_san 2>&1 && /tmp/_ci_san"

echo ""
if [ "$FAIL" -eq 0 ]; then
    green "=== All CI checks passed ==="
else
    red "=== Some CI checks failed ==="
fi
exit $FAIL
