cmake_minimum_required(VERSION 3.14)
project(dyadic LANGUAGES CXX VERSION 1.0.0)

# Require C++20, disable compiler extensions
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Header-only interface library for IDE discovery
add_library(dyadic INTERFACE)
target_include_directories(dyadic INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(dyadic INTERFACE cxx_std_20)

# Option for exhaustive compile-time proofs
option(DYADIC_HEAVY_PROOFS "Enable exhaustive (but expensive) compile-time proofs" OFF)
if(DYADIC_HEAVY_PROOFS)
    target_compile_definitions(dyadic INTERFACE DYADIC_HEAVY_PROOFS)
endif()

# Tests (only if not being included as a subdirectory)
option(DYADIC_BUILD_TESTS "Build test executables" ON)
# PROJECT_IS_SUBMODULE requires CMake 3.21+; on older CMake assume not a subproject
if(DYADIC_BUILD_TESTS AND (NOT DEFINED PROJECT_IS_SUBMODULE OR NOT PROJECT_IS_SUBMODULE))
    enable_testing()

    function(add_dyadic_test name source)
        add_executable(${name} ${source})
        target_link_libraries(${name} PRIVATE dyadic)
        # Increase constexpr evaluation budget for compile-time proofs
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            target_compile_options(${name} PRIVATE -fconstexpr-ops-limit=200000000)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(${name} PRIVATE -fconstexpr-steps=50000000)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            target_compile_options(${name} PRIVATE /constexpr:steps50000000)
        endif()
        add_test(NAME ${name} COMMAND ${name})
    endfunction()

    add_dyadic_test(test_core      test/test_core.cpp)
    add_dyadic_test(test_verify    test/test_verify.cpp)
    add_dyadic_test(test_property  test/test_property.cpp)
    add_dyadic_test(test_full      test/test_full.cpp)
    add_dyadic_test(test_negatives test/test_negatives.cpp)  # POSIX-only (fork); Windows stub returns 0

    # Benchmark (not part of default CTest — run manually)
    add_executable(benchmark test/benchmark.cpp)
    target_link_libraries(benchmark PRIVATE dyadic)

    # Sanitized build target (ASan + UBSan) — GCC/Clang only
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        add_executable(test_verify_san test/test_verify.cpp)
        target_link_libraries(test_verify_san PRIVATE dyadic)
        target_compile_options(test_verify_san PRIVATE
            -fsanitize=address,undefined -O1)
        target_link_options(test_verify_san PRIVATE -fsanitize=address,undefined)
        if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            target_compile_options(test_verify_san PRIVATE
                -fconstexpr-ops-limit=200000000 -fconstexpr-depth=1000000)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(test_verify_san PRIVATE
                -fconstexpr-steps=50000000)
        endif()
    endif()
endif()

# Install header-only library
include(GNUInstallDirs)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.h")
install(TARGETS dyadic EXPORT dyadic-targets)
install(EXPORT dyadic-targets
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/dyadic)
