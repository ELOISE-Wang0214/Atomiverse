/// @file test_pw_si.cpp
/// @brief Integration test for pw::solve_si().

#include <gtest/gtest.h>
#include <cmath>
#include "atomiverse/cpu/pw/solve_si.h"

using namespace ATV::cpu::pw;

TEST(PWSi, SolveSiCompilesAndRuns) {
    // Run with toy potential (no UPF file), low E_cut for speed
    double e = solve_si("", 10.0, 2, 1e-6, 30);
    // Just verify it returns a finite number
    EXPECT_TRUE(std::isfinite(e)) << "Total energy should be finite, got: " << e;
}

TEST(PWSi, SolveSiRunsWithoutCrash) {
    // With toy potential and no UPF file, energy is not physically meaningful.
    // This test just verifies the integration plumbing works end-to-end.
    double e = solve_si("", 10.0, 2, 1e-6, 20);
    EXPECT_TRUE(std::isfinite(e)) << "Energy should be finite, got: " << e;
}
