/// @file test_crystal.cpp
#include <gtest/gtest.h>
#include "atomiverse/cpu/pw/internal/crystal.h"
#include <cmath>

using namespace ATV;
using namespace ATV::cpu::pw::internal;

TEST(Crystal, SiDiamond) {
    auto f = build_si_diamond(5.43);
    EXPECT_EQ(f.n_atoms(), 2);
    EXPECT_EQ(f.type_id[0], 14);
    EXPECT_EQ(f.type_id[1], 14);
    EXPECT_NEAR(static_cast<double>(f.mass[0]), 28.0855, 1e-4);
}

TEST(Crystal, Positions) {
    auto f = build_si_diamond(5.43);
    Real a4 = Real(5.43) / Real(4);
    EXPECT_NEAR(static_cast<double>(f.x[0]), 0.0, 1e-12);
    EXPECT_NEAR(static_cast<double>(f.y[0]), 0.0, 1e-12);
    EXPECT_NEAR(static_cast<double>(f.z[0]), 0.0, 1e-12);
    EXPECT_NEAR(static_cast<double>(f.x[1]), static_cast<double>(a4), 1e-12);
    EXPECT_NEAR(static_cast<double>(f.y[1]), static_cast<double>(a4), 1e-12);
    EXPECT_NEAR(static_cast<double>(f.z[1]), static_cast<double>(a4), 1e-12);
}
