/// @file tests/cpu/gpw/test_gpw_method.cpp
#include <gtest/gtest.h>
#include "atomiverse/cpu/gpw/gpw_solver.h"
using namespace ATV::cpu::gpw;
TEST(GpwMethod, DefaultMembersConstructible) {
    GpwDft m;
    EXPECT_GT(m.cutoff_ry, 0.0);
    EXPECT_EQ(m.ngrids, 1);
    EXPECT_FALSE(m.use_ot);
}
