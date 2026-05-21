/// @file test_g_vector.cpp
#include "atomiverse/core/box.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/vec3.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include <gtest/gtest.h>

using namespace ATV;
using namespace ATV::cpu::pw::internal;

TEST(GVector, SiFCC_30Ry_GCount) {
    Real a = Real(5.43) * kAngstrom2Bohr; // 10.26 Bohr
    auto box = Box::cubic(a);
    Real ecut = Real(15); // 30 Ry = 15 Hartree
    auto gv = enumerate_g_vectors(box, ecut, Vec3::zeros());
    EXPECT_GE(gv.count(), 1000); // ~3000 for cubic 10.26Bohr cell at 15Ht
    // All G vectors satisfy |G|^2/2 <= ecut
    for (int i = 0; i < gv.count(); ++i)
        EXPECT_LE(gv.ekin[i], ecut + Real(1e-10));
}

TEST(GVector, ZeroEcut) {
    auto box = Box::cubic(Real(10));
    auto gv = enumerate_g_vectors(box, Real(0), Vec3::zeros());
    EXPECT_EQ(gv.count(), 1); // only G=0
}

TEST(GVector, SortedByNorm) {
    auto box = Box::cubic(Real(10));
    auto gv = enumerate_g_vectors(box, Real(2), Vec3::zeros());
    for (int i = 1; i < gv.count(); ++i)
        EXPECT_LE(gv.g[i - 1].norm_squared(),
                  gv.g[i].norm_squared() + Real(1e-12));
}
