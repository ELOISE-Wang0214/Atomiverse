/// @file test_ewald.cpp
#include <gtest/gtest.h>
#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/core/box.h"
#include "atomiverse/core/vec3.h"

using namespace ATV;
using namespace ATV::cpu::pw::internal;

TEST(Ewald, SiTwoAtoms) {
    Box box = Box::cubic(10.26); // Bohr
    std::vector<Vec3> pos = {{0,0,0}, {2.565, 2.565, 2.565}}; // Bohr
    std::vector<Real> charges = {4, 4};
    Real ecut = Real(15); // Hartree

    EwaldSummation e1(pos, charges, box, ecut);
    EwaldSummation e2(pos, charges, box, ecut*Real(2));
    Real diff = std::abs(e1.total_energy() - e2.total_energy());
    EXPECT_LT(static_cast<double>(diff), 1e-4); // Hartree
}

TEST(Ewald, ZeroCharges) {
    Box box = Box::cubic(10);
    std::vector<Vec3> pos = {{0,0,0}, {2,2,2}};
    std::vector<Real> charges = {0, 0};
    EwaldSummation e(pos, charges, box, Real(15));
    EXPECT_NEAR(static_cast<double>(e.total_energy()), 0.0, 1e-12);
}
