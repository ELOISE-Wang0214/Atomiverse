/// @file test_units.cpp
/// @brief Unit tests for new unit conversion constants.

#include <gtest/gtest.h>

#include "atomiverse/core/real.h"
#include "atomiverse/core/units.h"

using namespace ATV;

TEST(Units, AngstromToBohr) {
    Real bohr = Real(5.43) * kAngstrom2Bohr;
    EXPECT_NEAR(static_cast<double>(bohr),
                5.43 * static_cast<double>(kAngstrom2Bohr), 1e-2);
}

TEST(Units, RyToHartreeExact) { EXPECT_EQ(kRy2Hartree, Real(0.5)); }

TEST(Units, RyToEV) { EXPECT_EQ(kRy2eV, Real(13.605693)); }

TEST(Units, RoundtripAngstromBohr) {
    Real a = 5.43;
    Real b = a * kAngstrom2Bohr;
    Real back = b * kBohr2Ang;
    EXPECT_NEAR(static_cast<double>(back), static_cast<double>(a), 1e-6);
}
