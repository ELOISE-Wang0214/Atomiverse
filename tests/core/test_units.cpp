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

TEST(Units, ThermalDeBroglieHydrogenAt300K) {
    // Λ = h/sqrt(2π m kB T): H (1.008 amu) at 300 K → 1.0037 Å
    // (h = 4.135667696e-15 eV·s, kB = 8.617333262e-5 eV/K, CODATA 2018).
    const double lambda = units::lambda_thermal(1.008, 300.0);
    EXPECT_NEAR(lambda, 1.0038, 2e-3);
    // Heavier + hotter shrinks it: Λ ~ 1/sqrt(mT).
    EXPECT_NEAR(units::lambda_thermal(4.032, 300.0), lambda / 2.0, 1e-12);
    EXPECT_NEAR(units::lambda_thermal(1.008, 1200.0), lambda / 2.0, 1e-12);
}
