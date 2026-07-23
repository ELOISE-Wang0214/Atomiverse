/// @file tests/cpu/gpw/test_cpu_gpw_dftd4.cpp
/// @brief Periodic + isolated DFT-D4 E_disp vs offline dftd4 3.7.0 pins.

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/scf/dftd4.h"

using ATV::cpu::scf::dftd4_energy;

namespace {

constexpr double kAngToBohr = 1.0 / 0.52917720859;

// Water dimer (Å) matching gpw-pbe-d3bj-ice / gpw-d4-water REFERENCE.
std::vector<std::array<double, 3>> water_dimer_bohr() {
    const double ang[6][3] = {
        {0.000000, 0.000000, 0.000000},  // O
        {0.758602, 0.000000, 0.504284},  // H
        {0.758602, 0.000000, -0.504284}, // H
        {3.000000, 0.000000, 0.000000},  // O
        {3.758602, 0.000000, 0.504284},  // H
        {3.758602, 0.000000, -0.504284}, // H
    };
    std::vector<std::array<double, 3>> pos(6);
    for (int i = 0; i < 6; ++i)
        for (int d = 0; d < 3; ++d)
            pos[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] =
                ang[i][d] * kAngToBohr;
    return pos;
}

const std::vector<int> kZ = {8, 1, 1, 8, 1, 1};

} // namespace

// Isolated two-body: dftd4 3.7.0 --mbdscale 0.0 on XYZ.
TEST(GpwDftD4, WaterDimerIsolatedMatchesDftd4) {
    const auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const double e = dftd4_energy(kZ, pos, cell, /*periodic=*/false);
    const double pin = -9.8816586433966e-04;
    std::printf("[D4 iso] E=%.14e  pin=%.14e  diff=%.2e Ha\n", e, pin, e - pin);
    EXPECT_NEAR(e, pin, 1e-8);
}

// Periodic cubic 12 Å two-body: dftd4 3.7.0 --mbdscale 0.0 on POSCAR.
TEST(GpwDftD4, WaterDimerPeriodicMatchesDftd4) {
    const auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    const double e = dftd4_energy(kZ, pos, cell, /*periodic=*/true);
    const double pin = -9.9458842745847e-04;
    std::printf("[D4 per] E=%.14e  pin=%.14e  diff=%.2e Ha\n", e, pin, e - pin);
    EXPECT_NEAR(e, pin, 1e-8);
    // Images add attraction vs isolated.
    const double e_iso = dftd4_energy(kZ, pos, cell, /*periodic=*/false);
    EXPECT_LT(e, e_iso);
}
