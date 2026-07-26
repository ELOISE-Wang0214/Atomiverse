/// @file tests/cpu/gpw/test_cpu_gpw_dftd4.cpp
/// @brief Periodic + isolated DFT-D4 E_disp and forces vs offline dftd4 3.7.0.

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/scf/dftd4.h"

using ATV::cpu::scf::D4Result;
using ATV::cpu::scf::dftd4;
using ATV::cpu::scf::dftd4_energy;

namespace {

constexpr double kAngToBohr = 1.0 / 0.52917720859;

// Water dimer (Å) matching gpw-pbe-d3bj-ice / gpw-d4-water REFERENCE.
// Atom order: O H H O H H (matches isolated XYZ pin, not POSCAR O O H H H H).
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

// Offline dftd4 3.7.0 -g on XYZ (order O H H O H H): gradient dE/dR [Ha/Bohr].
// ATV force = −gradient.
const double kGradIso[6][3] = {
    {-5.367487773544261e-05, 0.0, 6.468403283702467e-20},
    {-5.984939684718166e-05, 0.0, -1.0621704377196226e-05},
    {-5.9849396847181536e-05, 0.0, 1.0621704377196153e-05},
    {1.371416991978683e-04, 0.0, 4.748814682847428e-21},
    {1.811598611596879e-05, 0.0, -8.567147309218983e-06},
    {1.8115986115968777e-05, 0.0, 8.567147309218974e-06},
};

// Offline dftd4 3.7.0 -g on POSCAR (order O O H H H H) remapped to O H H O H H.
// POSCAR atoms: O0, O1, H0, H1, H2, H3 → ATV [O0, H0, H1, O1, H2, H3].
const double kGradPer[6][3] = {
    {-5.3738961271913915e-05, -1.5591535100726387e-24, -6.1987516039752966e-21},
    {-5.9550148958160894e-05, -3.0412121452621079e-25, -1.0810592129949165e-05},
    {-5.9550148958160921e-05, -5.7104654456581212e-25, 1.0810592129949182e-05},
    {1.3609823922940861e-04, 6.8060330704498828e-25, -8.4063699110321910e-20},
    {1.8370509979413634e-05, 5.7191630572539285e-25, -8.5497124425277334e-06},
    {1.8370509979413539e-05, 1.2193640572572847e-24, 8.5497124425278096e-06},
};

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

TEST(GpwDftD4, IsolatedForcesMatchDftd4Grad) {
    const auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const D4Result r = dftd4(kZ, pos, cell, /*periodic=*/false);
    EXPECT_NEAR(r.energy, -9.8816586433966e-04, 1e-8);
    double max_err = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 3; ++k) {
            // Force = −gradient
            const double want = -kGradIso[i][k];
            const double got = r.forces[static_cast<std::size_t>(i)]
                                       [static_cast<std::size_t>(k)];
            max_err = std::max(max_err, std::abs(got - want));
        }
    }
    std::printf("[D4 iso F] max|F+grad|=%.3e Ha/Bohr\n", max_err);
    EXPECT_LT(max_err, 1e-7);
}

TEST(GpwDftD4, PeriodicForcesMatchDftd4Grad) {
    const auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    const D4Result r = dftd4(kZ, pos, cell, /*periodic=*/true);
    EXPECT_NEAR(r.energy, -9.9458842745847e-04, 1e-8);
    double max_err = 0.0;
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 3; ++k) {
            const double want = -kGradPer[i][k];
            const double got = r.forces[static_cast<std::size_t>(i)]
                                       [static_cast<std::size_t>(k)];
            max_err = std::max(max_err, std::abs(got - want));
        }
    }
    std::printf("[D4 per F] max|F+grad|=%.3e Ha/Bohr\n", max_err);
    EXPECT_LT(max_err, 1e-7);
}

TEST(GpwDftD4, ForcesMatchEnergyCentralDifference) {
    const auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    const D4Result r = dftd4(kZ, pos, cell, /*periodic=*/true);
    // Independent FD with a different step confirms internal consistency.
    constexpr double h = 2.0e-4;
    double max_err = 0.0;
    auto p = pos;
    for (int i = 0; i < 6; ++i) {
        for (int k = 0; k < 3; ++k) {
            const double x0 = p[static_cast<std::size_t>(i)]
                               [static_cast<std::size_t>(k)];
            p[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                x0 + h;
            const double ep = dftd4_energy(kZ, p, cell, true);
            p[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
                x0 - h;
            const double em = dftd4_energy(kZ, p, cell, true);
            p[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] = x0;
            const double f_fd = -(ep - em) / (2.0 * h);
            max_err = std::max(
                max_err, std::abs(r.forces[static_cast<std::size_t>(i)]
                                          [static_cast<std::size_t>(k)]
                                  - f_fd));
        }
    }
    std::printf("[D4 self-FD] max|F-Ffd|=%.3e Ha/Bohr\n", max_err);
    EXPECT_LT(max_err, 1e-7);
}

// Isolated ATM: dftd4 3.7.0 -f pbe --mbdscale 1.0 on XYZ.
TEST(GpwDftD4, WaterDimerIsolatedAtmMatchesDftd4) {
    const auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const double e = dftd4_energy(kZ, pos, cell, /*periodic=*/false,
                                  /*s6=*/1.0, /*s8=*/0.95948085,
                                  /*a1=*/0.38574991, /*a2=*/4.80688534,
                                  /*ga=*/3.0, /*gc=*/2.0, /*wf=*/6.0,
                                  /*charge=*/0.0, /*cn_cut=*/30.0,
                                  /*disp2_cut=*/60.0, /*s9=*/1.0);
    const double pin = -9.8819741533958208e-04;
    const double e0 = dftd4_energy(kZ, pos, cell, /*periodic=*/false);
    std::printf("[D4 iso ATM] E=%.14e  pin=%.14e  dE_atm=%.2e Ha\n", e, pin,
                e - e0);
    EXPECT_NEAR(e, pin, 1e-8);
    // ATM contribution must be nonzero on this dimer.
    EXPECT_GT(std::abs(e - e0), 1e-10);
}

// Periodic ATM: dftd4 3.7.0 -f pbe (default s9=1) on POSCAR.
TEST(GpwDftD4, WaterDimerPeriodicAtmMatchesDftd4) {
    const auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    const double e = dftd4_energy(kZ, pos, cell, /*periodic=*/true,
                                  /*s6=*/1.0, /*s8=*/0.95948085,
                                  /*a1=*/0.38574991, /*a2=*/4.80688534,
                                  /*ga=*/3.0, /*gc=*/2.0, /*wf=*/6.0,
                                  /*charge=*/0.0, /*cn_cut=*/30.0,
                                  /*disp2_cut=*/60.0, /*s9=*/1.0);
    const double pin = -9.9554475834240288e-04;
    const double e0 = dftd4_energy(kZ, pos, cell, /*periodic=*/true);
    std::printf("[D4 per ATM] E=%.14e  pin=%.14e  dE_atm=%.2e Ha\n", e, pin,
                e - e0);
    EXPECT_NEAR(e, pin, 1e-8);
    EXPECT_GT(std::abs(e - e0), 1e-10);
}
