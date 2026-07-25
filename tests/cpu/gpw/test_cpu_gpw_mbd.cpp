/// @file tests/cpu/gpw/test_cpu_gpw_mbd.cpp
/// @brief Validate free-atom and periodic MBD@rsSCS vs offline libMBD pins.

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/mbd_dispersion.h"
#include "atomiverse/molrs.hpp"

using ATV::cpu::gpw::mbd_rsscs_pbe;
using ATV::cpu::gpw::MbdResult;

namespace {

constexpr double kAngToBohr = 1.0 / 0.52917720859;

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
            pos[i][d] = ang[i][d] * kAngToBohr;
    return pos;
}

const std::vector<int> kZ = {8, 1, 1, 8, 1, 1};

// Offline libMBD Fortran cmbd_mbd_scs_energy (git 88069c6, ATV Å→Bohr).
constexpr double kLibMbdIsolated = -1.4928804171407961e-03;

// Offline pure-Python pymbd.mbd_energy, ATV Å→Bohr, lattice 12 Å cubic,
// k_grid=(2,2,2), β=0.83 free-atom.
constexpr double kLibMbdPeriodic222 = -1.4965720704989138e-03;

} // namespace

TEST(GpwMbd, WaterDimerIsolatedMatchesLibMbd) {
    const auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const MbdResult r = mbd_rsscs_pbe(kZ, pos, cell, /*periodic=*/false);

    std::printf("[MBD] isolated E=%.16e  libMBD=%.16e  diff=%.2e Ha\n", r.energy,
                kLibMbdIsolated, r.energy - kLibMbdIsolated);
    EXPECT_NEAR(r.energy, kLibMbdIsolated, 1e-12);
    EXPECT_LT(r.energy, 0.0);
}

TEST(GpwMbd, WaterDimerPeriodicEwaldMatchesPymbd) {
    const auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    const int kmesh[3] = {2, 2, 2};
    const MbdResult r =
        mbd_rsscs_pbe(kZ, pos, cell, /*periodic=*/true, 0.83, 6.0, {}, kmesh);

    std::printf("[MBD] periodic 2x2x2 E=%.16e  pymbd=%.16e  diff=%.2e Ha\n",
                r.energy, kLibMbdPeriodic222, r.energy - kLibMbdPeriodic222);
    EXPECT_NEAR(r.energy, kLibMbdPeriodic222, 1e-10);
    EXPECT_LT(r.energy, 0.0);
    // Periodic correction is small but nonzero vs isolated on this box.
    EXPECT_GT(std::abs(r.energy - kLibMbdIsolated), 1e-8);
}

TEST(GpwMbd, GpwDftReportsEDispIsolated) {
    using ATV::cpu::gpw::GaussianBasis;
    using ATV::cpu::gpw::GpwAtom;
    using ATV::cpu::gpw::GpwDft;
    using ATV::cpu::gpw::GthPotential;
    using ATV::cpu::gpw::XcKind;

    const double ang[6][3] = {
        {0.000000, 0.000000, 0.000000}, {0.758602, 0.000000, 0.504284},
        {0.758602, 0.000000, -0.504284}, {3.000000, 0.000000, 0.000000},
        {3.758602, 0.000000, 0.504284}, {3.758602, 0.000000, -0.504284},
    };
    const GaussianBasis bO(molrs::Element::O, "DZVP-MOLOPT-SR-GTH");
    const GaussianBasis bH(molrs::Element::H, "DZVP-MOLOPT-SR-GTH");
    const GthPotential pO(molrs::Element::O, "GTH-PBE-q6");
    const GthPotential pH(molrs::Element::H, "GTH-PBE-q1");
    std::vector<GpwAtom> atoms;
    for (int i = 0; i < 6; ++i) {
        const double pos[3] = {ang[i][0] * kAngToBohr, ang[i][1] * kAngToBohr,
                               ang[i][2] * kAngToBohr};
        if (i == 0 || i == 3)
            atoms.emplace_back(pos, bO, pO, 8);
        else
            atoms.emplace_back(pos, bH, pH, 1);
    }

    GpwDft gpw;
    const double L = 12.0 * kAngToBohr;
    gpw.cell[0] = gpw.cell[1] = gpw.cell[2] = L;
    gpw.cutoff_ry = 200.0;
    gpw.n_occ = 8;
    gpw.max_iter = 6;
    gpw.tol = 1e-5;
    gpw.mix = 0.4;
    gpw.xc = XcKind::PBE;
    gpw.mbd = true;
    // default mbd_k_mesh={0,0,0} → free-atom isolated

    gpw.build(atoms);
    const auto& res = gpw.run(atoms);
    EXPECT_NEAR(res.e_disp, kLibMbdIsolated, 1e-12);
    EXPECT_LT(res.e_disp, 0.0);
    EXPECT_TRUE(std::isfinite(res.e_total));
}

TEST(GpwMbd, GpwDftReportsEDispPeriodic) {
    using ATV::cpu::gpw::GaussianBasis;
    using ATV::cpu::gpw::GpwAtom;
    using ATV::cpu::gpw::GpwDft;
    using ATV::cpu::gpw::GthPotential;
    using ATV::cpu::gpw::XcKind;

    const double ang[6][3] = {
        {0.000000, 0.000000, 0.000000}, {0.758602, 0.000000, 0.504284},
        {0.758602, 0.000000, -0.504284}, {3.000000, 0.000000, 0.000000},
        {3.758602, 0.000000, 0.504284}, {3.758602, 0.000000, -0.504284},
    };
    const GaussianBasis bO(molrs::Element::O, "DZVP-MOLOPT-SR-GTH");
    const GaussianBasis bH(molrs::Element::H, "DZVP-MOLOPT-SR-GTH");
    const GthPotential pO(molrs::Element::O, "GTH-PBE-q6");
    const GthPotential pH(molrs::Element::H, "GTH-PBE-q1");
    std::vector<GpwAtom> atoms;
    for (int i = 0; i < 6; ++i) {
        const double pos[3] = {ang[i][0] * kAngToBohr, ang[i][1] * kAngToBohr,
                               ang[i][2] * kAngToBohr};
        if (i == 0 || i == 3)
            atoms.emplace_back(pos, bO, pO, 8);
        else
            atoms.emplace_back(pos, bH, pH, 1);
    }

    GpwDft gpw;
    const double L = 12.0 * kAngToBohr;
    gpw.cell[0] = gpw.cell[1] = gpw.cell[2] = L;
    gpw.cutoff_ry = 200.0;
    gpw.n_occ = 8;
    gpw.max_iter = 6;
    gpw.tol = 1e-5;
    gpw.mix = 0.4;
    gpw.xc = XcKind::PBE;
    gpw.mbd = true;
    gpw.mbd_k_mesh[0] = gpw.mbd_k_mesh[1] = gpw.mbd_k_mesh[2] = 2;

    gpw.build(atoms);
    const auto& res = gpw.run(atoms);
    EXPECT_NEAR(res.e_disp, kLibMbdPeriodic222, 1e-10);
    EXPECT_LT(res.e_disp, 0.0);
    EXPECT_TRUE(std::isfinite(res.e_total));
}
