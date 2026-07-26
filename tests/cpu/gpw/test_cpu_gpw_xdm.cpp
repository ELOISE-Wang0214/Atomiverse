/// @file tests/cpu/gpw/test_cpu_gpw_xdm.cpp
/// @brief Validate free-atom XDM vs offline dual-Python pins and FD forces.

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/xdm_dispersion.h"
#include "atomiverse/molrs.hpp"

using ATV::cpu::gpw::xdm_pbe;
using ATV::cpu::gpw::XdmResult;

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

} // namespace

// Isolated free-atom XDM vs offline dual pure-Python lattice sum
// (tools/cp2k_oracle/inputs/gpw-xdm-ice/REFERENCE pin).
TEST(GpwXdm, WaterDimerIsolatedMatchesOfflinePin) {
    const auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const XdmResult r = xdm_pbe(kZ, pos, cell, /*periodic=*/false);

    const double pin = -2.1020435998516220e-03;
    std::printf("[XDM] isolated E=%.16e  pin=%.16e  diff=%.2e Ha\n", r.energy,
                pin, r.energy - pin);
    EXPECT_NEAR(r.energy, pin, 1e-12);
}

// Periodic pure real-space rcut = 50 Å pin (same dual-Python lattice sum).
TEST(GpwXdm, WaterDimerPeriodicEnergy) {
    const auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    constexpr double kRcut = 50.0 * kAngToBohr;
    const XdmResult r = xdm_pbe(kZ, pos, cell, /*periodic=*/true, kRcut);

    const double pin = -2.1148157612879029e-03;
    std::printf("[XDM] periodic E=%.16e  pin=%.16e  diff=%.2e Ha\n", r.energy,
                pin, r.energy - pin);
    EXPECT_NEAR(r.energy, pin, 1e-12);
}

TEST(GpwXdm, AnalyticForcesMatchFiniteDifference) {
    auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const XdmResult r = xdm_pbe(kZ, pos, cell, /*periodic=*/false);

    // Offline dual-Python force sample (REFERENCE gpw-xdm-forces).
    constexpr double kFIsoPin[6][3] = {
        {1.6875829936003448e-04, 0.0, 0.0},
        {6.9736686381575340e-05, 0.0, -1.3172272218329521e-05},
        {6.9736686381575340e-05, 0.0, 1.3172272218329521e-05},
        {-1.1378113807647229e-04, 0.0, 0.0},
        {-9.7225267023356434e-05, 0.0, -1.5448899346549854e-05},
        {-9.7225267023356434e-05, 0.0, 1.5448899346549854e-05},
    };
    double max_pin = 0.0;
    double sum[3] = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < 6; ++i) {
        for (int d = 0; d < 3; ++d) {
            max_pin =
                std::max(max_pin, std::abs(r.forces[i][d] - kFIsoPin[i][d]));
            sum[d] += r.forces[i][d];
        }
    }
    EXPECT_LT(max_pin, 1e-12);
    EXPECT_LT(std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]),
              1e-12);

    const double h = 1e-5;
    double max_err = 0.0;
    for (std::size_t i = 0; i < pos.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            const double save = pos[i][d];
            pos[i][d] = save + h;
            const double ep = xdm_pbe(kZ, pos, cell, false).energy;
            pos[i][d] = save - h;
            const double em = xdm_pbe(kZ, pos, cell, false).energy;
            pos[i][d] = save;
            const double fd = -(ep - em) / (2.0 * h);
            max_err = std::max(max_err, std::abs(fd - r.forces[i][d]));
        }
    }
    std::printf("[XDM] max |F_analytic - F_fd| = %.2e Ha/Bohr (pin %.2e)\n",
                max_err, max_pin);
    EXPECT_LT(max_err, 1e-7);
}

TEST(GpwXdm, PeriodicAnalyticForcesMatchFiniteDifference) {
    auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    constexpr double kRcut = 50.0 * kAngToBohr;
    const XdmResult r = xdm_pbe(kZ, pos, cell, /*periodic=*/true, kRcut);

    const double h = 1e-5;
    double max_err = 0.0;
    double sum[3] = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < pos.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            const double save = pos[i][d];
            pos[i][d] = save + h;
            const double ep =
                xdm_pbe(kZ, pos, cell, true, kRcut).energy;
            pos[i][d] = save - h;
            const double em =
                xdm_pbe(kZ, pos, cell, true, kRcut).energy;
            pos[i][d] = save;
            const double fd = -(ep - em) / (2.0 * h);
            max_err = std::max(max_err, std::abs(fd - r.forces[i][d]));
            sum[d] += r.forces[i][d];
        }
    }
    const double sum_norm =
        std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
    std::printf("[XDM per] max |F_analytic - F_fd| = %.2e  |sum F|=%.2e\n",
                max_err, sum_norm);
    EXPECT_LT(max_err, 1e-7);
    EXPECT_LT(sum_norm, 1e-10);
}

TEST(GpwXdm, GpwDftReportsEDisp) {
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

    constexpr double kEDispRef = -2.1148157612879029e-03;
    GpwDft gpw;
    const double L = 12.0 * kAngToBohr;
    gpw.cell[0] = gpw.cell[1] = gpw.cell[2] = L;
    gpw.cutoff_ry = 200.0;
    gpw.n_occ = 8;
    gpw.max_iter = 8;
    gpw.tol = 1e-6;
    gpw.mix = 0.4;
    gpw.xc = XcKind::PBE;
    gpw.nonperiodic = false;
    gpw.xdm = true;

    gpw.build(atoms);
    const auto& res = gpw.run(atoms);
    EXPECT_NEAR(res.e_disp, kEDispRef, 1e-10);
    EXPECT_TRUE(std::isfinite(res.e_total));
    EXPECT_LT(res.e_disp, 0.0);
}
