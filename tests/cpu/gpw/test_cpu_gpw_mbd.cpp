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
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"
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

// Product forces must match an independent central FD of the same energy
// (h_product = 1e-4; check with h = 1e-5).
TEST(GpwMbd, ProductForcesMatchIndependentFiniteDifference) {
    auto pos = water_dimer_bohr();
    const double cell[3] = {0.0, 0.0, 0.0};
    const MbdResult r = mbd_rsscs_pbe(kZ, pos, cell, /*periodic=*/false);

    ASSERT_EQ(r.forces.size(), pos.size());
    EXPECT_NEAR(r.energy, kLibMbdIsolated, 1e-12);

    // Independent FD with a tighter step (energy-only re-eval via public API).
    const double h = 1e-5;
    double max_err = 0.0;
    double sum[3] = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < pos.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            const double save = pos[i][d];
            pos[i][d] = save + h;
            const double ep = mbd_rsscs_pbe(kZ, pos, cell, false).energy;
            pos[i][d] = save - h;
            const double em = mbd_rsscs_pbe(kZ, pos, cell, false).energy;
            pos[i][d] = save;
            const double fd = -(ep - em) / (2.0 * h);
            max_err = std::max(max_err, std::abs(fd - r.forces[i][d]));
            sum[d] += r.forces[i][d];
        }
    }
    std::printf("[MBD] max |F_product - F_fd(h=1e-5)| = %.2e Ha/Bohr\n",
                max_err);
    EXPECT_LT(max_err, 1e-6);
    const double sum_norm =
        std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
    std::printf("[MBD] |sum F| = %.2e Ha/Bohr\n", sum_norm);
    EXPECT_LT(sum_norm, 1e-8);
    // Net dispersion force on O0 should be nonzero (dimer attraction along x).
    EXPECT_GT(std::abs(r.forces[0][0]), 1e-6);
}

TEST(GpwMbd, PeriodicForcesFiniteAndTranslationallyInvariant) {
    auto pos = water_dimer_bohr();
    const double L = 12.0 * kAngToBohr;
    const double cell[3] = {L, L, L};
    const int kmesh[3] = {2, 2, 2};
    const MbdResult r =
        mbd_rsscs_pbe(kZ, pos, cell, /*periodic=*/true, 0.83, 6.0, {}, kmesh);
    ASSERT_EQ(r.forces.size(), pos.size());
    EXPECT_NEAR(r.energy, kLibMbdPeriodic222, 1e-10);
    double sum[3] = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < r.forces.size(); ++i) {
        for (int d = 0; d < 3; ++d) {
            EXPECT_TRUE(std::isfinite(r.forces[i][d]));
            sum[d] += r.forces[i][d];
        }
    }
    const double sum_norm =
        std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
    EXPECT_LT(sum_norm, 1e-7);
}

// Multi-k force assembly must fold geometry-only MBD (was skipped before
// gpw-mbd-kpoint-forces). H2 1×1×2: one SCF, two assemblies vs geometry kernel.
TEST(GpwMbd, KpointProductForcesFoldMbd) {
    using ATV::cpu::gpw::GaussianBasis;
    using ATV::cpu::gpw::GpwAtom;
    using ATV::cpu::gpw::GpwDft;
    using ATV::cpu::gpw::GthPotential;
    using ATV::cpu::gpw::XcKind;
    using ATV::cpu::gpw::gpw_forces_from_scf;

    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.5, 6.0, 6.0};
    std::vector<GpwAtom> atoms = {GpwAtom(rA, basis, pot, 1),
                                  GpwAtom(rB, basis, pot, 1)};

    auto make_cfg = [](bool mbd) {
        GpwDft m;
        m.cell[0] = m.cell[1] = m.cell[2] = 12.0;
        m.cutoff_ry = 100.0;
        m.n_occ = 1;
        m.tol = 1e-8;
        m.max_iter = 60;
        m.nonperiodic = false;
        m.local_pp_on_grid = true;
        m.xc = XcKind::PADE;
        m.eps_default = 1e-8;
        m.k_mesh[0] = 1;
        m.k_mesh[1] = 1;
        m.k_mesh[2] = 2;
        m.mbd = mbd;
        return m;
    };

    GpwDft cfg = make_cfg(true);
    cfg.build(atoms);
    const auto& scf = cfg.run(atoms);
    ASSERT_TRUE(scf.converged);
    // 1×1×2 + time-reversal often yields n_kpts==1 irreducible point, but the
    // SCF still uses the multi-k path (C_occ_k payload; Γ has n_kpts==0).
    EXPECT_GE(scf.n_kpts, 1);
    EXPECT_FALSE(scf.C_occ_k.empty());
    EXPECT_TRUE(std::isfinite(scf.e_disp));
    EXPECT_LT(scf.e_disp, 0.0);

    // Force assembly only reads public flags + cell from cfg (GpwDft non-copyable).
    const auto F_mbd = gpw_forces_from_scf(atoms, cfg, scf);
    GpwDft plain = make_cfg(false);
    const auto F_plain = gpw_forces_from_scf(atoms, plain, scf);
    ASSERT_EQ(F_mbd.size(), 2u);
    ASSERT_EQ(F_plain.size(), 2u);

    std::vector<int> Z = {1, 1};
    std::vector<std::array<double, 3>> pos = {
        {rA[0], rA[1], rA[2]}, {rB[0], rB[1], rB[2]}};
    const double cell0[3] = {0.0, 0.0, 0.0}; // free-atom MBD (default mesh)
    const MbdResult geom = mbd_rsscs_pbe(Z, pos, cell0, /*periodic=*/false);

    double max_d = 0.0;
    for (int a = 0; a < 2; ++a) {
        for (int k = 0; k < 3; ++k) {
            const double d =
                (F_mbd[static_cast<std::size_t>(a)][k]
                 - F_plain[static_cast<std::size_t>(a)][k])
                - geom.forces[static_cast<std::size_t>(a)][k];
            max_d = std::max(max_d, std::abs(d));
        }
    }
    std::printf("[MBD k-forces] e_disp=%.6e  max|F_mbd-F_plain-F_geom|=%.2e\n",
                scf.e_disp, max_d);
    EXPECT_LT(max_d, 1e-8)
        << "multi-k assemble must add geometry-only MBD forces";
    // Non-trivial MBD contribution on the compressed H–H bond.
    EXPECT_GT(std::abs(geom.forces[0][0]), 1e-8);
}
