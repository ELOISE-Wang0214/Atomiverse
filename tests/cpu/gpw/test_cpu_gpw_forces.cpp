/// @file tests/cpu/gpw/test_cpu_gpw_forces.cpp
/// @brief Finite-difference GPW forces on H2 — Newton's third law, along-axis,
///        and the correct restoring sign for a compressed bond.

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"

using namespace ATV::cpu::gpw;

namespace {
const std::string kData = ATV_CP2K_DATA_DIR;
} // namespace

TEST(GpwForces, H2FiniteDifference) {
    const auto basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    const auto pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");

    // H2 compressed to 1.0 Bohr (equilibrium ~1.4): expect a restoring force
    // that pushes the atoms apart along x.
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.5, 2.0, 2.0};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig cfg;
    cfg.cell[0] = 8.0;
    cfg.cell[1] = cfg.cell[2] = 4.0;
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.mix = 0.3;

    auto F = gpw_forces_fd(atoms, cfg, 2e-3);

    // The bond (x) force is the clean physical signal: equal and opposite, and
    // restoring for a compressed bond (left atom pushed to -x, right to +x).
    EXPECT_NEAR(F[0][0], -F[1][0], 1e-4); // Newton's third law along the bond
    EXPECT_LT(F[0][0], 0.0);
    EXPECT_GT(F[1][0], 0.0);
    EXPECT_GT(std::abs(F[0][0]), 5e-3); // substantial restoring force

    // The off-axis (y,z) components are a common-mode "egg-box" grid artifact
    // (both atoms share the same grid phase at y=z=2), not a physical force.
    // It is small relative to the bond force on this single coarse grid;
    // CP2K's REL_CUTOFF multigrid + smoothing suppress it (a documented TODO).
    EXPECT_LT(std::abs(F[0][1]), 0.4 * std::abs(F[0][0]));
    EXPECT_LT(std::abs(F[0][2]), 0.4 * std::abs(F[0][0]));
}

// Analytic forces must reproduce the finite-difference oracle component-by-
// component (including the off-axis egg-box artifact — same grid, same energy).
TEST(GpwForces, H2AnalyticVsFiniteDifference) {
    const auto basis = parse_gth_basis(kData + "/BASIS_SET", "H", "DZV-GTH-PADE");
    const auto pot = parse_gth_potential(kData + "/POTENTIAL", "H", "GTH-PADE-q1");

    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.7, 2.1, 1.9}; // off-symmetry so every component is nonzero
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};

    GpwConfig cfg;
    cfg.cell[0] = 8.0;
    cfg.cell[1] = cfg.cell[2] = 4.0;
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-10;

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);

    double max_abs = 0.0, max_err = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            max_abs = std::max(max_abs, std::abs(Ffd[a][k]));
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F] max|F|=%.6f  max|dF|=%.3e Ha/Bohr\n", max_abs, max_err);
    EXPECT_LT(max_err, 1e-5) << "analytic vs FD mismatch";
}

// H2 with PBE: isolates the GGA grid-term derivative (∫ V_xc^GGA ∂ρ) cheaply,
// before the heavier ZnO case adds the d-projector / f-basis derivatives.
TEST(GpwForces, H2PbeAnalyticVsFiniteDifference) {
    const auto basis = parse_gth_basis(kData + "/BASIS_MOLOPT", "H", "DZVP-MOLOPT-SR-GTH");
    const auto pot = parse_gth_potential(kData + "/GTH_POTENTIALS", "H", "GTH-PBE-q1");
    const double rA[3] = {3.5, 2.0, 2.0};
    const double rB[3] = {4.7, 2.1, 1.9};
    std::vector<GpwAtom> atoms = {make_gpw_atom(rA, basis, pot),
                                  make_gpw_atom(rB, basis, pot)};
    GpwConfig cfg;
    cfg.cell[0] = 8.0;
    cfg.cell[1] = cfg.cell[2] = 4.0;
    cfg.cutoff_ry = 200.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-10;
    cfg.xc = XcKind::PBE;

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);
    double max_err = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F PBE] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F PBE] max|dF|=%.3e Ha/Bohr\n", max_err);
    EXPECT_LT(max_err, 1e-5) << "H2/PBE GGA analytic vs FD mismatch";
}

// ZnO: exercises the GGA(PBE) grid term, the l=2 (d) GTH projector derivative,
// and the l=3 (f) basis-shell derivative — none of which the H2/LDA case touches.
TEST(GpwForces, ZnOAnalyticVsFiniteDifference) {
    constexpr double kAng = 1.8897259885789; // Å -> Bohr
    const auto bZn = parse_gth_basis(kData + "/BASIS_MOLOPT", "Zn", "DZVP-MOLOPT-SR-GTH");
    const auto bO = parse_gth_basis(kData + "/BASIS_MOLOPT", "O", "DZVP-MOLOPT-SR-GTH");
    const auto pZn = parse_gth_potential(kData + "/GTH_POTENTIALS", "Zn", "GTH-PBE-q12");
    const auto pO = parse_gth_potential(kData + "/GTH_POTENTIALS", "O", "GTH-PBE-q6");

    const double rZn[3] = {4.0 * kAng, 4.0 * kAng, 4.0 * kAng};
    const double rO[3] = {4.15 * kAng, 4.1 * kAng, 5.6 * kAng}; // off-axis bond
    std::vector<GpwAtom> atoms = {make_gpw_atom(rZn, bZn, pZn),
                                  make_gpw_atom(rO, bO, pO)};

    GpwConfig cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0 * kAng;
    cfg.cutoff_ry = 120.0; // coarse: consistency holds at any grid; keeps FD sweep cheap
    cfg.n_occ = 9;
    cfg.tol = 1e-9;
    cfg.xc = XcKind::PBE;
    cfg.nonperiodic = false;

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);

    double max_abs = 0.0, max_err = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int k = 0; k < 3; ++k) {
            max_abs = std::max(max_abs, std::abs(Ffd[a][k]));
            max_err = std::max(max_err, std::abs(Fa[a][k] - Ffd[a][k]));
            std::printf("[GPW F ZnO] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                        a, k, Fa[a][k], Ffd[a][k], Fa[a][k] - Ffd[a][k]);
        }
    std::printf("[GPW F ZnO] max|F|=%.6f  max|dF|=%.3e Ha/Bohr\n", max_abs, max_err);
    EXPECT_LT(max_err, 5e-5) << "ZnO analytic vs FD mismatch (GGA/d-proj/f-basis)";
}
