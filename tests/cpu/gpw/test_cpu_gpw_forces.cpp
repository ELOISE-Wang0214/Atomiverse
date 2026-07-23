/// @file tests/cpu/gpw/test_cpu_gpw_forces.cpp
/// @brief Unit: analytic GPW forces match the FD energy gradient (tiny H2).
///
/// E2E CP2K FORCE_EVAL parity lives in regressions/gpw/gpw-forces-water/.
/// This file only checks that gpw_forces_analytic is the gradient of the same
/// GpwDft energy (including local_pp_on_grid consistency).

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/gpw_forces.h"

using namespace ATV::cpu::gpw;

// Analytic forces must reproduce the finite-difference oracle component-by-
// component on a cheap closed-shell system (H2 / PADE, isolated).
TEST(GpwForces, H2AnalyticVsFiniteDifference) {
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");

    // Off-symmetry so every Cartesian component is nonzero.
    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.7, 6.1, 5.9};
    std::vector<GpwAtom> atoms = {GpwAtom(rA, basis, pot, 1),
                                  GpwAtom(rB, basis, pot, 1)};

    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0;
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-10;
    cfg.max_iter = 100;
    cfg.nonperiodic = true;
    cfg.local_pp_on_grid = true; // product default

    const auto Fa = gpw_forces_analytic(atoms, cfg);
    const auto Ffd = gpw_forces_fd(atoms, cfg, 1e-3);

    ASSERT_EQ(Fa.size(), 2u);
    ASSERT_EQ(Ffd.size(), 2u);

    double max_abs = 0.0;
    double max_err = 0.0;
    for (int a = 0; a < 2; ++a) {
        for (int k = 0; k < 3; ++k) {
            max_abs = std::max(max_abs, std::abs(Ffd[static_cast<std::size_t>(a)][k]));
            max_err = std::max(
                max_err, std::abs(Fa[static_cast<std::size_t>(a)][k]
                                  - Ffd[static_cast<std::size_t>(a)][k]));
            std::printf(
                "[GPW F] atom %d axis %d  analytic % .8f  fd % .8f  d % .2e\n",
                a, k, Fa[static_cast<std::size_t>(a)][k],
                Ffd[static_cast<std::size_t>(a)][k],
                Fa[static_cast<std::size_t>(a)][k]
                    - Ffd[static_cast<std::size_t>(a)][k]);
        }
    }
    std::printf("[GPW F] max|F|=%.6f  max|dF|=%.3e Ha/Bohr\n", max_abs, max_err);
    EXPECT_LT(max_err, 5e-5) << "analytic vs FD mismatch (Ha/Bohr)";
}

// Compressed H2: bond force is restoring and equal-and-opposite (Newton 3).
TEST(GpwForces, H2RestoringBond) {
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");

    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.5, 6.0, 6.0}; // 1.0 Bohr (eq ~1.4)
    std::vector<GpwAtom> atoms = {GpwAtom(rA, basis, pot, 1),
                                  GpwAtom(rB, basis, pot, 1)};

    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0;
    cfg.cutoff_ry = 140.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.max_iter = 100;
    cfg.nonperiodic = true;

    const auto F = gpw_forces_analytic(atoms, cfg);
    ASSERT_EQ(F.size(), 2u);
    EXPECT_NEAR(F[0][0], -F[1][0], 1e-4);
    EXPECT_LT(F[0][0], 0.0);
    EXPECT_GT(F[1][0], 0.0);
    EXPECT_GT(std::abs(F[0][0]), 5e-3);
}
