/// @file tests/cpu/gpw/test_cpu_gpw_stress.cpp
/// @brief Unit gates for GPW cell stress (analytic fixed-P + multi-SCF FD oracle).
///
/// E2E CP2K / analytic-vs-FD water twin lives in
/// regressions/gpw/gpw-analytic-stress/ and regressions/gpw/gpw-stress-water/.

#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/evaluate/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"

using namespace ATV;
using namespace ATV::cpu::gpw;

namespace {

GpwDft cheap_h2_cfg() {
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0; // Bohr
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.max_iter = 80;
    cfg.nonperiodic = true;
    cfg.local_pp_on_grid = true;
    cfg.xc = XcKind::PADE;
    return cfg;
}

std::vector<GpwAtom> cheap_h2_atoms() {
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.7, 6.1, 5.9};
    return {GpwAtom(rA, basis, pot, 1), GpwAtom(rB, basis, pot, 1)};
}

} // namespace

TEST(GpwStress, AnalyticOffDiagonalZeroOrthorhombic) {
    const auto atoms = cheap_h2_atoms();
    const auto s = gpw_stress(atoms, cheap_h2_cfg(), 3e-3);
    for (int k = 0; k < 6; ++k)
        EXPECT_TRUE(std::isfinite(s[static_cast<std::size_t>(k)])) << k;
    EXPECT_DOUBLE_EQ(s[3], 0.0);
    EXPECT_DOUBLE_EQ(s[4], 0.0);
    EXPECT_DOUBLE_EQ(s[5], 0.0);
    const double max_diag =
        std::max({std::abs(s[0]), std::abs(s[1]), std::abs(s[2])});
    std::printf(
        "[GPW stress analytic] Voigt Ha/Bohr3: %.6e %.6e %.6e | 0 0 0\n", s[0],
        s[1], s[2]);
    EXPECT_GT(max_diag, 1e-8);
}

TEST(GpwStress, AnalyticMatchesFdDiagonal) {
    // Hellmann–Feynman fixed-P product vs multi-SCF strain-FD oracle on the
    // same cheap H₂ fixture. Residual is SCF residual / O(h²); gate is unit-
    // tier loose (e2e water twin uses 1e-5).
    const auto atoms = cheap_h2_atoms();
    const auto cfg = cheap_h2_cfg();
    const double h = 3e-3;
    const auto san = gpw_stress(atoms, cfg, h);
    const auto sfd = gpw_stress_fd(atoms, cfg, h);
    double max_err = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double err = std::abs(san[static_cast<std::size_t>(a)]
                                    - sfd[static_cast<std::size_t>(a)]);
        max_err = std::max(max_err, err);
        std::printf("  axis %d  an=%.6e  fd=%.6e  |d|=%.3e\n", a,
                    san[static_cast<std::size_t>(a)],
                    sfd[static_cast<std::size_t>(a)], err);
    }
    EXPECT_LT(max_err, 5e-5) << "analytic vs multi-SCF FD diagonal";
}

namespace {

/// Toy single-s-primitive periodic H₂ — keeps the axes-mask gate at unit-tier
/// cost (small grid, seconds) unlike the DZV MT fixture above.
std::vector<GpwAtom> toy_h2_atoms() {
    auto s_prim = [](double alpha, const double pos[3]) {
        GridBasisFunction bf;
        GridPrimitive p;
        p.alpha = alpha;
        p.center[0] = pos[0];
        p.center[1] = pos[1];
        p.center[2] = pos[2];
        p.l[0] = p.l[1] = p.l[2] = 0;
        p.coeff = 1.0;
        bf.push_back(p);
        return bf;
    };
    auto hydrogen = [&](const double pos[3]) {
        GpwAtom a;
        a.z_val = 1.0;
        a.z = 1;
        a.pos[0] = pos[0];
        a.pos[1] = pos[1];
        a.pos[2] = pos[2];
        a.basis = {s_prim(0.5, pos)};
        a.rloc = 0.2;
        a.cloc = {-4.17890044, 0.72446331};
        return a;
    };
    const double rA[3] = {3.3, 4.0, 4.0};
    const double rB[3] = {4.7, 4.1, 3.9};
    return {hydrogen(rA), hydrogen(rB)};
}

GpwDft toy_h2_cfg() {
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 8.0; // Bohr, periodic
    cfg.cutoff_ry = 80.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.max_iter = 80;
    cfg.xc = XcKind::PADE;
    return cfg;
}

} // namespace

TEST(GpwStress, AxesMaskSkipsFrozenAxesAnalyticAndFd) {
    const auto atoms = toy_h2_atoms();
    const auto cfg = toy_h2_cfg();
    const double h = 3e-3;
    const auto full_an = gpw_stress(atoms, cfg, h);
    const auto x_an = gpw_stress(atoms, cfg, h, /*axes*/ {true, false, false});
    EXPECT_NEAR(x_an[0], full_an[0], 1e-12);
    EXPECT_DOUBLE_EQ(x_an[1], 0.0);
    EXPECT_DOUBLE_EQ(x_an[2], 0.0);

    const auto full_fd = gpw_stress_fd(atoms, cfg, h);
    const auto x_fd =
        gpw_stress_fd(atoms, cfg, h, /*axes*/ {true, false, false});
    EXPECT_NEAR(x_fd[0], full_fd[0], 1e-12);
    EXPECT_DOUBLE_EQ(x_fd[1], 0.0);
    EXPECT_DOUBLE_EQ(x_fd[2], 0.0);

    EXPECT_THROW(gpw_stress(atoms, cfg, h, {false, false, false}),
                 std::invalid_argument);
    EXPECT_THROW(gpw_stress_fd(atoms, cfg, h, {false, false, false}),
                 std::invalid_argument);
}

TEST(GpwStress, EvaluatorWritesKStressAnalyticDefault) {
    cpu::Frame frame(2);
    frame.span(kAtom_type)[0] = 1;
    frame.span(kAtom_type)[1] = 1;
    constexpr double b2a = 1.0 / 1.8897261245650618;
    frame.span(kAtom_x)[0] = 5.5 * b2a;
    frame.span(kAtom_y)[0] = 6.0 * b2a;
    frame.span(kAtom_z)[0] = 6.0 * b2a;
    frame.span(kAtom_x)[1] = 6.7 * b2a;
    frame.span(kAtom_y)[1] = 6.1 * b2a;
    frame.span(kAtom_z)[1] = 5.9 * b2a;
    frame.simbox = cpu::Box::cubic(12.0 * b2a);

    const Kind H{molrs::Element::H,
                 GaussianBasis(molrs::Element::H, "DZV-GTH-PADE"),
                 GthPotential(molrs::Element::H, "GTH-PADE-q1")};

    GPWEvaluator eval(/*kinds*/ {H}, /*cutoff_ry*/ 120.0, /*max_iter*/ 80,
                      /*tol*/ 1e-8, /*mix*/ 0.3, /*xc*/ XcKind::PADE,
                      /*nonperiodic*/ true, /*compute_forces*/ false);
    eval.compute_stress = true;
    eval.stress_use_fd = false; // product default
    eval.stress_strain_h = 3e-3;

    cpu::Evaluate driver;
    driver.on(evaluate::Phase::evaluate).add_ref(eval);
    driver.enable_timing(false);
    driver.build();
    driver.run(frame);

    EXPECT_TRUE(eval.last_result().converged);
    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
    const auto& st = frame.meta.get(kStress);
    for (int k = 0; k < 6; ++k)
        EXPECT_TRUE(std::isfinite(st[static_cast<std::size_t>(k)])) << k;
    EXPECT_NEAR(st[3], 0.0, 1e-14);
    EXPECT_NEAR(st[4], 0.0, 1e-14);
    EXPECT_NEAR(st[5], 0.0, 1e-14);
    const auto& sh = eval.last_stress_ha_bohr3();
    EXPECT_GT(std::max({std::abs(sh[0]), std::abs(sh[1]), std::abs(sh[2])}),
              1e-8);
}

// ── Multi-k analytic stress ────────────────────────────────────────────────

TEST(GpwStress, KpointAnalyticMatchesFdDiagonal) {
    // Periodic H2, mesh 1x1x2 (TR → 1 k at κ_z = +1/4). Hellmann–Feynman
    // fixed-P_k product vs multi-SCF strain-FD oracle.
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.5, 6.0, 6.0};
    const std::vector<GpwAtom> atoms = {GpwAtom(rA, basis, pot, 1),
                                        GpwAtom(rB, basis, pot, 1)};

    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0;
    cfg.cutoff_ry = 100.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-8;
    cfg.max_iter = 60;
    cfg.nonperiodic = false;
    cfg.local_pp_on_grid = true;
    cfg.xc = XcKind::PADE;
    cfg.eps_default = 1e-8;
    cfg.k_mesh[0] = 1;
    cfg.k_mesh[1] = 1;
    cfg.k_mesh[2] = 2;

    const double h = 3e-3;
    const auto san = gpw_stress(atoms, cfg, h);
    const auto sfd = gpw_stress_fd(atoms, cfg, h);
    double max_err = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double err = std::abs(san[static_cast<std::size_t>(a)]
                                    - sfd[static_cast<std::size_t>(a)]);
        max_err = std::max(max_err, err);
        std::printf("  k-mesh axis %d  an=%.6e  fd=%.6e  |d|=%.3e\n", a,
                    san[static_cast<std::size_t>(a)],
                    sfd[static_cast<std::size_t>(a)], err);
        EXPECT_TRUE(std::isfinite(san[static_cast<std::size_t>(a)]));
        EXPECT_TRUE(std::isfinite(sfd[static_cast<std::size_t>(a)]));
    }
    EXPECT_DOUBLE_EQ(san[3], 0.0);
    EXPECT_DOUBLE_EQ(san[4], 0.0);
    EXPECT_DOUBLE_EQ(san[5], 0.0);
    EXPECT_LT(max_err, 5e-5) << "k-mesh analytic vs multi-SCF FD diagonal";
}

TEST(GpwStress, KpointNoThrowAndFinite) {
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double rA[3] = {5.5, 6.0, 6.0};
    const double rB[3] = {6.7, 6.1, 5.9};
    const std::vector<GpwAtom> atoms = {GpwAtom(rA, basis, pot, 1),
                                        GpwAtom(rB, basis, pot, 1)};
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = 12.0;
    cfg.cutoff_ry = 80.0;
    cfg.n_occ = 1;
    cfg.tol = 1e-7;
    cfg.max_iter = 50;
    cfg.xc = XcKind::PADE;
    cfg.local_pp_on_grid = true;
    cfg.k_mesh[0] = cfg.k_mesh[1] = cfg.k_mesh[2] = 2; // MP 2x2x2 TR
    EXPECT_NO_THROW({
        const auto s = gpw_stress(atoms, cfg, 4e-3);
        for (int k = 0; k < 6; ++k)
            EXPECT_TRUE(std::isfinite(s[static_cast<std::size_t>(k)])) << k;
        EXPECT_DOUBLE_EQ(s[3], 0.0);
        EXPECT_DOUBLE_EQ(s[4], 0.0);
        EXPECT_DOUBLE_EQ(s[5], 0.0);
    });
}

// ── Γ-UKS analytic stress ──────────────────────────────────────────────────

TEST(GpwStress, UksGammaAnalyticMatchesFdDiagonal) {
    // High-spin H₂ (Nα=2, Nβ=0) Γ-UKS: fixed-P spin Lagrangian vs multi-SCF FD.
    const auto atoms = cheap_h2_atoms();
    GpwDft cfg = cheap_h2_cfg();
    cfg.n_occ = 0;
    cfg.n_alpha = 2;
    cfg.n_beta = 0;
    cfg.xc = XcKind::PADE;
    const double h = 3e-3;
    const auto san = gpw_stress(atoms, cfg, h);
    const auto sfd = gpw_stress_fd(atoms, cfg, h);
    double max_err = 0.0;
    for (int a = 0; a < 3; ++a) {
        const double err = std::abs(san[static_cast<std::size_t>(a)]
                                    - sfd[static_cast<std::size_t>(a)]);
        max_err = std::max(max_err, err);
        std::printf("  UKS axis %d  an=%.6e  fd=%.6e  |d|=%.3e\n", a,
                    san[static_cast<std::size_t>(a)],
                    sfd[static_cast<std::size_t>(a)], err);
        EXPECT_TRUE(std::isfinite(san[static_cast<std::size_t>(a)]));
        EXPECT_TRUE(std::isfinite(sfd[static_cast<std::size_t>(a)]));
    }
    EXPECT_DOUBLE_EQ(san[3], 0.0);
    EXPECT_DOUBLE_EQ(san[4], 0.0);
    EXPECT_DOUBLE_EQ(san[5], 0.0);
    const double max_diag =
        std::max({std::abs(san[0]), std::abs(san[1]), std::abs(san[2])});
    EXPECT_GT(max_diag, 1e-8);
    EXPECT_LT(max_err, 5e-5) << "UKS analytic vs multi-SCF FD diagonal";
}

TEST(GpwStress, UksKpointThrows) {
    const auto atoms = cheap_h2_atoms();
    GpwDft cfg = cheap_h2_cfg();
    cfg.n_occ = 0;
    cfg.n_alpha = 2;
    cfg.n_beta = 0;
    cfg.nonperiodic = false;
    cfg.k_mesh[0] = cfg.k_mesh[1] = cfg.k_mesh[2] = 2;
    EXPECT_THROW(gpw_stress(atoms, cfg, 3e-3), std::invalid_argument);
}
