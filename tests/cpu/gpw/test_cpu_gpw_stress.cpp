/// @file tests/cpu/gpw/test_cpu_gpw_stress.cpp
/// @brief Unit gates for GPW cell stress (orthorhombic strain-FD v1).
///
/// E2E CP2K STRESS_TENSOR parity lives in regressions/gpw/gpw-stress-water/.

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

TEST(GpwStress, OffDiagonalZeroOrthorhombicV1) {
    const auto atoms = cheap_h2_atoms();
    const auto s = gpw_stress(atoms, cheap_h2_cfg(), 3e-3);
    for (int k = 0; k < 6; ++k)
        EXPECT_TRUE(std::isfinite(s[static_cast<std::size_t>(k)])) << k;
    // v1: shear components identically zero (no triclinic cell).
    EXPECT_DOUBLE_EQ(s[3], 0.0);
    EXPECT_DOUBLE_EQ(s[4], 0.0);
    EXPECT_DOUBLE_EQ(s[5], 0.0);
    // Diagonal must be finite and nonzero for a confined molecule.
    const double max_diag =
        std::max({std::abs(s[0]), std::abs(s[1]), std::abs(s[2])});
    std::printf("[GPW stress] Voigt Ha/Bohr3: %.6e %.6e %.6e | 0 0 0\n", s[0],
                s[1], s[2]);
    EXPECT_GT(max_diag, 1e-8);
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

TEST(GpwStress, AxesMaskSkipsFrozenAxes) {
    // Masked axes report exactly 0 (no FD solves); the unmasked axis is the
    // same central FD as the full call — identical solver path, so equal to
    // tight tolerance.
    const auto atoms = toy_h2_atoms();
    const auto full = gpw_stress(atoms, toy_h2_cfg(), 3e-3);
    const auto x_only = gpw_stress(atoms, toy_h2_cfg(), 3e-3,
                                   /*axes*/ {true, false, false});
    EXPECT_NEAR(x_only[0], full[0], 1e-12);
    EXPECT_GT(std::abs(full[1]) + std::abs(full[2]), 0.0); // full solves y,z
    EXPECT_DOUBLE_EQ(x_only[1], 0.0);
    EXPECT_DOUBLE_EQ(x_only[2], 0.0);
    EXPECT_DOUBLE_EQ(x_only[3], 0.0);
    EXPECT_DOUBLE_EQ(x_only[4], 0.0);
    EXPECT_DOUBLE_EQ(x_only[5], 0.0);
    // All-masked is a caller error.
    EXPECT_THROW(gpw_stress(atoms, toy_h2_cfg(), 3e-3,
                            {false, false, false}),
                 std::invalid_argument);
}

TEST(GpwStress, EvaluatorWritesKStress) {
    cpu::Frame frame(2);
    frame.span(kAtom_type)[0] = 1;
    frame.span(kAtom_type)[1] = 1;
    // Bohr coords from cheap_h2 → Å
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
