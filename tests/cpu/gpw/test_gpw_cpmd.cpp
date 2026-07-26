/// @file tests/cpu/gpw/test_gpw_cpmd.cpp
/// @brief Unit tests: GPWEvaluator CPMD build() guards + the multi-k
///        fixed-density H_k/S_k export the CPMD path consumes.
///
/// The Cpmd / CpmdCplx propagator classes themselves live in
/// cpu/electronic/ and are pinned by tests/cpu/electronic/test_cpmd.cpp.
/// No multi-step MD / full water SCF trajectory here (that is the e2e
/// regressions/dynamics/cpmd-water and cpmd-rattle-water pins). build()
/// validation and the tiny inline H2 payload only.

#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"

// ── GPWEvaluator cpmd guards ────────────────────────────────────────────

namespace {

ATV::cpu::Frame water_frame() {
    ATV::cpu::Frame f(3);
    const int z[3] = {8, 1, 1};
    const double x[3] = {5.0000, 5.7586, 4.2414};
    const double y[3] = {5.0000, 5.0000, 5.0000};
    const double zc[3] = {5.0000, 5.5043, 5.5043};
    for (int i = 0; i < 3; ++i) {
        f.span(ATV::kAtom_type)[i] = z[i];
        f.span(ATV::kAtom_x)[i] = x[i];
        f.span(ATV::kAtom_y)[i] = y[i];
        f.span(ATV::kAtom_z)[i] = zc[i];
    }
    f.box = ATV::cpu::Box::cubic(10.0);
    return f;
}

ATV::cpu::gpw::GPWEvaluator water_eval() {
    const ATV::cpu::gpw::Kind O{
        molrs::Element::O,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
    const ATV::cpu::gpw::Kind H{
        molrs::Element::H,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};
    return ATV::cpu::gpw::GPWEvaluator(
        {O, H}, /*cutoff_ry*/ 200.0, /*max_iter*/ 60, /*tol*/ 1e-5,
        /*mix*/ 0.3, ATV::cpu::gpw::XcKind::PBE, /*nonperiodic*/ true,
        /*compute_forces*/ true);
}

} // namespace

TEST(GpwEvaluatorCpmd, BuildAcceptsCpmdWithUks) {
    // Γ-UKS CPMD is product (dual α/β propagators); O2-like occupations.
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorCpmd, BuildRejectsUksCpmdWithZeroAlpha) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 0;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithOt) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithXlBomd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.xl_bomd = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsNonPositiveMuOrDt) {
    auto frame = water_frame();
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_mu = 0.0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_dt = -0.1;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
}

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithoutForces) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.compute_forces = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// ── Nosé–Hoover electronic thermostat guards ────────────────────────────

TEST(GpwEvaluatorCpmd, BuildRejectsThermostatWithoutCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = false;
    eval.cpmd_thermostat = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsBadThermostatParams) {
    auto frame = water_frame();
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_thermostat = true;
        eval.cpmd_ekin_target_ha = 0.0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_thermostat = true;
        eval.cpmd_tau_e = -5.0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
}

TEST(GpwEvaluatorCpmd, BuildAcceptsCpmdWithThermostat) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.cpmd_thermostat = true;
    eval.cpmd_ekin_target_ha = 1e-5;
    eval.cpmd_tau_e = 50.0;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

// ── Multi-k Bloch CPMD evaluator guards ─────────────────────────────────

TEST(GpwEvaluatorCpmd, BuildAcceptsMultiKRksCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.k_mesh[0] = 1;
    eval.k_mesh[1] = 1;
    eval.k_mesh[2] = 2;
    eval.k_time_reversal = false;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorCpmd, BuildRejectsMultiKUksCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.k_mesh[2] = 2;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsMultiKCpmdThermostat) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.cpmd_thermostat = true;
    eval.k_mesh[2] = 2;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// Multi-k fixed-density payload (tiny H2 GpwDft).
namespace {

ATV::cpu::gpw::GridBasisFunction s_basis_cp(double alpha, const double c[3]) {
    ATV::cpu::gpw::GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = c[0];
    p.center[1] = c[1];
    p.center[2] = c[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = std::pow(2.0 * alpha / M_PI, 0.75);
    return ATV::cpu::gpw::GridBasisFunction{p};
}

ATV::cpu::gpw::GpwAtom hydrogen_cp(const double pos[3]) {
    ATV::cpu::gpw::GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis_cp(0.5, pos)};
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

std::vector<ATV::cpu::gpw::GpwAtom> small_h2_cp(double box = 8.0) {
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    return {hydrogen_cp(rA), hydrogen_cp(rB)};
}

ATV::cpu::gpw::GpwDft small_cfg_cp(double box = 8.0) {
    ATV::cpu::gpw::GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = box;
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 1;
    cfg.max_iter = 80;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;
    return cfg;
}

} // namespace

TEST(GpwSolverFixedDensityK, MultiKFixedDensityExportsHkSk) {
    const auto atoms = small_h2_cp(/*box=*/8.0);
    ATV::cpu::gpw::GpwDft cold = small_cfg_cp(/*box=*/8.0);
    cold.k_mesh[0] = 1;
    cold.k_mesh[1] = 1;
    cold.k_mesh[2] = 2;
    cold.k_time_reversal = false;
    cold.nonperiodic = false;
    cold.build(atoms);
    const auto ref = cold.run(atoms);
    ASSERT_TRUE(ref.converged);
    ASSERT_GE(ref.n_kpts, 2);
    ASSERT_EQ(static_cast<int>(ref.P_k.size()), ref.n_kpts);

    ATV::cpu::gpw::GpwDft fix = small_cfg_cp(/*box=*/8.0);
    fix.k_mesh[0] = 1;
    fix.k_mesh[1] = 1;
    fix.k_mesh[2] = 2;
    fix.k_time_reversal = false;
    fix.nonperiodic = false;
    fix.max_iter = 0;
    fix.build(atoms);
    const auto& r = fix.run(atoms, nullptr, nullptr, nullptr, &ref.P_k);
    EXPECT_TRUE(r.converged);
    EXPECT_EQ(r.iters, 0);
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_EQ(static_cast<int>(r.H_ks_k.size()), r.n_kpts);
    EXPECT_EQ(static_cast<int>(r.S_k.size()), r.n_kpts);
    EXPECT_EQ(static_cast<int>(r.P_k.size()), r.n_kpts);
    // Energy near the cold SCF (same density).
    EXPECT_NEAR(r.e_total, ref.e_total, 1e-6);
}

TEST(GpwSolverFixedDensityK, MultiKFixedDensityWithoutGuessThrows) {
    const auto atoms = small_h2_cp();
    ATV::cpu::gpw::GpwDft fix = small_cfg_cp();
    fix.k_mesh[2] = 2;
    fix.k_time_reversal = false;
    fix.nonperiodic = false;
    fix.max_iter = 0;
    fix.build(atoms);
    EXPECT_THROW(fix.run(atoms), std::invalid_argument);
}
