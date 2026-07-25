/// @file tests/cpu/gpw/test_gpw_cpmd.cpp
/// @brief Unit tests: Cpmd electronic propagator + GPWEvaluator CPMD guards.
///
/// No multi-step MD / full water SCF trajectory here (that is the e2e
/// regressions/dynamics/cpmd-water pin). Matrix-level and build() validation only.

#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>
#include <xtensor/generators/xbuilder.hpp>

#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/cpmd.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/linalg/matrix.h"

using ATV::cpu::gpw::Cpmd;
using ATV::cpu::linalg::RealMatrix;

// ── Cpmd propagator ─────────────────────────────────────────────────────

TEST(Cpmd, ConstructsWithPositiveMass) {
    const Cpmd cp(400.0);
    EXPECT_DOUBLE_EQ(cp.mu(), 400.0);
    EXPECT_FALSE(cp.ready());
}

TEST(Cpmd, NonPositiveMassThrows) {
    EXPECT_THROW(Cpmd(0.0), std::invalid_argument);
    EXPECT_THROW(Cpmd(-10.0), std::invalid_argument);
}

TEST(Cpmd, SeedThenReadyAtRest) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({4, 2});
    C(0, 0) = 1.0;
    C(1, 1) = 1.0;
    cp.seed(C);
    ASSERT_TRUE(cp.ready());
    EXPECT_EQ(cp.C().shape(0), 4u);
    EXPECT_EQ(cp.C().shape(1), 2u);
    EXPECT_NEAR(cp.kinetic_ha(), 0.0, 1e-15);
}

TEST(Cpmd, AdvanceBeforeSeedThrows) {
    Cpmd cp(400.0);
    RealMatrix H = xt::zeros<double>({2, 2});
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = 1.0;
    S(1, 1) = 1.0;
    EXPECT_THROW(cp.advance(H, S, 0.1), std::logic_error);
}

TEST(Cpmd, AdvanceWithZeroHKeepsRest) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    RealMatrix H = xt::zeros<double>({2, 2});
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = 1.0;
    S(1, 1) = 1.0;
    cp.advance(H, S, 0.1);
    EXPECT_NEAR(cp.C()(0, 0), 1.0, 1e-12);
    EXPECT_NEAR(cp.C()(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(cp.kinetic_ha(), 0.0, 1e-12);
}

TEST(Cpmd, EmptySeedThrows) {
    Cpmd cp(400.0);
    RealMatrix empty;
    EXPECT_THROW(cp.seed(empty), std::invalid_argument);
}

TEST(Cpmd, ResetDropsReady) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    cp.reset();
    EXPECT_FALSE(cp.ready());
    EXPECT_THROW(cp.C(), std::logic_error);
}

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
    f.simbox = ATV::cpu::Box::cubic(10.0);
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

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithUks) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 5;
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
