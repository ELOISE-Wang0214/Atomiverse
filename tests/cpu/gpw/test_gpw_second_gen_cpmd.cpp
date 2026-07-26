/// @file tests/cpu/gpw/test_gpw_second_gen_cpmd.cpp
/// @brief Unit tests: GPWEvaluator second-generation CP (ASPC + partial SCF)
///        build() guards only. No multi-step MD here (e2e is
///        regressions/dynamics/second-gen-cpmd-water and second-gen-cpmd-uks).

#include <stdexcept>

#include <gtest/gtest.h>

#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gth_potential.h"

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

TEST(GpwEvaluatorSecondGenCp, BuildAcceptsWithUks) {
    // Γ-UKS dual ASPC (product path, wave-7 second-gen-cpmd-uks).
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorSecondGenCp, BuildAcceptsUksWithOt) {
    // Γ-UKS OT + dual ASPC corrector (product path, second-gen-cpmd-uks-ot).
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorSecondGenCp, BuildAcceptsWithOt) {
    // RKS OT partial-SCF corrector is product (second-gen-cpmd-ot).
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorSecondGenCp, BuildRejectsWithCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.cpmd = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorSecondGenCp, BuildRejectsWithXlBomd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.xl_bomd = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorSecondGenCp, BuildRejectsWithoutForces) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.compute_forces = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorSecondGenCp, BuildRejectsBadOrderOrBudget) {
    auto frame = water_frame();
    {
        auto eval = water_eval();
        eval.second_gen_cp = true;
        eval.second_gen_cp_k = -1;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
    {
        auto eval = water_eval();
        eval.second_gen_cp = true;
        eval.second_gen_cp_max_iter = 0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
}

TEST(GpwEvaluatorSecondGenCp, BuildAcceptsDefaultMembers) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.second_gen_cp = true;
    eval.second_gen_cp_k = 3;
    eval.second_gen_cp_max_iter = 5;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorSecondGenCp, BuildAcceptsMultiK) {
    // Multi-k RKS ASPC on complex P_k (product path, second-gen-cpmd-kpoints).
    auto frame = water_frame();
    // Periodic cell required for a k-mesh.
    frame.box = ATV::cpu::Box::cubic(10.0);
    auto eval = water_eval();
    eval.nonperiodic = false;
    eval.second_gen_cp = true;
    eval.k_mesh[0] = 1;
    eval.k_mesh[1] = 1;
    eval.k_mesh[2] = 2;
    eval.k_time_reversal = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorSecondGenCp, BuildRejectsMultiKWithOt) {
    auto frame = water_frame();
    frame.box = ATV::cpu::Box::cubic(10.0);
    auto eval = water_eval();
    eval.nonperiodic = false;
    eval.second_gen_cp = true;
    eval.use_ot = true;
    eval.k_mesh[0] = 1;
    eval.k_mesh[1] = 1;
    eval.k_mesh[2] = 2;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorSecondGenCp, BuildAcceptsMultiKUks) {
    auto frame = water_frame();
    frame.box = ATV::cpu::Box::cubic(10.0);
    auto eval = water_eval();
    eval.nonperiodic = false;
    eval.second_gen_cp = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.k_mesh[0] = 1;
    eval.k_mesh[1] = 1;
    eval.k_mesh[2] = 2;
    eval.k_time_reversal = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}
