/// @file tests/cpu/pw/test_pw_forces.cpp
/// @brief Unit gates for plane-wave FD ionic forces (cheap toy settings).

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/pw/pw_evaluator.h"
#include "atomiverse/cpu/pw/pw_solver.h"
#include "atomiverse/facade/stage.h"

using namespace ATV;
using namespace ATV::cpu::pw;

namespace {

/// Single H atom in a cubic cell — force must vanish by lattice symmetry.
cpu::Frame make_h_cubic(double L) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.simbox = cpu::Box::cubic(L);
    return f;
}

PWConfig cheap_cfg() {
    PWConfig cfg;
    cfg.ecut_ry = 5.0;
    cfg.k_mesh = 1;
    cfg.tol = 1e-4;
    cfg.max_iter = 6;
    cfg.mix_beta = 0.5;
    return cfg;
}

} // namespace

TEST(PWForces, RejectsBadStep) {
    auto f = make_h_cubic(6.0);
    EXPECT_THROW(pw_forces_fd(f, cheap_cfg(), 0.0), std::invalid_argument);
    EXPECT_THROW(pw_forces_fd(f, cheap_cfg(), -1e-3), std::invalid_argument);
}

TEST(PWForces, SingleAtomForceVanishes) {
    // 6 SCF solves; cheap settings keep this under ~30 s on a laptop.
    auto f = make_h_cubic(6.0);
    const auto F = pw_forces_fd(f, cheap_cfg(), 5e-3);
    ASSERT_EQ(F.size(), 1u);
    const double mag =
        std::sqrt(F[0][0] * F[0][0] + F[0][1] * F[0][1] + F[0][2] * F[0][2]);
    // Toy PP + loose SCF: force is numerical noise around zero.
    EXPECT_LT(mag, 0.5) << " |F|=" << mag;
}

TEST(PWForces, EvaluatorWritesForcesWhenRequested) {
    auto frame = make_h_cubic(6.0);
    for (int i = 0; i < frame.nlocal(); ++i) {
        frame.span(kAtom_fx)[i] = 99.0;
        frame.span(kAtom_fy)[i] = 99.0;
        frame.span(kAtom_fz)[i] = 99.0;
    }

    PWEvaluator pwe{cheap_cfg()};
    pwe.with_gradient(true).set_fd_step(5e-3);

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(pwe);
    eval.enable_timing(false);
    eval.build();
    eval.run(frame);

    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
    // Forces overwritten (not left at 99).
    EXPECT_NE(frame.span(kAtom_fx)[0], AccReal(99.0));
    const double mag = std::sqrt(static_cast<double>(
        frame.span(kAtom_fx)[0] * frame.span(kAtom_fx)[0]
        + frame.span(kAtom_fy)[0] * frame.span(kAtom_fy)[0]
        + frame.span(kAtom_fz)[0] * frame.span(kAtom_fz)[0]));
    EXPECT_LT(mag, 0.5);
}

TEST(PWForces, EvaluatorEnergyOnlyDoesNotSpendFdForces) {
    auto frame = make_h_cubic(6.0);
    PWEvaluator pwe{cheap_cfg()};
    EXPECT_FALSE(pwe.compute_forces()); // default off — no 6N SCF force path

    cpu::Evaluate eval;
    eval.on(evaluate::Phase::evaluate).add_ref(pwe);
    eval.enable_timing(false);
    eval.build();
    eval.run(frame);
    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
}
