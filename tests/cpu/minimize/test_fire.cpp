/// @file tests/cpu/minimize/test_fire.cpp
/// @brief Unit tests for cpu::minimizer::FIRE (single-function / method cases).

#include "atomiverse/cpu/driver/criteria.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/minimize/fire.h"
#include "atomiverse/cpu/minimize/minimize.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"

#include <cmath>
#include <gtest/gtest.h>

using namespace ATV;

namespace {

cpu::Frame lj_dimer(double r) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = 39.948;
        f.span(kAtom_type)[i] = 0;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_vx)[i] = 0.0;
        f.span(kAtom_vy)[i] = 0.0;
        f.span(kAtom_vz)[i] = 0.0;
    }
    f.simbox = cpu::Box::cubic(40.0);
    return f;
}

} // namespace

TEST(TestFIRE, DefaultsMatchAseFire10) {
    cpu::minimizer::FIRE fire;
    EXPECT_DOUBLE_EQ(fire.dt0, 0.1);
    EXPECT_DOUBLE_EQ(fire.dt_max, 1.0);
    EXPECT_EQ(fire.n_min, 5);
    EXPECT_DOUBLE_EQ(fire.f_inc, 1.1);
    EXPECT_DOUBLE_EQ(fire.f_dec, 0.5);
    EXPECT_DOUBLE_EQ(fire.alpha0, 0.1);
    EXPECT_DOUBLE_EQ(fire.f_alpha, 0.99);
    EXPECT_DOUBLE_EQ(fire.max_step, 0.2);
}

TEST(TestFIRE, DefaultBindingsSearchAndConvergenceOnly) {
    cpu::minimizer::FIRE fire;
    const auto b = fire.default_bindings();
    ASSERT_EQ(b.size(), 2u);
    EXPECT_EQ(b[0].phase, ATV::slot_of(minimizer::Phase::search_direction));
    EXPECT_EQ(b[1].phase, ATV::slot_of(minimizer::Phase::convergence_check));
    for (const auto& e : b)
        EXPECT_NE(e.phase, ATV::slot_of(minimizer::Phase::line_search));
}

TEST(TestFIRE, NameAndOptimizerKind) {
    cpu::minimizer::FIRE fire;
    EXPECT_EQ(fire.name(), "fire");
    cpu::Optimizer* as_opt = &fire;
    EXPECT_NE(dynamic_cast<cpu::Optimizer*>(as_opt), nullptr);
}

TEST(TestFIRE, ShortDescentImprovesForceAndLeavesVelocitiesZero) {
    const double eps = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 2.5 * sigma;
    const double r_star = sigma * std::pow(2.0, 1.0 / 6.0);
    auto frame = lj_dimer(1.5 * r_star);

    cpu::minimizer::Minimize opt;
    cpu::minimizer::FIRE fire;
    fire.dt0 = 1.0;
    fire.dt_max = 20.0;
    opt.set_optimizer_ref(fire);
    opt.add_criterion(cpu::ForceConverged{1e-3, 5e-4, 0.2, 0.1, 1e-6});
    opt.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    opt.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.3});
    opt.build();
    opt.run(frame, 500);

    ASSERT_FALSE(fire.trace().empty());
    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
    EXPECT_TRUE(opt.converged());
    EXPECT_LT(fire.max_force(), 1e-3);

    const double dx = frame.span(kAtom_x)[1] - frame.span(kAtom_x)[0];
    const double r = std::abs(dx);
    EXPECT_NEAR(r, r_star, 0.05);

    for (int i = 0; i < 2; ++i) {
        EXPECT_DOUBLE_EQ(frame.span(kAtom_vx)[i], 0.0);
        EXPECT_DOUBLE_EQ(frame.span(kAtom_vy)[i], 0.0);
        EXPECT_DOUBLE_EQ(frame.span(kAtom_vz)[i], 0.0);
    }
}

TEST(TestFIRE, MaxStepClamp) {
    const double eps = 0.0104;
    const double sigma = 3.4;
    const double cutoff = 2.5 * sigma;
    const double r_star = sigma * std::pow(2.0, 1.0 / 6.0);
    auto frame = lj_dimer(1.5 * r_star);

    const double x0 = frame.span(kAtom_x)[1];

    cpu::minimizer::Minimize opt;
    cpu::minimizer::FIRE fire;
    fire.dt0 = 5.0;
    fire.dt_max = 10.0;
    fire.max_step = 0.05;
    opt.set_optimizer_ref(fire);
    // One outer iteration is enough to exercise the clamp.
    opt.add_criterion(cpu::ForceConverged{1e-30, 1e-30, 1e-30, 1e-30, 1e-30});
    opt.add_potential(cpu::PairLJCut{eps, sigma, cutoff});
    opt.set_neighbor(cpu::md::NeighborList<>{cutoff, 0.3});
    opt.build();
    opt.run(frame, 2);

    const double dx = std::abs(frame.span(kAtom_x)[1] - x0);
    // Per-atom clamp: after a few steps total move can exceed max_step, but
    // each individual FIRE step is clamped; with budget 2 expect modest move.
    EXPECT_LE(dx, 0.05 * 2 + 1e-9);
}
