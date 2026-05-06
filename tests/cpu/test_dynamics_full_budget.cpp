/// @file tests/cpu/test_dynamics_full_budget.cpp
/// @brief Phase 1.8 — DynamicsDriver runs the full step budget.
///
/// Two cases:
///   A) A ConvergenceCheck that sets `converged=true` every step is *ignored*
///      by Dynamics — the loop runs all `n_steps`.
///   B) A ConvergenceCheck that sets `stop_requested=true` on a specific step
///      *does* break the loop early.
///
/// Acceptance: ac-005 (Dynamics half) in
/// `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#include <gtest/gtest.h>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace {

struct StubIntegrator : ATV::cpu::ComponentBase {
    int step_count = 0;
    std::string name() const override { return "stub_integrator"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Integrator; }
    void launch(const ATV::cpu::LaunchContext&) override { ++step_count; }
};

struct StubForceEvaluator : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_force"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Evaluator; }
    ATV::Capability   capabilities() const override {
        return ATV::Capability::Energy | ATV::Capability::Force;
    }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

struct AlwaysConverged : ATV::cpu::ComponentBase {
    std::string name() const override { return "always_converged"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::ConvergenceCheck;
    }
    void launch(const ATV::cpu::LaunchContext& ctx) override {
        ctx.runtime.outer_status.converged = true;
    }
};

struct StopAtStep : ATV::cpu::ComponentBase {
    int stop_at;
    explicit StopAtStep(int s) : stop_at(s) {}
    std::string name() const override { return "stop_at_step"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::ConvergenceCheck;
    }
    void launch(const ATV::cpu::LaunchContext& ctx) override {
        if (ctx.step >= stop_at) {
            ctx.runtime.outer_status.stop_requested = true;
            ctx.runtime.outer_status.reason = "stop-requested";
        }
    }
};

ATV::cpu::Frame make_trivial_frame() {
    ATV::cpu::Frame f(1);
    f.atoms.x[0] = 0.0;
    f.atoms.y[0] = 0.0;
    f.atoms.z[0] = 0.0;
    f.atoms.mass[0] = 1.0;
    f.atoms.type[0] = 0;
    f.box = ATV::Mat3::cubic(10.0);
    return f;
}

} // namespace

TEST(DynamicsFullBudget, IgnoresConvergedFlag) {
    using namespace ATV;
    cpu::DynamicsDriver d;
    StubIntegrator integ;
    StubForceEvaluator force;
    AlwaysConverged conv;

    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(force);
    d.on(md::Phase::end_of_step).add(conv);
    d.build();

    auto frame = make_trivial_frame();
    d.run(frame, /*n_steps=*/5);
    EXPECT_EQ(integ.step_count, 5)
        << "Dynamics must ignore outer_status.converged and run all 5 steps";
}

TEST(DynamicsFullBudget, BreaksOnStopRequested) {
    using namespace ATV;
    cpu::DynamicsDriver d;
    StubIntegrator integ;
    StubForceEvaluator force;
    StopAtStep stopper{/*stop_at=*/3};

    d.on(md::Phase::initial_integrate).add(integ);
    d.on(md::Phase::force).add(force);
    d.on(md::Phase::end_of_step).add(stopper);
    d.build();

    auto frame = make_trivial_frame();
    d.run(frame, /*n_steps=*/10);
    // Driver runs step 0,1,2,3 — at step 3 the stopper flips
    // stop_requested=true, the driver checks after the step and breaks.
    EXPECT_EQ(integ.step_count, 4) << "Expected break after step 3 fires stopper";
}
