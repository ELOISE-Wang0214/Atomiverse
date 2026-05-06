/// @file tests/cpu/test_optimize_driver_early_stop.cpp
/// @brief Phase 1.7 — OptimizeDriver breaks the loop when a
///        ConvergenceCheck-tagged component sets
///        `outer_status.converged = true`. Returns the iteration count
///        actually run.
///
/// Acceptance: ac-005 (Optimize half) in
/// `.claude/specs/cpu-driver-taxonomy.acceptance.md`.

#include <gtest/gtest.h>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/optimize_driver.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace {

struct StubOptimizer : ATV::cpu::ComponentBase {
    int step_count = 0;
    std::string name() const override { return "stub_optimizer"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Optimizer; }
    void launch(const ATV::cpu::LaunchContext&) override { ++step_count; }
};

struct StubEnergyEvaluator : ATV::cpu::ComponentBase {
    std::string name() const override { return "stub_energy"; }
    ATV::ComponentRole role() const override { return ATV::ComponentRole::Evaluator; }
    ATV::Capability   capabilities() const override { return ATV::Capability::Energy; }
    void launch(const ATV::cpu::LaunchContext&) override {}
};

/// Flips outer_status.converged to true on the iteration whose 0-based
/// step index equals `converge_at_step`. Uses `LaunchContext::step` so
/// the trigger lines up with the driver's iteration counter.
struct ConvergeAtStep : ATV::cpu::ComponentBase {
    int converge_at_step;
    int launch_count = 0;
    explicit ConvergeAtStep(int s) : converge_at_step(s) {}
    std::string name() const override { return "converge_at_step"; }
    ATV::ComponentRole role() const override {
        return ATV::ComponentRole::ConvergenceCheck;
    }
    void launch(const ATV::cpu::LaunchContext& ctx) override {
        ++launch_count;
        if (ctx.step >= converge_at_step) {
            ctx.runtime.outer_status.converged = true;
            ctx.runtime.outer_status.reason = "stub-converged";
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

TEST(OptimizeEarlyStop, BreaksAtConvergeStep) {
    using namespace ATV;
    cpu::OptimizeDriver opt;
    StubOptimizer optimizer;
    StubEnergyEvaluator energy;
    ConvergeAtStep converge{/*step=*/3};

    opt.on(minimizer::Phase::search_direction).add(optimizer);
    opt.on(minimizer::Phase::force).add(energy);
    opt.on(minimizer::Phase::convergence_check).add(converge);
    opt.build();

    auto frame = make_trivial_frame();
    int iters = opt.run(frame, /*max_iter=*/100);

    // 0-based step indices 0,1,2,3 — the convergence flips at step==3
    // and the driver breaks AFTER step 3 completes, so iters == 4.
    EXPECT_EQ(iters, 4) << "Expected 4 iterations (0,1,2,3 then break)";
    EXPECT_EQ(optimizer.step_count, 4);
    EXPECT_EQ(converge.launch_count, 4);
}

TEST(OptimizeEarlyStop, RunsToBudgetWhenNeverConverged) {
    using namespace ATV;
    cpu::OptimizeDriver opt;
    StubOptimizer optimizer;
    StubEnergyEvaluator energy;
    ConvergeAtStep never{/*step=*/9999};  // far beyond max_iter

    opt.on(minimizer::Phase::search_direction).add(optimizer);
    opt.on(minimizer::Phase::force).add(energy);
    opt.on(minimizer::Phase::convergence_check).add(never);
    opt.build();

    auto frame = make_trivial_frame();
    int iters = opt.run(frame, /*max_iter=*/7);
    EXPECT_EQ(iters, 7);
    EXPECT_EQ(optimizer.step_count, 7);
}
