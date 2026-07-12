/// @file tests/cpu/driver/test_sub_driver.cpp
/// @brief A Driver is an executable node — nesting works, and an inner loop cannot stop
///        its parent.
///
/// Acceptance: ac-003 (status belongs to the loop) and ac-005 (a Driver is a node) in
/// `.claude/specs/driver-handler.acceptance.md`.
///
/// The bug being pinned here shipped TWICE:
///   - probe 2: `LBFGS` wrote `outer_status.converged = true` unconditionally, because
///     the optimizer owned the stopping criterion.
///   - probe 3: so the first NEB image to converge would have stopped the entire band,
///     and `ReplicaSet` had to grow a private convergence table to work around a status
///     it could not trust.
/// Both had the same cause: ONE `IterationStatus`, living on the Context, reachable by
/// anybody. Now each loop owns its own and the Context merely points at the active one,
/// so `NestedConvergenceDoesNotStopTheParent` is not testing a rule — it is testing that
/// there is no longer a pointer from the child to the parent's status.

#include <gtest/gtest.h>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/driver/sub_driver.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/facade/stage.h"

#include "cpu/driver_test_helpers.h"

using namespace ATV;

namespace {

/// @brief A minimal loop over md::Phase — the PARENT. Counts its own iterations.
class OuterLoop : public cpu::CPUDriver<ATV::md::Phase> {
  public:
    int run(cpu::Frame& frame, int budget) {
        setup_pass(frame);
        int iters = 0;
        for (int s = 0; s < budget; ++s) {
            mutable_outer_status().reset();
            step_once(s);
            ++iters;
            if (poll_handlers_(s) != Verdict::Continue)
                break;
        }
        teardown_pass();
        return iters;
    }

  protected:
    void validate_components() const override {} // a bare loop requires nothing
};

/// @brief A minimal loop over minimizer::Phase — the CHILD. A DIFFERENT phase enum, which
///        the parent never sees: erasure is at the node boundary, not the stage boundary.
class InnerLoop : public cpu::CPUDriver<ATV::minimizer::Phase> {
  protected:
    void validate_components() const override {}
};

/// @brief Counts build/launch/destroy, and reports convergence on a chosen iteration.
///
/// It writes into `ctx.runtime.loop.metrics` — the ONLY channel a component has. It
/// cannot write `converged`: there is no such field reachable from here, which is the
/// point. A handler turns the metric into the verdict.
class ConvergeAt : public cpu::ComponentBase {
  public:
    explicit ConvergeAt(int at) : at_(at) {}
    std::string name() const override { return "converge_at"; }
    void build(const cpu::BuildContext&) override { ++builds; }
    void destroy(const cpu::DestroyContext&) override { ++destroys; }
    void launch(const cpu::LaunchContext& ctx) override {
        // Report a huge force until the target iteration, then a zero one. ForceConverged
        // judges the metric; this component never judges anything.
        const bool done = (calls >= at_);
        ++calls;
        if (auto* m = ctx.runtime.loop.metrics) {
            m->max_force = done ? 0.0 : 1e3;
            m->rms_force = done ? 0.0 : 1e3;
            m->energy_change = done ? 0.0 : 1e3;
            m->has_displacement = false;
        }
    }
    int builds = 0;
    int destroys = 0;
    int calls = 0;

  private:
    int at_;
};

/// @brief Publishes a fixed force metric every iteration. The PARENT's measurement.
class PublishForce : public cpu::ComponentBase {
  public:
    explicit PublishForce(double f) : f_(f) {}
    std::string name() const override { return "publish_force"; }
    void launch(const cpu::LaunchContext& ctx) override {
        if (auto* m = ctx.runtime.loop.metrics) {
            m->max_force = f_;
            m->rms_force = f_;
            m->energy_change = f_;
            m->has_displacement = false;
        }
    }

  private:
    double f_;
};

/// @brief Records the Backend address it is launched with, and the DD pointer.
class RuntimeProbe : public cpu::ComponentBase {
  public:
    std::string name() const override { return "runtime_probe"; }
    void build(const cpu::BuildContext& ctx) override {
        ++builds;
        build_runtime = ctx.runtime.exec.runtime;
    }
    void launch(const cpu::LaunchContext& ctx) override {
        ++calls;
        launch_runtime = ctx.runtime.exec.runtime;
        launch_dd = ctx.runtime.md.dd;
        frame_addr = &ctx.frame;
    }
    int builds = 0;
    int calls = 0;
    cpu::Backend* build_runtime = nullptr;
    cpu::Backend* launch_runtime = nullptr;
    const void* launch_dd = nullptr;
    const void* frame_addr = nullptr;
};

cpu::Frame two_atoms() {
    cpu::Frame f(2);
    for (int i = 0; i < 2; ++i) {
        f.atoms.x[i] = 0.5 * i;
        f.atoms.y[i] = 0.0;
        f.atoms.z[i] = 0.0;
        f.atoms.mass[i] = 1.0;
        f.atoms.type[i] = 1;
    }
    f.box = Mat3::cubic(10.0);
    return f;
}

} // namespace

// ── ac-003: the inner loop's convergence does NOT stop the outer loop ─────────

TEST(SubDriver, NestedConvergenceDoesNotStopTheParent) {
    cpu::Frame frame = two_atoms();

    // The CHILD converges: from its 3rd inner iteration it publishes a zero force, and
    // its own ForceConverged handler stops it.
    InnerLoop inner;
    ConvergeAt converge{/*at=*/2};
    inner.on(ATV::minimizer::Phase::convergence_check).add(converge);
    inner.add_handler(cpu::ForceConverged{});

    // The PARENT does NOT converge: it publishes a huge force every iteration. And,
    // crucially, it HAS a handler — so it is a loop that *could* be stopped. Without a
    // handler the parent could never stop for any reason and this test would prove
    // nothing.
    OuterLoop outer;
    PublishForce big{1e3};
    cpu::SubDriver node{inner, /*budget=*/10};
    outer.on(ATV::md::Phase::force).add(big);
    outer.on(ATV::md::Phase::end_of_step).add(node); // the child runs AFTER the parent
    outer.add_handler(cpu::ForceConverged{});        //   measures — so it writes last
    outer.build();

    const int outer_iters = outer.run(frame, /*budget=*/5);

    // THE ASSERTION, and it bites. The child publishes `max_force = 0` into
    // `ctx.loop.metrics`, in the end_of_step phase — i.e. it is the LAST thing written
    // before the parent polls its handlers. If the Context still pointed at ONE shared
    // LoopMetrics, those zeros would land on the parent's metrics, the parent's
    // ForceConverged would read `max_force == 0`, and the outer loop would stop after
    // iteration 1 having moved nothing.
    //
    // It runs all 5. `run_nested` re-points `ctx.loop` at the child's own status and
    // metrics for the duration of the inner loop and restores the parent's on exit, so
    // the child's measurements are the child's. This is probe 2 and probe 3's bug, and
    // it is now unwritable rather than merely forbidden.
    EXPECT_EQ(outer_iters, 5) << "the inner loop stopped the outer loop — the bug that "
                                 "shipped in probe 2 and again in probe 3";
    EXPECT_FALSE(outer.converged())
        << "the parent adopted its child's verdict; status is not per-loop";

    // ...and the child really did converge, on its own status. So the test is not
    // passing merely because nothing converged anywhere.
    EXPECT_TRUE(inner.converged());
    EXPECT_EQ(inner.stop_reason(), "force_converged");
}

// ── ac-005: one Runtime, one build, nothing allocated per outer iteration ─────

TEST(SubDriver, RunsOnTheParentsRuntime) {
    cpu::Frame frame = two_atoms();

    InnerLoop inner;
    RuntimeProbe inner_probe;
    inner.on(ATV::minimizer::Phase::force).add(inner_probe);

    OuterLoop outer;
    RuntimeProbe outer_probe;
    cpu::SubDriver node{inner, /*budget=*/3};
    outer.on(ATV::md::Phase::force).add(outer_probe);
    outer.on(ATV::md::Phase::end_of_step).add(node);
    outer.build();

    outer.run(frame, /*budget=*/4);

    // ONE Runtime. No second Backend was constructed for the child: it executed its own
    // Plan on its parent's. Same Context, same domain decomposition, same frame.
    EXPECT_NE(inner_probe.launch_runtime, nullptr);
    EXPECT_EQ(inner_probe.launch_runtime, outer_probe.launch_runtime);
    EXPECT_EQ(inner_probe.launch_dd, outer_probe.launch_dd);
    EXPECT_EQ(inner_probe.frame_addr, outer_probe.frame_addr);
    EXPECT_EQ(inner_probe.build_runtime, outer_probe.build_runtime);
}

TEST(SubDriver, InnerBuildRunsOnceAtOuterBuild) {
    cpu::Frame frame = two_atoms();

    InnerLoop inner;
    ConvergeAt inner_comp{/*at=*/999}; // never converges; always runs the full budget
    inner.on(ATV::minimizer::Phase::convergence_check).add(inner_comp);

    OuterLoop outer;
    cpu::SubDriver node{inner, /*budget=*/3};
    outer.on(ATV::md::Phase::end_of_step).add(node);

    outer.build();
    EXPECT_EQ(inner_comp.builds, 0) << "build() compiles Plans; it must not run them";

    outer.run(frame, /*budget=*/4);

    // Built ONCE, at the outer setup — not once per outer iteration. This is the
    // "allocates nothing per outer iteration" contract: component build() is where
    // allocation happens, and it happened once.
    EXPECT_EQ(inner_comp.builds, 1);
    EXPECT_EQ(inner_comp.destroys, 1);

    // ...while launch() ran 4 outer iterations x 3 inner = 12 times.
    EXPECT_EQ(inner_comp.calls, 12);
    EXPECT_EQ(node.last_iterations(), 3);
}

// ── The child keeps its own phase enum ────────────────────────────────────────

TEST(SubDriver, ChildKeepsItsOwnPhaseEnum) {
    // The parent loops over md::Phase; the child over minimizer::Phase. Neither template
    // is instantiated with the other's enum. Type erasure happens where the node crosses
    // into the Plan (as a ComponentBase*), never at the stage boundary — stages stay
    // typed enums all the way down.
    static_assert(
        std::is_base_of_v<cpu::CPUDriver<ATV::md::Phase>, OuterLoop>,
        "parent loops over md::Phase");
    static_assert(
        std::is_base_of_v<cpu::CPUDriver<ATV::minimizer::Phase>, InnerLoop>,
        "child loops over minimizer::Phase");
    static_assert(
        std::is_base_of_v<cpu::ComponentBase, cpu::SubDriver<ATV::minimizer::Phase>>,
        "a SubDriver IS a node — that is what lets a Plan slot hold a driver");
    SUCCEED();
}
