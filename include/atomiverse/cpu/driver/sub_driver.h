/// @file cpu/driver/sub_driver.h
/// @brief SubDriver — a Driver, in a Plan slot. Nesting, in one class.
///
/// ## "Can a driver call a driver?"
///
/// The question dissolves. A Driver has the same lifecycle as a component —
/// `build()` → `launch()` × N → `destroy()` — so it simply IS a node, and a Plan slot
/// holds a component or a sub-driver indifferently. `SubDriver` is the thin adapter that
/// says so to the type system.
///
/// @code
///   cpu::minimizer::Minimize inner;          // its own phases, handlers, optimizer
///   inner.set_optimizer(cpu::minimizer::LBFGS{...});
///   inner.on(minimizer::Phase::force).add(forces);
///
///   cpu::md::MD outer;                       // a different phase enum entirely
///   outer.set_integrator(...);
///   outer.on(md::Phase::end_of_step).add(cpu::SubDriver{inner, /*budget=*/50});
///   outer.build();                           // builds BOTH; inner allocates once
///   outer.run(frame, 1000);                  // 1000 outer steps, each running inner
/// @endcode
///
/// ## Type erasure sits at the NODE boundary, not the stage boundary
///
/// The child above loops over `minimizer::Phase`; the parent loops over `md::Phase`. They
/// never have to agree, and neither ever sees the other's enum. What crosses the boundary
/// is a `ComponentBase*` — a node — and stages stay typed enums all the way down. This is
/// the whole reason `LoopHandler` is not a template and `slot_of` is the only cast in the
/// codebase.
///
/// ## One Runtime
///
/// The child runs its Plan on the PARENT's Backend, reached through `ctx.exec.runtime`.
/// No second Backend, no second Context, no rebuilt domain decomposition, no re-bound
/// frame, no split communicator. Nesting costs a loop, not an engine.
///
/// ## And the inner loop cannot stop the outer one
///
/// `Driver::run_nested` re-points `ctx.loop` at the CHILD's `IterationStatus` for the
/// duration of the inner loop and restores the parent's on exit. When the child
/// converges, the verdict lands on the child's status — which the parent's loop never
/// reads. The bug that shipped twice is now structurally unwritable.
///
/// ## Known ceiling
///
/// **Components cannot be cloned.** Running the same inner Driver over N replicas would
/// need N independent component sets (N L-BFGS curvature histories), and there is no
/// `clone()` on `ComponentBase`. Nesting makes that problem SHARPER, not softer — it does
/// not solve it, and this spec deliberately does not decide it. Today a `SubDriver` holds
/// ONE inner driver with ONE component set. @see .claude/NOTES.md
#pragma once

#include <string>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/backend.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Adapter that puts a Driver into another Driver's Plan slot.
///
/// Non-owning: the caller owns the inner driver and may inspect it after the run
/// (`inner.converged()`, `inner.stop_reason()`), which is the same lvalue-registration
/// contract `Driver::add()` already offers for components.
///
/// @tparam InnerPhase The child's Phase enum — distinct from the parent's, and never
///         visible to it.
template <typename InnerPhase>
class SubDriver : public cpu::ComponentBase {
  public:
    /// @brief Wrap an inner driver as a node.
    /// @param inner The child loop. Borrowed; must outlive the parent run.
    /// @param budget Iteration budget for ONE run of the child, per outer iteration.
    /// @param label Name reported in plans and timing summaries.
    SubDriver(CPUDriver<InnerPhase>& inner, int budget,
              std::string label = "sub_driver")
        : inner_(&inner), budget_(budget), label_(std::move(label)) {}

    std::string name() const override { return label_; }

    /// @brief Structural role: none. A sub-driver is not an Integrator or an Evaluator —
    ///        it is a loop. What it CONTAINS has roles; the container does not.
    ComponentRole role() const override { return ComponentRole::None; }

    /// @brief Compile the child's Plan and build its components — ONCE, here, at the
    ///        parent's `build()`.
    ///
    /// Nothing is allocated per outer iteration: `launch()` only loops.
    void build(const cpu::BuildContext& ctx) override {
        inner_->compile(); // the child's Plan: install_defaults + validate + compile
        ctx.runtime.exec.runtime->build_plan(inner_->plan());
        built_ = true;
    }

    /// @brief Run the child's whole loop, once, as one operation of the parent's.
    void launch(const cpu::LaunchContext& ctx) override {
        last_iterations_ = inner_->run_nested(ctx.runtime, ctx.frame, budget_);
    }

    /// @brief Destroy the child's components, in reverse Plan order.
    void destroy(const cpu::DestroyContext& ctx) override {
        if (!built_)
            return;
        ctx.runtime.exec.runtime->destroy_plan(inner_->plan());
        built_ = false;
    }

    /// @brief Iterations the child ran on its most recent `launch()`.
    int last_iterations() const { return last_iterations_; }

  private:
    CPUDriver<InnerPhase>* inner_; ///< The child loop (non-owning).
    int budget_;                   ///< Child iteration budget per outer iteration.
    std::string label_;            ///< Name in plans / timing.
    int last_iterations_ = 0;      ///< Iterations the child last ran.
    bool built_ = false;           ///< Guards a destroy() without a build().
};

/// @brief Deduce the child's Phase enum from the driver — `cpu::SubDriver{inner, 50}`.
template <typename InnerPhase>
SubDriver(CPUDriver<InnerPhase>&, int) -> SubDriver<InnerPhase>;

/// @brief Deduction guide for the labelled form.
template <typename InnerPhase>
SubDriver(CPUDriver<InnerPhase>&, int, std::string) -> SubDriver<InnerPhase>;

} // namespace cpu
} // namespace ATV
