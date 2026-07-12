/// @file cpu/driver/optimize_driver.h
/// @brief cpu::OptimizeDriver — energy-minimizer driver.
///
/// `OptimizeDriver` is the only current category that treats convergence
/// as the normal successful termination path. It iterates up to
/// `max_iter` times and breaks when `outer_status.converged` is set —
/// in practice by the Optimizer component itself during the
/// `convergence_check` phase (LBFGS / CG own their stopping criteria;
/// there is no separate ConvergenceCheck-tagged component) — or when
/// `stop_requested` is set by error / interrupt.
///
/// Unlike `DynamicsDriver::run` (continuable: consecutive calls extend
/// one trajectory), each `run()` here is a complete, self-contained
/// optimization: `setup_pass` re-runs component `build()`, wiping the
/// optimizer's curvature history and step counter, and `teardown_pass`
/// runs before returning. That is deliberate — a converged geometry has
/// no useful "continuation state": re-running means re-optimizing
/// (possibly after the caller perturbed the system), so stale L-BFGS
/// history would be wrong to keep.
///
/// Validation rule (build-time):
///   - ≥ 1 component with `role()==Optimizer`
///   - ≥ 1 component with `role()==Evaluator` AND `capabilities()` carrying
///     either `Capability::Energy` or `Capability::Force`.
///
/// Concrete user-facing class: `cpu::minimizer::Minimize`.
#pragma once

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/driver/cpu_driver.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/frame/materialize.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

class OptimizeDriver : public CPUDriver<ATV::minimizer::Phase> {
  public:
    /// @brief Iterate to convergence or up to `max_iter` iterations.
    /// @param frame Host-resident system state, updated in place.
    /// @param max_iter Iteration budget.
    /// @return Number of iterations actually executed (≤ `max_iter`).
    int run(cpu::Frame& frame, int max_iter) {
        setup_pass(frame);
        // Publish the re-entrant force re-eval handle so an Optimizer
        // component can evaluate trial geometries during the line_search
        // phase. Replays the force phase (zeroing → evaluator → reduce) on
        // the frame's current geometry. Captures `this`; valid for the run.
        this->backend_->context().opt.reeval_forces =
            [this](bool /*need_grad*/) {
                this->backend_->run_stage(Traits::eval_phase, /*step_idx=*/0);
            };
        int iters = 0;
        for (int s = 0; s < max_iter; ++s) {
            mutable_outer_status().reset();
            step_once(s);
            ++iters;
            // The optimizer PUBLISHED its metrics during the step; the handlers judge
            // them here. Neither LBFGS nor CG decides its own convergence any more —
            // which is why a nested optimizer can no longer stop this loop.
            const Verdict v = poll_handlers_(s);
            if (v != Verdict::Continue)
                break;
        }
        teardown_pass();
        return iters;
    }

    /// @brief User-facing entry — optimize the system stored in a molrs.Frame.
    ///
    /// Materializes the molrs.Frame into a transient `cpu::Frame`, runs the
    /// minimization loop on it, and writes the optimized geometry back into
    /// the *same* molrs.Frame. The borrow window is strictly this one call.
    ///
    /// @param frame molrs.Frame handle (input geometry + output sink).
    /// @param max_iter Iteration budget.
    /// @return Number of iterations actually executed (≤ `max_iter`).
    int run(molrs::FrameRef& frame, int max_iter) {
        cpu::Frame work = cpu::materialize(frame);
        const int iters = run(work, max_iter);
        cpu::write_back(work, frame);
        return iters;
    }

  protected:
    void validate_components() const override {
        bool has_optimizer = false;
        bool has_energy_or_force_evaluator = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (c.role() == ComponentRole::Optimizer)
                has_optimizer = true;
            if (c.role() == ComponentRole::Evaluator
                && (has_capability(c.capabilities(), Capability::Energy)
                    || has_capability(c.capabilities(), Capability::Force))) {
                has_energy_or_force_evaluator = true;
            }
        });
        if (!has_optimizer)
            throw_missing("Optimizer");
        if (!has_energy_or_force_evaluator)
            throw_missing(
                "Evaluator with Capability::Energy or Capability::Force");
    }
};

} // namespace cpu
} // namespace ATV
