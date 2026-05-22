/// @file cpu/optimize_driver.h
/// @brief cpu::OptimizeDriver — energy-minimizer driver.
///
/// `OptimizeDriver` is the only current category that treats convergence
/// as the normal successful termination path. It iterates up to
/// `max_iter` times and breaks when `outer_status.converged` is set by
/// a `ConvergenceCheck`-tagged component (or when `stop_requested` is
/// set by error / interrupt).
///
/// Validation rule (build-time):
///   - ≥ 1 component with `role()==Optimizer`
///   - ≥ 1 component with `role()==Evaluator` AND `capabilities()` carrying
///     either `Capability::Energy` or `Capability::Force`.
///
/// Concrete user-facing class: `cpu::minimizer::Minimize` (Phase 2
/// renames the existing `cpu::minimizer::LBFGS` driver to this and
/// extracts LBFGS / CG / FIRE as `Optimizer` components).
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/cpu_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/materialize.h"
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
        int iters = 0;
        for (int s = 0; s < max_iter; ++s) {
            mutable_outer_status().reset();
            step_once(s);
            ++iters;
            if (outer_status().converged || outer_status().stop_requested)
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
