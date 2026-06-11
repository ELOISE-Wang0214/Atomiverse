/// @file cuda/driver/optimize_driver.cuh
/// @brief cuda::OptimizeDriver — CUDA energy-minimizer driver.
///
/// `OptimizeDriver` is the only current category that treats convergence
/// as the normal successful termination path. It iterates up to
/// `max_iter` times and breaks when `outer_status.converged` is set by a
/// `ConvergenceCheck`-tagged component (or when `stop_requested` is set
/// by error / interrupt).
///
/// Validation rule (build-time):
///   - ≥ 1 component with `role()==Optimizer`
///   - ≥ 1 component with `role()==Evaluator` AND `capabilities()`
///     carrying either `Capability::Energy` or `Capability::Force`.
///
/// Mirrors `cpu::OptimizeDriver`. The concrete CUDA `minimizer::Minimize`
/// driver is added in a later spec; this base provides the
/// `run(frame, max_iter)` machinery and the role/capability contract.
#pragma once

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/component/component.cuh"
#include "atomiverse/cuda/driver/cuda_driver.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {

/// @brief CUDA energy-minimizer category base.
class OptimizeDriver : public CUDADriver<ATV::minimizer::Phase> {
  public:
    /// @brief Iterate to convergence or up to `max_iter` iterations.
    /// @param frame Host-resident system state; transferred H2D, minimized
    ///        on the device, retrieved D2H, updated in place.
    /// @param max_iter Iteration budget [dimensionless iteration count].
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
        teardown_pass(frame);
        return iters;
    }

  protected:
    /// @brief Require at least one `Optimizer`-role component and at least
    ///        one `Evaluator`-role component carrying `Energy` or `Force`.
    void validate_components() const override {
        bool has_optimizer = false;
        bool has_energy_or_force_evaluator = false;
        for_each_pending([&](const cuda::ComponentBase& c) {
            if (c.role() == ComponentRole::Optimizer)
                has_optimizer = true;
            if (c.role() == ComponentRole::Evaluator
                && (has_capability(c.capabilities(), Capability::Energy)
                    || has_capability(c.capabilities(), Capability::Force)))
                has_energy_or_force_evaluator = true;
        });
        if (!has_optimizer)
            throw_missing("Optimizer");
        if (!has_energy_or_force_evaluator)
            throw_missing(
                "Evaluator with Capability::Energy or Capability::Force");
    }
};

} // namespace cuda
} // namespace ATV
