/// @file cuda/evaluate_driver.cuh
/// @brief cuda::EvaluateDriver — one-pass CUDA evaluator driver.
///
/// `EvaluateDriver` runs the registered evaluator components exactly once
/// against the current device-resident frame. It is the CUDA outer-task
/// category for "compute one or more quantities at this geometry" —
/// single-point energy, force, density, dipole, charge, etc.
///
/// Validation rule (build-time): at least one component with
/// `role()==Evaluator` AND non-empty `capabilities()`. Energy-only,
/// property-only, and force-providing evaluators are all valid here.
///
/// Mirrors `cpu::EvaluateDriver`. The concrete CUDA `Evaluate` driver
/// users construct is added in a later spec; this base provides the
/// `run(frame)` machinery and the role/capability contract.
#pragma once

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/cuda_driver.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {

/// @brief CUDA one-pass evaluator category base.
class EvaluateDriver : public CUDADriver<ATV::evaluate::Phase> {
  public:
    /// @brief Run a single evaluation pass through the compiled Plan.
    /// @param frame Host-resident system state; transferred H2D, evaluated
    ///        on the device, retrieved D2H, updated in place.
    void run(cpu::Frame& frame) {
        setup_pass(frame);
        step_once(0);
        teardown_pass(frame);
    }

  protected:
    /// @brief Require at least one `Evaluator`-role component with a
    ///        non-empty `capabilities()` bitset.
    void validate_components() const override {
        bool found_evaluator = false;
        for_each_pending([&](const cuda::ComponentBase& c) {
            if (c.role() == ComponentRole::Evaluator
                && !is_none(c.capabilities()))
                found_evaluator = true;
        });
        if (!found_evaluator)
            throw_missing("Evaluator with non-empty capabilities()");
    }
};

} // namespace cuda
} // namespace ATV
