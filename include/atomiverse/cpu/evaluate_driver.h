/// @file cpu/evaluate_driver.h
/// @brief cpu::EvaluateDriver — one-pass evaluator driver.
///
/// `EvaluateDriver` runs the registered evaluator components exactly once
/// against the current `cpu::Frame`. It is the outer-task category for
/// "compute one or more quantities at this geometry" — single-point
/// energy, force, density, dipole, charge, etc.
///
/// Validation rule (build-time): at least one component with
/// `role()==Evaluator` AND non-empty `capabilities()`. Energy-only,
/// property-only, and force-providing evaluators are all valid here.
///
/// `cpu::Evaluate` (in cpu/evaluate.h, Phase 5 task) is the concrete
/// driver users construct; this base provides the `run(frame)` machinery.
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/cpu_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

class EvaluateDriver : public CPUDriver<ATV::evaluate::Phase> {
  public:
    /// @brief Run a single evaluation pass through the compiled Plan.
    /// @param frame Host-resident system state, updated in place.
    void run(cpu::Frame& frame) {
        setup_pass(frame);
        step_once(0);
        teardown_pass();
    }

  protected:
    void validate_components() const override {
        bool found_evaluator = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (c.role() == ComponentRole::Evaluator
                && !is_none(c.capabilities())) {
                found_evaluator = true;
            }
        });
        if (!found_evaluator)
            throw_missing("Evaluator with non-empty capabilities()");
    }
};

} // namespace cpu
} // namespace ATV
