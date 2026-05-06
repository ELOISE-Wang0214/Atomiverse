/// @file cpu/dynamics_driver.h
/// @brief cpu::DynamicsDriver — time-advance driver.
///
/// `DynamicsDriver` runs the requested step budget. It does **not**
/// interpret `outer_status.converged` — time evolution doesn't have a
/// "convergence" success path. It does honour `outer_status.stop_requested`
/// (NaN detector, future StopCondition components, etc.) and breaks out
/// early when set.
///
/// Validation rule (build-time):
///   - ≥ 1 component with `role()==Integrator`
///   - ≥ 1 component with `role()==Evaluator` AND `capabilities()` carrying
///     `Capability::Force` (since dynamics needs forces to integrate).
///
/// Concrete user-facing class: `cpu::md::MD` (Phase 2 migrates its
/// parent from `Driver<md::Phase, cpu::device_tag>` to `DynamicsDriver`).
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/cpu_driver.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

class DynamicsDriver : public CPUDriver<ATV::md::Phase> {
public:
    /// @brief Advance the system for `n_steps` MD steps.
    /// Honours `outer_status.stop_requested` (early break); ignores
    /// `outer_status.converged` (Dynamics has no convergence success path).
    /// @param frame Host-resident system state, updated in place.
    /// @param n_steps Step budget.
    void run(cpu::Frame& frame, int n_steps) {
        setup_pass(frame);
        for (int s = 0; s < n_steps; ++s) {
            mutable_outer_status().reset();
            step_once(s);
            if (outer_status().stop_requested) break;
            // Dynamics deliberately ignores outer_status().converged.
        }
        teardown_pass();
    }

protected:
    void validate_components() const override {
        bool has_integrator = false;
        bool has_force_writer = false;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (c.role() == ComponentRole::Integrator)
                has_integrator = true;
            // Force capability is the load-bearing predicate — *any* component
            // (Evaluator-role LJ, Thermostat-role Langevin) that writes
            // forces satisfies Dynamics's "the integrator has something to
            // integrate" contract. Role is structural; capability is what
            // the integrator actually needs.
            if (has_capability(c.capabilities(), Capability::Force))
                has_force_writer = true;
        });
        if (!has_integrator)
            throw_missing("Integrator");
        if (!has_force_writer)
            throw_missing("Component with Capability::Force");
    }
};

} // namespace cpu
} // namespace ATV
