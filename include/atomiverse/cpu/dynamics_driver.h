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
#include "atomiverse/cpu/cpu_driver.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

class DynamicsDriver : public CPUDriver<ATV::md::Phase> {
  public:
    /// @brief Flush component `destroy()` hooks if a run is still live.
    ///
    /// `run()` is continuable: `setup_pass()` (which calls component
    /// `build()`) fires exactly once, on the first `run()`, and the
    /// matching `destroy()` pass is deferred to here so the component
    /// `build → launch×N → destroy` lifecycle spans *all* `run()` calls.
    /// The base `~Driver()` shuts the backend down afterwards.
    ~DynamicsDriver() {
        if (live_ && this->backend_)
            this->backend_->destroy_components();
    }

    /// @brief Advance the system for `n_steps` MD steps.
    ///
    /// **Continuable.** Consecutive calls extend one trajectory — the
    /// standard `run(burn_in); run(sample); …` pattern. Component
    /// `build()` runs once on the first call only; rebuilding it every
    /// call would wipe the extended-variable state of stateful
    /// integrators (`NoseHooverNPT` / `NoseHoover` NHC chains,
    /// `FixLangevin` RNG, …) and re-randomise `InitVelocity`, silently
    /// corrupting any chunked run. The step index passed to components
    /// is cumulative across calls so frequency-gated components
    /// (`Thermo`, `AtomMigration`, recorders) keep a continuous cadence.
    ///
    /// Honours `outer_status.stop_requested` (early break); ignores
    /// `outer_status.converged` (Dynamics has no convergence success path).
    /// @param frame Host-resident system state, updated in place.
    /// @param n_steps Step budget.
    void run(cpu::Frame& frame, int n_steps) {
        // A different frame ⇒ a genuinely independent simulation: tear
        // the previous run down and start fresh. Same frame ⇒ continue.
        if (live_ && frame_ptr_ != &frame) {
            this->backend_->destroy_components();
            live_ = false;
        }
        if (!live_) {
            setup_pass(frame);
            frame_ptr_ = &frame;
            step_counter_ = 0;
            live_ = true;
        }
        for (int s = 0; s < n_steps; ++s) {
            mutable_outer_status().reset();
            step_once(step_counter_++);
            if (outer_status().stop_requested)
                break;
            // Dynamics deliberately ignores outer_status().converged.
        }
        // Drain async-analysis workers so every snapshot from this
        // segment reaches its subscribers before `run()` returns — the
        // per-call teardown used to do this; continuation defers full
        // teardown, so the barrier is made explicit here.
        this->backend_->drain_analysis();
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

  private:
    /// @brief True between the first `run()` and driver destruction —
    /// the backend + components are live and must not be rebuilt.
    bool live_ = false;
    /// @brief Frame bound by the live run; a different frame restarts.
    cpu::Frame* frame_ptr_ = nullptr;
    /// @brief Cumulative MD step index across all `run()` calls.
    int step_counter_ = 0;
};

} // namespace cpu
} // namespace ATV
