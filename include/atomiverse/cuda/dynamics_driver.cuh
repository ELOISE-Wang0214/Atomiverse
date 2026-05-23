/// @file cuda/dynamics_driver.cuh
/// @brief cuda::DynamicsDriver — CUDA time-advance driver.
///
/// `DynamicsDriver` runs the requested step budget. It does **not**
/// interpret `outer_status.converged` — time evolution doesn't have a
/// "convergence" success path. It does honour
/// `outer_status.stop_requested` (NaN detector, future StopCondition
/// components, etc.) and breaks out early when set.
///
/// Validation rule (build-time):
///   - ≥ 1 component with `role()==Integrator`
///   - ≥ 1 component with `capabilities()` carrying `Capability::Force`
///     (since dynamics needs forces to integrate).
///
/// Mirrors `cpu::DynamicsDriver`. The concrete CUDA `md::MD` driver is
/// reparented onto this base in a later spec; this base provides the
/// `run(frame, n_steps)` machinery and the role/capability contract.
#pragma once

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/cuda_driver.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {

/// @brief CUDA time-advance category base.
class DynamicsDriver : public CUDADriver<ATV::md::Phase> {
  public:
    /// @brief Advance the system for `n_steps` MD steps.
    ///
    /// **Continuable.** Consecutive `run()` calls on the same `frame`
    /// extend one trajectory — components are built once on the first
    /// call, and the step index passed to components is cumulative
    /// across calls. This matches `cpu::DynamicsDriver` so that:
    ///   - `FixLangevin`'s per-step RNG keying (which mixes the
    ///     cumulative `ctx.step` into `base_seed`) does NOT repeat the
    ///     same z-stream every `run()` call — a critical FDT-correctness
    ///     requirement; if `s` reset to 0 each call, the Gaussian
    ///     impulse sequence repeats with period `n_steps`, breaking
    ///     equipartition (T_eq diverges to ~3.5× target observed in the
    ///     `Thermalises` test).
    ///   - `NoseHoover`/`NoseHooverNPT` extended variables persist.
    ///   - Recorder cadence keys (modulo `frequency`) stay continuous.
    /// A different `frame` ⇒ genuinely independent simulation: tear the
    /// previous run down and restart the counter.
    ///
    /// Honours `outer_status.stop_requested` (early break); ignores
    /// `outer_status.converged` (Dynamics has no convergence success
    /// path).
    /// @param frame Host-resident system state; transferred H2D, advanced
    ///        on the device, retrieved D2H, updated in place.
    /// @param n_steps Step budget [dimensionless step count].
    void run(cpu::Frame& frame, int n_steps) {
        if (live_ && frame_ptr_ != &frame) {
            teardown_pass(*frame_ptr_);
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
        // Pull device state back to the host frame so the caller can
        // inspect intermediate observables between chunks (matches the
        // pre-continuation behaviour of `teardown_pass` for the result
        // retrieval portion). Components remain alive on the device for
        // subsequent `run()` calls.
        this->backend_->retrieve(frame);
    }

    ~DynamicsDriver() {
        if (live_ && frame_ptr_)
            teardown_pass(*frame_ptr_);
    }

  protected:
    /// @brief Require at least one `Integrator`-role component and at
    ///        least one component carrying `Capability::Force`.
    void validate_components() const override {
        bool has_integrator = false;
        bool has_force_writer = false;
        for_each_pending([&](const cuda::ComponentBase& c) {
            if (c.role() == ComponentRole::Integrator)
                has_integrator = true;
            // Force capability is the load-bearing predicate — any
            // component that writes forces satisfies Dynamics's "the
            // integrator has something to integrate" contract.
            if (has_capability(c.capabilities(), Capability::Force))
                has_force_writer = true;
        });
        if (!has_integrator)
            throw_missing("Integrator");
        if (!has_force_writer)
            throw_missing("Component with Capability::Force");
    }

  private:
    /// @brief True between first `run()` and driver destruction —
    /// components are live on the device, must not be rebuilt.
    bool live_ = false;
    /// @brief Frame bound by the live run; a different frame restarts.
    cpu::Frame* frame_ptr_ = nullptr;
    /// @brief Cumulative MD step index across all `run()` calls.
    int step_counter_ = 0;
};

} // namespace cuda
} // namespace ATV
