/// @file cpu/cpu_driver.h
/// @brief cpu::CPUDriver — the CPU-backend layer of the driver hierarchy.
///
/// Hierarchy (two orthogonal axes — backend & category):
///
/// @code
///   Driver<PhaseEnum, DeviceTag>          ← plan / dispatch primitive
///     ├─ cpu::CPUDriver<PhaseEnum>         ← CPU backend layer  (this file)
///     │    ├─ cpu::EvaluateDriver
///     │    ├─ cpu::DynamicsDriver
///     │    └─ cpu::OptimizeDriver
///     └─ cuda::CUDADriver<PhaseEnum>       ← CUDA backend layer (planned —
///          ├─ cuda::EvaluateDriver           cuda-driver-taxonomy.md)
///          ├─ cuda::DynamicsDriver
///          └─ cuda::OptimizeDriver
/// @endcode
///
/// `CPUDriver<P>` inherits the underlying `Driver<P, cpu::device_tag>`
/// (which holds the Plan and Backend ownership) and adds the
/// CPU-specific machinery the three category bases share:
///
///   1. A `build()` override that runs `validate_components()` before
///      the Plan is compiled. Subclasses override
///      `validate_components()` to enforce their required role +
///      capability set.
///   2. Protected helpers `setup_pass`, `step_once`, `teardown_pass`
///      that subclass `run(...)` implementations use to drive the
///      per-iteration loop with category-correct termination.
///   3. Read/write access to the outer-driver `IterationStatus`
///      (`cpu::Context::outer_status`), which `OptimizeDriver` polls
///      for `converged` and `DynamicsDriver` polls for
///      `stop_requested`.
///
/// User code does not subclass this directly — the three concrete CPU
/// category bases (`cpu::EvaluateDriver` / `DynamicsDriver` /
/// `OptimizeDriver`) do.
#pragma once

#include <stdexcept>
#include <string>

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/iteration_status.h"

namespace ATV {
namespace cpu {

/// @brief CPU-backend driver layer — shared by all three CPU category bases.
template <typename PhaseEnum>
class CPUDriver : public Driver<PhaseEnum, cpu::device_tag> {
public:
    using Base = Driver<PhaseEnum, cpu::device_tag>;

    /// @brief Compile the Plan after validating required role + capability set.
    void build() {
        validate_components();
        Base::build();
    }

protected:
    /// @brief Subclass hook — validate the component set and throw on
    /// missing roles / capabilities. Called *before* the Plan is built,
    /// so `this->pending_` still carries the registered components.
    virtual void validate_components() const = 0;

    /// @brief Walk every (component, phase, order) pending binding.
    /// Components with multiple `default_bindings()` entries appear once
    /// per binding — duplicate-tolerant predicates expected.
    template <typename Fn>
    void for_each_pending(Fn&& fn) const {
        for (const auto& pb : this->pending_) {
            fn(*pb.component);
        }
    }

    /// @brief One-time setup — backend init, plan load, set frame, build components.
    /// Resets outer_status before the run so a stale flag from a prior
    /// invocation cannot leak into this one.
    void setup_pass(cpu::Frame& frame) {
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->init();
        this->backend_->load(this->plan_);
        this->backend_->set_frame(frame);
        this->backend_->build_components();
        mutable_outer_status().reset();
    }

    /// @brief Single iteration through the compiled Plan.
    void step_once(int iter_index) {
        this->backend_->step(iter_index);
    }

    /// @brief Tear down — destroy components + shut down backend.
    void teardown_pass() {
        this->backend_->destroy_components();
        this->backend_->shutdown();
    }

    /// @brief Read the outer-driver iteration status (set by ConvergenceCheck
    /// components during the most recent `step_once`).
    const IterationStatus& outer_status() const {
        return this->backend_->context().outer_status;
    }

    /// @brief Mutable handle to outer status (for per-iteration `reset()`).
    IterationStatus& mutable_outer_status() {
        return this->backend_->context().outer_status;
    }

    /// @brief Convenience: throw a clear exception when a required role or
    /// capability is absent at `validate_components()` time.
    [[noreturn]] static void throw_missing(const char* what) {
        throw std::runtime_error(std::string("driver build: missing required ") + what);
    }
};

} // namespace cpu
} // namespace ATV

// Explicit Driver<> instantiations live in src/cpu/driver.cpp alongside
// the underlying Driver<> template instantiations.
