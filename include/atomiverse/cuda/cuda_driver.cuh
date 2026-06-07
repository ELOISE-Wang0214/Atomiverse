/// @file cuda/cuda_driver.cuh
/// @brief cuda::CUDADriver — the CUDA-backend layer of the driver hierarchy.
///
/// Hierarchy (two orthogonal axes — backend & category):
///
/// @code
///   Driver<PhaseEnum, DeviceTag>          ← plan / dispatch primitive
///     ├─ cpu::CPUDriver<PhaseEnum>         ← CPU backend layer
///     │    ├─ cpu::EvaluateDriver
///     │    ├─ cpu::DynamicsDriver
///     │    └─ cpu::OptimizeDriver
///     └─ cuda::CUDADriver<PhaseEnum>       ← CUDA backend layer (this file)
///          ├─ cuda::EvaluateDriver
///          ├─ cuda::DynamicsDriver
///          └─ cuda::OptimizeDriver
/// @endcode
///
/// `CUDADriver<P>` inherits the underlying `Driver<P, cuda::device_tag>`
/// (which holds the Plan and Backend ownership) and adds the
/// CUDA-specific machinery the three category bases share:
///
///   1. A `build()` override that runs `validate_components()` before the
///      Plan is compiled. Subclasses override `validate_components()` to
///      enforce their required role + capability set.
///   2. Protected helpers `setup_pass`, `step_once`, `teardown_pass` that
///      subclass `run(...)` implementations use to drive the
///      per-iteration loop with category-correct termination. The CUDA
///      flavour of `setup_pass` materializes the host frame onto the
///      device (H2D) and `teardown_pass` retrieves it back (D2H).
///   3. Read/write access to the outer-driver `IterationStatus`
///      (`cuda::Context::outer_status`), which `OptimizeDriver` polls for
///      `converged` and `DynamicsDriver` polls for `stop_requested`.
///
/// Unlike `cpu::CPUDriver`, this layer deliberately provides **no**
/// categorical-setter API (`set_pair`, `set_bond`, …) and no
/// `install_infra()` auto-injection — those are out of scope for the
/// CUDA driver-taxonomy chain.
///
/// User code does not subclass this directly — the three concrete CUDA
/// category bases (`cuda::EvaluateDriver` / `DynamicsDriver` /
/// `OptimizeDriver`) do.
#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/backend.cuh"
#include "atomiverse/cuda/backend_traits.cuh"
#include "atomiverse/cuda/component.cuh"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/iteration_status.h"

namespace ATV {
namespace cuda {

/// @brief CUDA-backend driver layer — shared by all three CUDA category
///        bases.
///
/// @tparam PhaseEnum Phase enum type for the driver category (must end
///         with a `COUNT` sentinel).
template <typename PhaseEnum>
class CUDADriver : public Driver<PhaseEnum, cuda::device_tag> {
  public:
    using Base = Driver<PhaseEnum, cuda::device_tag>;

    /// @brief Compile the Plan after validating the required role +
    ///        capability set.
    ///
    /// Installs no infrastructure (CUDA auto-injection is out of scope
    /// for this chain). Order: `validate_components()` first (so a build
    /// with a missing required component throws before the Plan is
    /// compiled), `Base::build()` second.
    void build() {
        validate_components();
        Base::build();
    }

  protected:
    /// @brief Subclass hook — validate the component set and throw on
    ///        missing roles / capabilities.
    ///
    /// Called *before* the Plan is built, so `this->pending_` still
    /// carries the registered components.
    virtual void validate_components() const = 0;

    /// @brief Walk every pending component binding.
    ///
    /// Components with multiple `default_bindings()` entries appear once
    /// per binding — duplicate-tolerant predicates are expected.
    /// @tparam Fn Callable invoked as `fn(const cuda::ComponentBase&)`.
    /// @param fn Visitor applied to each pending component.
    template <typename Fn>
    void for_each_pending(Fn&& fn) const {
        for (const auto& pb : this->pending_)
            fn(*pb.component);
    }

    /// @brief Mutable variant of for_each_pending — visitor may mutate the
    /// component (e.g. `install_infra_` gating per-step virial via
    /// `set_stress_required`). @tparam Fn Callable `fn(cuda::ComponentBase&)`.
    template <typename Fn>
    void for_each_pending_mut(Fn&& fn) {
        for (auto& pb : this->pending_)
            fn(*pb.component);
    }

    // ── Driver loop helpers (protected — subclasses only) ──────────

    /// @brief One-time setup — backend init, plan load, H2D materialize,
    ///        component build.
    ///
    /// The CUDA flavour transfers the host frame onto the device
    /// (`materialize`, H2D) before building components. Resets
    /// `outer_status` so a stale flag from a prior invocation cannot
    /// leak into this run.
    /// @param frame Host-resident system state, transferred to the device.
    void setup_pass(cpu::Frame& frame) {
        if (!this->backend_)
            this->backend_ = std::make_unique<typename Base::BackendT>();
        this->backend_->init();
        this->backend_->load(this->plan_);
        this->backend_->materialize(frame);
        this->backend_->build_components();
        mutable_outer_status().reset();
    }

    /// @brief Single iteration through the compiled Plan.
    /// @param iter_index Zero-based iteration index passed to components.
    void step_once(int iter_index) { this->backend_->step(iter_index); }

    /// @brief Tear down — destroy components, retrieve results (D2H),
    ///        shut down the backend.
    /// @param frame Host-resident frame populated from the device.
    void teardown_pass(cpu::Frame& frame) {
        this->backend_->destroy_components();
        this->backend_->retrieve(frame);
        this->backend_->shutdown();
    }

    /// @brief Read the outer-driver iteration status (set by
    ///        `ConvergenceCheck` components during the most recent
    ///        `step_once`).
    /// @return Const reference to the backend context's `outer_status`.
    const IterationStatus& outer_status() const {
        return this->backend_->context().outer_status;
    }

    /// @brief Mutable handle to the outer status (for per-iteration
    ///        `reset()`).
    /// @return Mutable reference to the backend context's `outer_status`.
    IterationStatus& mutable_outer_status() {
        return this->backend_->context().outer_status;
    }

    /// @brief Throw a clear exception when a required role or capability
    ///        is absent at `validate_components()` time.
    /// @param what Human-readable name of the missing requirement.
    [[noreturn]] static void throw_missing(const char* what) {
        throw std::runtime_error(std::string("driver build: missing required ")
                                 + what);
    }
};

} // namespace cuda
} // namespace ATV
