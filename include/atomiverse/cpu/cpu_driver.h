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

#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "atomiverse/core/device.h"
#include "atomiverse/core/real.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/pair_evaluator.h"
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/iteration_status.h"

namespace ATV {
namespace cpu {

/// @brief CPU-backend driver layer — shared by all three CPU category bases.
template <typename PhaseEnum>
class CPUDriver : public Driver<PhaseEnum, cpu::device_tag> {
public:
    using Base = Driver<PhaseEnum, cpu::device_tag>;

    /// @brief Compile the Plan after installing infrastructure and
    /// validating the required role + capability set.
    ///
    /// Order: `install_infra()` first (so auto-injected components are
    /// visible to validation and the compiled Plan), `validate_components()`
    /// second (so the auto-injected set counts toward role/capability
    /// requirements), `Base::build()` last.
    void build() {
        install_infra();
        validate_components();
        Base::build();
    }

protected:
    /// @brief Subclass hook — install opt-in infrastructure components
    /// derived from declared user state (e.g. `NoseHooverNPT::reads()
    /// == Stress` → inject `VirialStress`; pair evaluator + multi-rank
    /// build → inject MPI halo exchange). Default no-op. Subclasses
    /// override to apply category-specific rules. Called from `build()`
    /// before `validate_components()` so the auto-injected set counts
    /// toward validation.
    /// @see .claude/specs/driver-infra-auto-injection.md
    virtual void install_infra() {}

    /// @brief Subclass hook — validate the component set and throw on
    /// missing roles / capabilities. Called *before* the Plan is built,
    /// so `this->pending_` still carries the registered components
    /// (including any auto-injected by `install_infra()`).
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

    // ─── Categorical-setter helpers (shared across MD / Minimize / Evaluate) ───

    /// @brief Singular setter using the component's `default_bindings()`.
    /// Idempotent on same-name re-call; throws on different-name re-call.
    template <typename C>
    void set_singular_default_(C&& comp,
                               std::optional<std::string>& slot,
                               const char* setter_name) {
        const std::string n = comp.name();
        if (slot) {
            if (*slot == n) return;
            throw std::invalid_argument(
                std::string(setter_name) + " already configured with '" +
                *slot + "', rejected re-call with '" + n + "'");
        }
        slot = n;
        this->add(std::forward<C>(comp));
    }

    /// @brief Singular setter at an explicit phase / order. Used for
    /// components without `default_bindings()` (NeighborList) or where
    /// the categorical API pins a specific slot (`set_pair` at force@0,
    /// `set_thermo` at end_of_step@+20, etc.).
    template <typename C>
    void set_singular_at_(C&& comp,
                          std::optional<std::string>& slot,
                          const char* setter_name,
                          PhaseEnum phase, int order) {
        const std::string n = comp.name();
        if (slot) {
            if (*slot == n) return;
            throw std::invalid_argument(
                std::string(setter_name) + " already configured with '" +
                *slot + "', rejected re-call with '" + n + "'");
        }
        slot = n;
        this->on(phase).add(std::forward<C>(comp), order);
    }

    /// @brief Aggregate of pending-component facts that `install_infra()`
    /// uses to decide what to auto-inject and what to throw about. Single
    /// pass over `pending_` populates everything.
    struct InfraScan {
        bool needs_stress      = false;  ///< Some pending component reads Stress.
        bool has_virial_stress = false;  ///< Canonical Stress provider already there.
        bool has_neighbor_list = false;  ///< NeighborList already there.
        std::string pair_evaluator_name; ///< First component reporting `needs_neighbor_list()`.
        std::optional<Real> pair_evaluator_cutoff;  ///< Set when the above is a `PairEvaluator`.
    };

    /// @brief Walk the pending list once, collect everything install_infra
    /// rules might need.
    InfraScan scan_pending_infra_() const {
        InfraScan s;
        for_each_pending([&](const cpu::ComponentBase& c) {
            if (has_capability(c.reads(), Capability::Stress))
                s.needs_stress = true;
            const std::string n = c.name();
            if (n == "virial_stress") s.has_virial_stress = true;
            if (n == "neighbor_list") s.has_neighbor_list = true;
            if (c.needs_neighbor_list() && s.pair_evaluator_name.empty()) {
                s.pair_evaluator_name = n;
                if (auto* pe = dynamic_cast<const cpu::PairEvaluator*>(&c))
                    s.pair_evaluator_cutoff = pe->pair_cutoff();
            }
        });
        return s;
    }

    /// @brief Shared missing-NL check. `setter_help` is the
    /// driver-specific instruction printed in the throw message
    /// (`set_pair` is the same, but the variable name differs:
    /// `sim` / `opt` / `eval`).
    void throw_if_missing_neighbor_list_(const InfraScan& s,
                                         const char* driver_name,
                                         const char* setter_help) const {
        if (!s.pair_evaluator_name.empty() && !s.has_neighbor_list) {
            throw std::runtime_error(
                std::string(driver_name) + "::build(): component '" +
                s.pair_evaluator_name + "' reports needs_neighbor_list()==true "
                "but no neighbor list was registered. " + setter_help);
        }
    }
};

} // namespace cpu
} // namespace ATV
