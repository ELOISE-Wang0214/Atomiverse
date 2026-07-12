/// @file driver.h
/// @brief Driver -- user-facing simulation facade with fluent API.
///
/// Backend-agnostic template: DeviceTag selects CPU or CUDA at compile time.
/// All method bodies are defined inline (template class). Explicit
/// instantiations live in src/cpu/driver/driver.cpp and src/cuda/driver/driver.cu.
#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "atomiverse/core/backend_traits.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/core/device.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/device/timer.h"
#include "atomiverse/facade/plan.h"
#include "atomiverse/facade/stage.h"

namespace ATV {

/// @brief GPU profiler type held by `Driver` (CUDA-free PIMPL facade).
///
/// `cuda::Timer` has a complete CUDA-free definition (`cuda/timer.h`)
/// with all special members declared and defined out-of-line in
/// `src/cuda/device/timer.cu`. This lets every translation unit — CUDA or not —
/// hold and destroy a `std::unique_ptr<Timer>` member. On non-CUDA builds
/// the profiler is simply never constructed (`enable_timing()` is a
/// no-op for `cpu::device_tag`).
using Timer = cuda::Timer;

namespace detail {
/// @brief Unresolved component binding collected before `Driver::build()`.
///
/// The driver stores these lightweight records while the user is still calling
/// `add()` and `on(...).add(...)`. `PlanBuilder` later resolves them into
/// sorted `StageSlot` operations.
template <typename Component>
struct PendingBinding {
    int phase;            ///< Target phase index.
    int order;            ///< Ordering within the phase (lower runs earlier).
    Component* component; ///< Non-owning pointer to the registered component.
};
} // namespace detail

/// @brief Fluent phase-binding helper returned by Driver::on().
template <typename PhaseEnum, typename DeviceTag>
class PhaseBuilder {
    using Component = typename backend_traits<DeviceTag>::component_type;

  public:
    /// @brief Register a component explicitly on the current phase.
    ///
    /// Ownership follows the value category of the argument (symmetric with
    /// `Driver::add`): rvalue → moved into the driver; lvalue → stored as a
    /// non-owning pointer, caller retains lifetime.
    /// @tparam Comp Concrete component type.
    /// @param comp Component instance.
    /// @param order Ordering within the phase (lower runs earlier).
    /// @return `*this` for fluent chaining.
    template <typename Comp>
    PhaseBuilder& add(Comp&& comp, int order = 0) {
        using C = std::remove_reference_t<Comp>;
        if constexpr (std::is_lvalue_reference_v<Comp>) {
            bindings_->push_back({ATV::slot_of(phase_), order, &comp});
        } else {
            auto ptr = std::make_unique<C>(std::forward<Comp>(comp));
            auto* raw = ptr.get();
            owned_->push_back(std::move(ptr));
            bindings_->push_back({ATV::slot_of(phase_), order, raw});
        }
        return *this;
    }

  private:
    template <typename PE, typename DT>
    friend class Driver;
    PhaseBuilder(std::vector<detail::PendingBinding<Component>>* b,
                 std::vector<std::unique_ptr<Component>>* o, PhaseEnum p)
        : bindings_(b), owned_(o), phase_(p) {}
    std::vector<detail::PendingBinding<Component>>* bindings_;
    std::vector<std::unique_ptr<Component>>* owned_;
    PhaseEnum phase_;
};

/// @brief Simulation driver -- user-facing facade.
///
/// @tparam PhaseEnum Phase enum type (must end with COUNT sentinel).
/// @tparam DeviceTag Device tag (cpu::device_tag or cuda::device_tag).
template <typename PhaseEnum, typename DeviceTag>
class Driver {
  protected:
    using traits = backend_traits<DeviceTag>;
    using Component = typename traits::component_type;
    using BackendT = typename traits::backend_type;

  public:
    Driver() = default;
    ~Driver() { shutdown(); }
    Driver(const Driver&) = delete;
    Driver& operator=(const Driver&) = delete;

    Driver(Driver&& o) noexcept
        : pending_(std::move(o.pending_)),
          owned_components_(std::move(o.owned_components_)),
          plan_(std::move(o.plan_)), backend_(std::move(o.backend_)),
          timer_(std::move(o.timer_)) {}

    Driver& operator=(Driver&& o) noexcept {
        if (this != &o) {
            shutdown();
            pending_ = std::move(o.pending_);
            owned_components_ = std::move(o.owned_components_);
            plan_ = std::move(o.plan_);
            backend_ = std::move(o.backend_);
            timer_ = std::move(o.timer_);
        }
        return *this;
    }

    /// @brief Add a component using its `default_bindings()`.
    ///
    /// Two ownership modes selected by the value category of the argument:
    ///   - **rvalue** (e.g. `sim.add(VelocityVerlet{0.001})`): driver takes
    ///     ownership by moving the component into a `unique_ptr`.
    ///   - **lvalue** (e.g. `MyObserver obs; sim.add(obs)`): driver stores
    ///     a non-owning pointer; caller retains lifetime, allowing post-run
    ///     state introspection. This is the preferred form for components
    ///     that are non-copyable / non-movable (e.g. `std::mutex` members)
    ///     or whose state the caller wants to read after `run()` returns.
    ///
    /// @tparam Comp Concrete component type (deduced).
    /// @param comp Component instance.
    template <typename Comp>
    void add(Comp&& comp) {
        using C = std::remove_reference_t<Comp>;
        if constexpr (std::is_lvalue_reference_v<Comp>) {
            // Non-owning: caller retains lifetime.
            auto bindings = comp.default_bindings();
            if (bindings.empty())
                throw std::runtime_error(
                    "Component '" + comp.name()
                    + "' has no default bindings. Use .on(phase).add().");
            validate_binding_kinds_(comp.name(), bindings);
            for (auto& db : bindings)
                pending_.push_back({db.phase, db.order, &comp});
        } else {
            auto ptr = std::make_unique<C>(std::forward<Comp>(comp));
            auto bindings = ptr->default_bindings();
            if (bindings.empty())
                throw std::runtime_error(
                    "Component '" + ptr->name()
                    + "' has no default bindings. Use .on(phase).add().");
            validate_binding_kinds_(ptr->name(), bindings);
            auto* raw = ptr.get();
            owned_components_.push_back(std::move(ptr));
            for (auto& db : bindings)
                pending_.push_back({db.phase, db.order, raw});
        }
    }

    /// @brief Begin explicit phase binding.
    /// @param phase Phase value returned by the driver's `Phase` enum.
    /// @return A `PhaseBuilder` bound to that phase.
    PhaseBuilder<PhaseEnum, DeviceTag> on(PhaseEnum phase) {
        return PhaseBuilder<PhaseEnum, DeviceTag>(&pending_, &owned_components_,
                                                  phase);
    }

    /// @brief Read-only access to the compiled Plan (after `build()`).
    /// Useful for tests that need to compare Plans built two ways
    /// (e.g. Plan-equivalence between the categorical setter API and
    /// the raw `add()`/`on()` path).
    const Plan<Component>& plan() const { return plan_; }

    /// @brief Compile the execution Plan.
    void build() {
        PlanBuilder<Component> builder;
        for (auto& pb : pending_)
            builder.add({pb.phase, pb.order, pb.component});
        pending_.clear();
        plan_ = builder.build(ATV::slot_count<PhaseEnum>());
    }

    /// @brief Run the full simulation lifecycle.
    /// @param frame Host-side frame used as input and final output.
    /// @param n_steps Number of time steps or SCF iterations to execute.
    void run(cpu::Frame& frame, int n_steps);

    /// @brief Enable per-component timing (call before run).
    void enable_timing();

    /// @brief Timing summary (call after run).
    /// @return Human-readable per-phase timing report, or an empty string on CPU.
    std::string summary() const;

    /// @brief Shut down Backend and release all resources.
    void shutdown();

  protected:
    // Internals are protected (not private) so the backend loop node
    // (cpu::CPUDriver / cuda::CUDADriver) can drive the per-iteration loop itself.
    // User code does not subclass `Driver<>` directly — it uses a preset
    // (cpu::md::MD, cpu::minimizer::Minimize, …) or subclasses the loop node.
    std::vector<detail::PendingBinding<Component>> pending_;
    std::vector<std::unique_ptr<Component>> owned_components_;
    Plan<Component> plan_;
    std::unique_ptr<BackendT> backend_;
    std::unique_ptr<Timer> timer_;

  private:
    /// @brief Registration-time gate: reject default bindings built from
    ///        another workflow's Phase enum.
    ///
    /// The five Phase enums share overlapping numeric values (all start at
    /// 0), so an unchecked cross-workflow binding lands on a syntactically
    /// valid but semantically wrong slot (e.g. `evaluate::evaluate` == 1 ==
    /// `md::post_integrate`) and runs silently in the wrong lifecycle
    /// position. The `PhaseKind` discriminant carried by every binding (via
    /// `ATV::make_binding`) makes the mismatch detectable exactly where the
    /// user can fix it — the `add()` call.
    ///
    /// @param name Component name, for the error message.
    /// @param bindings Bindings returned by the component's
    ///        `default_bindings()`.
    /// @throws std::runtime_error naming both workflows and pointing at
    ///         `.on(phase).add()` when any binding's kind mismatches.
    void
    validate_binding_kinds_(const std::string& name,
                            const std::vector<DefaultBinding>& bindings) const {
        const PhaseKind expected = kind_of(PhaseEnum{});
        for (const auto& db : bindings) {
            if (db.kind != expected)
                throw std::runtime_error(
                    "Component '" + name
                    + "' declares default bindings for the "
                    + std::string(phase_kind_name(db.kind))
                    + " workflow, but this driver executes "
                    + std::string(phase_kind_name(expected))
                    + " phases. Bind explicitly with .on(phase).add().");
        }
    }
};

} // namespace ATV
