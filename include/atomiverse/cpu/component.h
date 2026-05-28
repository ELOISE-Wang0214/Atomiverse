/// @file cpu/component.h
/// @brief cpu::ComponentBase -- OOP base class for CPU simulation components.
///
/// Mirrors the CUDA ComponentBase with identical lifecycle semantics
/// but operates on cpu::Frame (xtensor) and cpu::Context directly.
/// No CUDA dependency.
///
/// @see cuda/component.cuh for the GPU counterpart
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/types.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

/// @brief Context for one-time component build on CPU.
/// @see cpu::ComponentBase::build()
struct BuildContext {
    cpu::Context& runtime; ///< CPU execution context.
    cpu::Frame& frame;     ///< Host-resident system state.
    void* pw_workspace =
        nullptr; ///< PWWorkspace* for plane-wave components (optional).
};

/// @brief Context for per-step computation on CPU.
/// @see cpu::ComponentBase::launch()
struct LaunchContext {
    cpu::Context& runtime; ///< CPU execution context.
    cpu::Frame& frame;     ///< Host-resident system state.
    int stage;             ///< Current phase index within the step.
    int step;              ///< Current simulation step index (0-based).
};

/// @brief Context for component resource cleanup on CPU.
/// @see cpu::ComponentBase::destroy()
struct DestroyContext {
    cpu::Context& runtime; ///< CPU execution context.
    cpu::Frame& frame;     ///< Host-resident system state.
};

/// @brief Abstract base for CPU simulation components.
///
/// Lifecycle: build() -> launch() x N -> destroy().
/// All data is host-resident (xtensor). OpenMP parallelism in launch().
class ComponentBase {
  public:
    virtual ~ComponentBase() = default;

    /// @brief Human-readable component name.
    /// @return Stable identifier used in summaries and diagnostics.
    virtual std::string name() const = 0;

    /// @brief Structural role within an outer driver loop.
    ///
    /// Drives validation in `cpu::EvaluateDriver` / `DynamicsDriver` /
    /// `OptimizeDriver` at `build()` time. Defaults to `None` for
    /// "structural infrastructure" components (Zero*, KineticEnergy, ...).
    /// Concrete components override to declare their loop position
    /// (Integrator / Optimizer / Evaluator / Thermostat / ...).
    /// @see facade/component_role.h
    virtual ComponentRole role() const { return ComponentRole::None; }

    /// @brief Bitset of frame fields this component may write.
    ///
    /// Used by the category drivers to validate that, e.g., a
    /// `DynamicsDriver` carries at least one evaluator with
    /// `Capability::Force`. Defaults to `None`.
    virtual Capability capabilities() const { return Capability::None; }

    /// @brief Bitset of frame fields this component READS each launch.
    ///
    /// Mirrors `capabilities()` (which is "what I write") on the
    /// consumer side. The CPU category drivers use this in
    /// `install_defaults()` to auto-provide canonical providers for
    /// otherwise-unmet read requirements (e.g. `NoseHooverNPT::reads()
    /// == Stress` triggers auto-injection of `VirialStress`). Defaults
    /// to `None` — most components read implicit frame fields
    /// (positions / forces) that no provider needs to be matched
    /// against.
    /// @see facade/component_role.h
    virtual Capability reads() const { return Capability::None; }

    /// @brief Whether this component requires a `NeighborList` upstream
    /// in the same step.
    ///
    /// Pair-cutoff evaluators (LJ, future Coulomb-real-space) override
    /// to `true`. The category drivers throw at `build()` time if any
    /// registered component reports `true` and no `NeighborList` was
    /// set via `set_neighbor()` / `add()`. LAMMPS analogue: needing
    /// `pair_style lj/cut` to be paired with a `neighbor` command.
    /// Defaults to `false`.
    virtual bool needs_neighbor_list() const { return false; }

    /// @brief Scheduling trait (blocking or deferrable).
    /// @return `Scheduling::blocking` for CPU components unless overridden.
    virtual Scheduling scheduling() const { return Scheduling::blocking; }

    /// @brief Execution frequency: 0 = every step, N = every Nth step.
    /// @return Step stride. `0` means launch on every step.
    virtual int frequency() const { return 0; }

    /// @brief Default stage bindings. Empty = must bind explicitly via on().add().
    /// @return Zero or more `(phase, order)` bindings for `Driver::add()`.
    virtual std::vector<DefaultBinding> default_bindings() const { return {}; }

    /// @brief One-time resource preparation.
    ///
    /// Allocate xtensor workspaces, per-thread buffers, cell lists.
    /// No allocation in launch().
    ///
    /// @param ctx Build context containing the CPU runtime and frame.
    virtual void build(const BuildContext& ctx) {}

    /// @brief Per-step computation.
    ///
    /// OpenMP parallelism. Use pre-allocated workspaces only.
    ///
    /// @param ctx Launch context containing runtime, frame, phase, and step.
    virtual void launch(const LaunchContext& ctx) = 0;

    /// @brief Resource cleanup.
    ///
    /// Release workspaces (xtensor RAII handles most).
    ///
    /// @param ctx Destroy context containing the CPU runtime and frame.
    virtual void destroy(const DestroyContext& ctx) {}
};

} // namespace cpu
} // namespace ATV
