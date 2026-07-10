/// @file cpu/component/component.h
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
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cpu {

class HaloParticipant; // opt-in halo protocol (halo_participant.h)

/// @brief Context for one-time component build on CPU.
/// @see cpu::ComponentBase::build()
struct BuildContext {
    cpu::Context& runtime; ///< CPU execution context.
    cpu::Frame& frame;     ///< Host-resident system state.
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
    /// == Stress` triggers auto-injection of `Virial`). Defaults
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

    /// @brief Whether this component *is* a neighbor structure (the supply
    /// side of the `needs_neighbor_list()` contract).
    ///
    /// Defaults to `false`. `cpu::NeighborList` overrides to `true`. Used
    /// by `CPUDriver::scan_pending_infra_()` to detect a registered
    /// neighbor structure without keying on `name() == "neighbor_list"` —
    /// a user-supplied neighbor structure participates in auto-injection
    /// by overriding this, whatever its `name()`. Mirrors
    /// `cuda::ComponentBase::provides_neighbor_list()`.
    virtual bool provides_neighbor_list() const { return false; }

    /// @brief Typed accessor for the opt-in halo protocol.
    ///
    /// A component that participates in the MPI ghost exchange (the
    /// `comm::HaloSend` / `comm::HaloRecv` endpoints, or a field
    /// participant declaring extra per-atom payload) inherits
    /// `cpu::HaloParticipant` alongside `ComponentBase` and overrides this
    /// to return `this`. The driver infrastructure scan and the build-time
    /// `bind_participants` pass probe this single hook — no `dynamic_cast`,
    /// no `name()` matching. Default: not a participant.
    /// @return The participant view of this component, or `nullptr`.
    virtual HaloParticipant* halo_participant() { return nullptr; }

    /// @brief Neighbor-interaction cutoff hint [Å] for driver auto-injection.
    ///
    /// Components that report `needs_neighbor_list()` may override this to
    /// expose their interaction cutoff; the category drivers use it to size
    /// the auto-injected MPI halo exchange (skin-padded). `cpu::PairEvaluator`
    /// forwards its `pair_cutoff()`. Return `0` (the default) when the
    /// component declares no cutoff — auto-injection then skips the halo.
    /// Replaces the former `dynamic_cast<PairEvaluator*>` probe.
    virtual Real neighbor_cutoff_hint() const { return Real(0); }

    /// @brief Scheduling trait (blocking or deferrable).
    /// @return `Scheduling::blocking` for CPU components unless overridden.
    virtual Scheduling scheduling() const { return Scheduling::blocking; }

    /// @brief Execution frequency: 0 = every step, N = every Nth step.
    /// @return Step stride. `0` means launch on every step.
    virtual int frequency() const { return 0; }

    /// @brief Heap bytes this component owns (workspaces, neighbor lists,
    /// per-thread buffers, comm buffers). Used by the built-in end-of-run
    /// driver summary to attribute peak memory to a category (the component's
    /// `name()` classifies it, mirroring the timing breakdown). Default 0 for
    /// stateless components (integrators, Zero*, …). Override to report owned
    /// allocations; sum the `capacity()` of every persistent buffer so the
    /// figure tracks the resident footprint. Frame / per-atom storage is
    /// accounted separately by the driver, not here.
    virtual std::size_t memory_bytes() const { return 0; }

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

    // The former in-vtable halo contract (comm_forward/reverse +
    // pack/unpack quartet, LAMMPS src/pair.h pattern) moved to the opt-in
    // `cpu::HaloParticipant` protocol (halo_participant.h) — declarative
    // CommField descriptors with engine-side packing, shared shape with
    // the CUDA backend. See arch-refactor-02-halo.
};

} // namespace cpu
} // namespace ATV
