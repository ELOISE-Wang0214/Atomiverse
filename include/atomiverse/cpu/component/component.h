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
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {

class HaloParticipant;  // opt-in halo protocol (halo_participant.h)
class ConstraintProvider; // opt-in DOF protocol (constraint_protocol.h)
class NeighborConsumer; // opt-in neighbor protocols (neighbor_protocol.h)
class NeighborProvider;

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
    int stage;             ///< Current Plan slot index (erased — compare via `at()`).
    int step;              ///< Current simulation step index (0-based).

    /// @brief Is the component running in `phase` right now?
    ///
    /// The typed way to ask, and the only way a component should: a component bound in
    /// several phases dispatches on this. The raw `stage` int is the erasure the Plan
    /// needs (one dense array, six different Phase enums), and it is *not* comparable
    /// across workflows — all six enums start at 0, so `stage == 1` is simultaneously
    /// `md::pre_force`, `minimizer::line_search` and `evaluate::evaluate`. Going through
    /// `slot_of` keeps the type on the comparison, where the reader can see it.
    ///
    /// @tparam PhaseEnum The workflow Phase enum the component is bound in.
    /// @param phase The phase to test against.
    /// @return True when the current stage is that phase.
    template <typename PhaseEnum>
    bool at(PhaseEnum phase) const {
        return stage == ATV::slot_of(phase);
    }
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

    /// @brief Structural role within a driver loop.
    ///
    /// Drives each preset's `validate_components()` at `build()` time. Defaults to
    /// `None` for "structural infrastructure" components (Zero*, KineticEnergy, ...).
    /// Concrete components override to declare their loop position
    /// (Integrator / Optimizer / Evaluator / Thermostat / ...).
    /// @see facade/component_role.h
    virtual ComponentRole role() const { return ComponentRole::None; }

    /// @brief Bitset of frame fields this component may write.
    ///
    /// Used by each preset's `validate_components()` to check that, e.g., `cpu::md::MD`
    /// carries at least one component with `Capability::Force`. Defaults to `None`.
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

    /// @brief Typed accessor for the consumer side of the opt-in neighbor
    /// protocol (neighbor_protocol.h).
    ///
    /// A component that requires a neighbor list upstream in the same
    /// step (LAMMPS analogue: `pair_style lj/cut` needing a `neighbor`
    /// command) inherits `cpu::NeighborConsumer` and overrides this to
    /// return `this` — `cpu::PairEvaluator` does so for every pair style.
    /// The category drivers throw at `build()` when a consumer is
    /// registered without a neighbor structure, and size the MPI halo
    /// from the consumer's cutoff — no `dynamic_cast`, no `name()`
    /// matching. Default: not a consumer.
    /// @return The consumer view of this component, or `nullptr`.
    virtual const NeighborConsumer* as_neighbor_consumer() const {
        return nullptr;
    }

    /// @brief Typed accessor for the supply side of the opt-in neighbor
    /// protocol. The neighbor structure (`cpu::md::NeighborList`, or a
    /// user-supplied equivalent, whatever its `name()`) overrides this to
    /// return `this`. Default: not a provider.
    /// @return The provider view of this component, or `nullptr`.
    virtual const NeighborProvider* as_neighbor_provider() const {
        return nullptr;
    }

    /// @brief Opt-in `ConstraintProvider` protocol (constraint_protocol.h) — a
    /// component that removes degrees of freedom returns `this`. The driver probes
    /// this to compute the global DOF count before any thermostat builds.
    virtual const ConstraintProvider* as_constraint() const { return nullptr; }

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
