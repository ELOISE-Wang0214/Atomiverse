/// @file cuda/component.cuh
/// @brief ComponentBase -- OOP base class for all simulation components.
///
/// Components follow a strict 3-phase lifecycle:
///   build()   -- one-time resource preparation (after Runtime + DeviceFrame exist)
///   launch()  -- per-step GPU work submission (pure GPU, zero D2H/sync/malloc)
///   destroy() -- resource cleanup
///
/// Two registration paths:
/// - `sim.add(comp)` -- uses default_bindings().
/// - `sim.on(phase).add(comp)` -- explicit phase override.
///
/// CUDA-free header.
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/core/types.h"
#include "atomiverse/cuda/context.cuh"
#include "atomiverse/cuda/frame.cuh"
#include "atomiverse/facade/component_role.h"

namespace ATV {
namespace cuda {

/// @brief Context for one-time component build (after Runtime + DeviceFrame exist).
///
/// Passed to ComponentBase::build(). Components may perform cudaMalloc,
/// cudaMemcpy, cudaStreamSynchronize during build -- any host/device work.
///
/// @see ComponentBase::build()
struct BuildContext {
    cuda::Context& runtime; ///< CUDA execution context (streams).
    cuda::Frame& frame;     ///< Device-resident system state + host metadata.
};

/// @brief Context for per-step GPU work submission.
///
/// Passed to ComponentBase::launch() every simulation step. launch() must
/// remain GPU-pure: only kernel launches, cudaMemsetAsync, cudaMemcpyAsync D2D.
///
/// @see ComponentBase::launch()
struct LaunchContext {
    cuda::Context& runtime; ///< CUDA execution context (streams).
    cuda::Frame& frame;     ///< Device-resident system state + host metadata.
    int stage;              ///< Current phase index within the step.
    int step;               ///< Current simulation step index (0-based).
};

/// @brief Context for component resource cleanup.
///
/// Passed to ComponentBase::destroy() after the simulation loop completes.
/// Components must release all GPU resources allocated in build().
///
/// @see ComponentBase::destroy()
struct DestroyContext {
    cuda::Context& runtime; ///< CUDA execution context (streams).
    cuda::Frame& frame;     ///< Device-resident system state + host metadata.
};

/// @brief Abstract base for simulation components.
///
/// A Component is a CPU-side facade for staged algorithmic behavior.
/// It owns configuration, persistent state, and lifecycle, but not the
/// global loop. Components follow a strict **build / launch / destroy**
/// lifecycle managed by the Backend.
///
/// Two registration paths:
/// @code{.cpp}
/// sim.add(comp);                  // uses default_bindings()
/// sim.on(phase).add(comp, order); // explicit phase override
/// @endcode
///
/// @see Driver for the user-facing registration API
/// @see Backend for lifecycle dispatch
class ComponentBase {
  public:
    virtual ~ComponentBase() = default;

    /// @brief Human-readable component name (used in timing reports).
    /// @return Component identifier string.
    virtual std::string name() const = 0;

    /// @brief Structural role within an outer driver loop.
    ///
    /// Mirrors `cpu::ComponentBase::role()` (Phase 5). Defaults to `None`
    /// for structural infrastructure components; concrete components
    /// override to declare their slot (Integrator / Evaluator / ...).
    /// @see facade/component_role.h
    virtual ComponentRole role() const { return ComponentRole::None; }

    /// @brief Bitset of frame fields this component may write.
    ///
    /// Mirrors `cpu::ComponentBase::capabilities()`. Defaults to `None`;
    /// concrete evaluators override to declare what they produce.
    /// @see facade/component_role.h
    virtual Capability capabilities() const { return Capability::None; }

    /// @brief Bitset of frame fields this component consumes.
    ///
    /// Mirrors `cpu::ComponentBase::reads()`. Defaults to `None`;
    /// concrete components override to declare their inputs. Used by
    /// `cuda::md::MD::install_infra_` for capability-typed auto-injection
    /// — e.g. a barostat that declares `reads()==Stress` auto-pulls in
    /// `cuda::VirialStress` if none is registered.
    /// @return Capability bitset of frame fields read (defaults to None).
    /// @see facade/component_role.h
    virtual Capability reads() const { return Capability::None; }

    /// @brief Whether this component requires a `cuda::NeighborList` in
    ///        the same plan.
    ///
    /// Defaults to `false`. Pair evaluators (`cuda::Pair` subclasses)
    /// override to return `true`; `cuda::md::MD::install_infra_` uses
    /// this to detect missing neighbor lists and throws with the
    /// offender's name.
    /// @return true if this component needs a `cuda::NeighborList`
    ///         registered in the same plan.
    virtual bool needs_neighbor_list() const { return false; }

    /// @brief Whether this component *is* a neighbor list (the supply
    ///        side of the `needs_neighbor_list()` contract).
    ///
    /// Defaults to `false`. `cuda::NeighborList` overrides to return
    /// `true`. Used by `cuda::md::MD::install_infra_` to detect a
    /// registered NeighborList without keying on `name() == "neighbor_list"`.
    /// @return true if this component supplies the neighbor list to
    ///         downstream consumers.
    virtual bool provides_neighbor_list() const { return false; }

    /// @brief Scheduling trait (blocking or deferrable).
    /// @return Scheduling::blocking (default) or Scheduling::deferrable.
    virtual Scheduling scheduling() const { return Scheduling::blocking; }

    /// @brief Execution frequency: 0 = every step, N = every Nth step.
    /// @return Step frequency. 0 means execute on every step.
    virtual int frequency() const { return 0; }

    /// @brief Default stage bindings. Empty = must bind explicitly via on().add().
    /// @return Vector of DefaultBinding (phase index + order pairs).
    virtual std::vector<DefaultBinding> default_bindings() const { return {}; }

    /// @brief One-time resource preparation.
    ///
    /// Called once after Runtime and DeviceFrame are materialized.
    /// May perform cudaMalloc, D2H reads, sync -- any host/device work.
    /// All persistent GPU resources for repeated launch() must be created here.
    ///
    /// @param ctx Build context with Runtime and DeviceFrame references.
    virtual void build(const BuildContext& ctx) {}

    /// @brief Per-step GPU work submission.
    ///
    /// Called every step during the simulation loop.
    /// Must remain GPU-pure: only kernel launches, cudaMemsetAsync,
    /// cudaMemcpyAsync D2D. Zero cudaMalloc, zero D2H, zero sync.
    ///
    /// @param ctx Launch context with Runtime, DeviceFrame, phase, and step.
    virtual void launch(const LaunchContext& ctx) = 0;

    /// @brief Resource cleanup.
    ///
    /// Called once after the simulation loop completes.
    /// Must release all resources allocated in build().
    /// Executed in reverse Plan order.
    ///
    /// @param ctx Destroy context with Runtime and DeviceFrame references.
    virtual void destroy(const DestroyContext& ctx) {}
};

} // namespace cuda
} // namespace ATV
