/// @file cpu/component/context.h
/// @brief cpu::Context -- CPU execution context: flat core state plus
/// per-workflow service slices (arch-refactor-05-context).
///
/// The core execution state (threads, MPI identity, timer, outer status)
/// stays flat — it is read by every component. Workflow-specific services
/// are grouped into nested plain-aggregate slices, each filled only by its
/// owning category driver: `md` (classical force pipeline), `mc`
/// (Metropolis sampling), `opt` (optimizer line search), `analysis`
/// (async analysis executor). Unfilled slice fields stay `nullptr` =
/// service unavailable; components check before use. Cross-slice reads
/// are legitimate (slices group by *publisher*, not by exclusive
/// consumer): `comm::Migrate` reads `opt.line_search_trial`, Minimize's
/// force phase reads `md.pair_list`.
///
/// The grouping is zero-cost: one aggregate, stack-constructible,
/// constant-offset member access — no virtual dispatch, no registry, no
/// heap (hot-path invariant of the arch-refactor chain).
#pragma once

#include <functional>
#include <random>
#include <type_traits>

#include "atomiverse/facade/iteration_status.h"

namespace ATV {

class TimerRegistry;

namespace cpu {

class AnalysisExecutor;

namespace md {
struct PairList;
namespace mpi {
struct DomainDecomposition;
} // namespace mpi
} // namespace md

namespace mc {
struct Transaction;
} // namespace mc

/// @brief Classical force-pipeline services (publishers: `Backend`,
///        `NeighborList`, `DynamicsDriver`).
struct DynamicsServices {
    /// @brief Non-owning handle to the MD spatial decomposition.
    /// Allocated by `Backend::init()`, populated by `Backend::set_frame()`.
    /// The full type is deliberately not included here so `cpu::Context`
    /// stays free of MD-specific includes (SCF / minimizer do not need
    /// MPI headers transitively).
    md::mpi::DomainDecomposition* dd = nullptr;

    /// @brief Non-owning handle to the active neighbor list.
    /// Published by `NeighborList::build()` (address stable across the
    /// run). Kept here — rather than on `cpu::Frame` — so the frame stays
    /// confined to fundamental physical quantities.
    const md::PairList* pair_list = nullptr;

    /// @brief Opt-in: publish a FULL neighbor list instead of the half list.
    ///
    /// A consumer that cannot work from a half list sets this in its `build()`;
    /// `NeighborList` then calls `CellList::refresh_full` and the published
    /// `pair_list` carries, for every owned `i`, *all* of its neighbours in both
    /// directions. Mirrors `cuda::DynamicsServices::flat_pairs_requested`.
    ///
    /// Who needs it: any style whose per-atom term is not a sum over *unordered*
    /// pairs — the EAM embedding density `rho_i = sum_j phi(r_ij)` (which must see
    /// every neighbour of `i`, not the half of them that happen to sit in `i`'s
    /// row), and every 3-body style (Stillinger-Weber, Tersoff), whose angular
    /// term enumerates `(j,k)` pairs drawn from `i`'s complete neighbour set.
    ///
    /// Read in `NeighborList::launch()`, not `build()`: component `build()` runs
    /// in Plan order, and the neighbor phase precedes the force phase, so the
    /// requester's `build()` runs *after* the provider's. The provider therefore
    /// treats a mode change as a rebuild trigger — same lazy resolution as the
    /// CUDA flat-pairs export.
    bool full_list_requested = false;

    /// @brief True only while the final step of the current `run(n_steps)`
    /// call is executing. Set by the Dynamics driver loop before the last
    /// `step_once`, cleared otherwise. Lets cadence-gated end-of-step
    /// collectives (e.g. `mpi::Reduce` under `set_scalar_reduce_every > 1`)
    /// guarantee globally-reduced scalars are valid when `run()` returns.
    bool is_last_step = false;
};

/// @brief Metropolis-sampling services (publisher: `SamplingDriver::run`).
struct SamplingServices {
    /// @brief Non-owning handle to the active MC trial-move changeset.
    /// Address-stable for the run; a `Mover` records its proposed
    /// displacement here each `propose` phase, and the driver reads
    /// `old_*` to roll back a rejected move. Null on non-sampling drivers
    /// (mirrors `md.pair_list` — kept off `cpu::Frame` for Frame purity).
    mc::Transaction* proposal = nullptr;

    /// @brief Non-owning handle to the sampling driver's single seeded RNG
    /// stream. A `Mover` draws atom picks and displacements from it so the
    /// whole trajectory is reproducible from the driver's one seed.
    std::mt19937_64* rng = nullptr;
};

/// @brief Optimizer line-search services (publishers: `OptimizeDriver::run`
///        and the active Optimizer component).
struct OptimizeServices {
    /// @brief Re-entrant force re-evaluation handle.
    ///
    /// Published by `cpu::OptimizeDriver::run` after `setup_pass`; null on
    /// every other driver. Invoking it re-runs the driver's force-phase
    /// bindings (ZeroForces/ZeroEnergy → evaluator → MPI reduce) on the
    /// frame's *current* geometry. An Optimizer component calls it during
    /// the `line_search` phase to evaluate trial geometries `x + alpha p`
    /// without owning a loop — the Driver still owns the outer loop.
    /// `need_grad` is reserved for a future energy-only fast path. Always
    /// check for non-null before calling.
    std::function<void(bool need_grad)> reeval_forces = nullptr;

    /// @brief True while an Optimizer is evaluating line-search trial
    /// geometries via `reeval_forces`. Migration-type components
    /// (`comm::Migrate`) skip when set: the optimizer's per-row buffers are
    /// keyed to the current owned-row order, so reshuffling rows mid-search
    /// would corrupt it. Halo exchange still runs, keeping trial forces
    /// correct for the bounded step.
    bool line_search_trial = false;
};

/// @brief Async-analysis services (publisher: `AnalysisExecutor::install`).
struct AnalysisServices {
    /// @brief Non-owning handle to the per-rank analysis executor.
    /// Allocated by `Backend::init()`; published so async observers /
    /// recorders can reach the executor from any lifecycle context.
    AnalysisExecutor* executor = nullptr;
};

/// @brief CPU execution context: flat core state + per-workflow slices.
struct Context {
    int n_threads = 0; ///< OpenMP threads (0 = auto-detect).
    int mpi_rank = 0;  ///< MPI rank within communicator (0 if no MPI).
    int mpi_size = 1;  ///< MPI communicator size (1 if no MPI).

    /// @brief Non-owning handle to the per-Backend timing registry.
    /// Published by `Backend::init()`. Components reach it (via
    /// `ATV_PROFILE_SCOPE`) to record sub-region timings; null / disabled by
    /// default so timing is opt-in and zero-cost when off. @see core/timer.h
    ATV::TimerRegistry* timer = nullptr;

    /// @brief Outer-driver-loop status — owned by the category driver,
    ///        reset before each iteration, consulted at the end of each
    ///        iteration to decide early stop.
    ///
    /// Only `ConvergenceCheck`-tagged components write to this. SCF
    /// internals manage their own internal convergence and never touch
    /// the outer status. See facade/iteration_status.h for the contract.
    IterationStatus outer_status{};

    DynamicsServices md;       ///< Classical force-pipeline slice.
    SamplingServices mc;       ///< Metropolis-sampling slice.
    OptimizeServices opt;      ///< Optimizer line-search slice.
    AnalysisServices analysis; ///< Async-analysis slice.
};

// Zero-cost-grouping guards: the slices stay plain aggregates and the
// Context stays non-polymorphic (constant-offset access, stack
// construction). OptimizeServices is exempt from trivial copyability —
// its std::function predates this grouping.
static_assert(std::is_trivially_copyable_v<DynamicsServices>);
static_assert(std::is_trivially_copyable_v<SamplingServices>);
static_assert(std::is_trivially_copyable_v<AnalysisServices>);
static_assert(!std::is_polymorphic_v<Context>);

} // namespace cpu
} // namespace ATV
