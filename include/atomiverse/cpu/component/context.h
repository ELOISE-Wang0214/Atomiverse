/// @file cpu/component/context.h
/// @brief cpu::Context -- CPU execution context with OpenMP, MPI, and
/// (MD-only) spatial-decomposition awareness.
#pragma once

#include <functional>

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

/// @brief CPU execution context holding thread, MPI, and MD-DD configuration.
///
/// Thread count: 0 = auto-detect via omp_get_max_threads().
/// MPI fields: rank/size are populated in Backend::init() when MPI is available.
/// `dd`: non-owning pointer to the MD spatial decomposition owned by
/// `cpu::Backend`. Non-null after `Backend::init()`. Non-MD drivers simply
/// ignore this field. The full `DomainDecomposition` type is deliberately
/// not included here to keep `cpu::Context` free of MD-specific includes
/// (SCF / minimizer do not need MPI headers transitively).
/// `pair_list`: non-owning pointer to the active MD neighbor list,
/// published by `NeighborList` during `build()` and consumed by downstream
/// force components (LJ, etc.). Kept here -- rather than on `cpu::Frame`
/// -- so that the frame stays confined to fundamental physical quantities
/// and non-MD drivers (minimizer, SCF) that never see a neighbor list
/// simply ignore the field.
/// `executor`: non-owning pointer to the analysis executor (if any) used
/// by async observers/recorders. May be null on drivers that do not host
/// an executor.
struct Context {
    int n_threads = 0; ///< OpenMP threads (0 = auto-detect).
    int mpi_rank = 0;  ///< MPI rank within communicator (0 if no MPI).
    int mpi_size = 1;  ///< MPI communicator size (1 if no MPI).

    /// @brief Non-owning handle to the MD spatial decomposition.
    /// Allocated by `Backend::init()`, populated by `Backend::set_frame()`.
    md::mpi::DomainDecomposition* dd = nullptr;

    /// @brief Non-owning handle to the active MD neighbor list.
    /// Published by `NeighborList::build()` (and re-published on each
    /// `launch()` as a no-op refresh -- address is stable across the run).
    /// Null on drivers / configurations that don't carry a neighbor list.
    const md::PairList* pair_list = nullptr;

    /// @brief Non-owning handle to the per-rank analysis executor.
    /// Allocated by `Backend::init()`; published so components can reach
    /// the executor from any `BuildContext` / `LaunchContext`.
    AnalysisExecutor* executor = nullptr;

    /// @brief Non-owning handle to the per-Backend timing registry.
    /// Published by `Backend::init()`. Components reach it (via
    /// `ATV_PROFILE_SCOPE`) to record sub-region timings; null / disabled by
    /// default so timing is opt-in and zero-cost when off. @see core/timer.h
    ATV::TimerRegistry* timer = nullptr;

    /// @brief True only while the final step of the current `run(n_steps)`
    /// call is executing. Set by the Dynamics driver loop before the last
    /// `step_once`, cleared otherwise. Lets cadence-gated end-of-step
    /// collectives (e.g. `mpi::Reduce` under `set_scalar_reduce_every > 1`)
    /// guarantee the globally-reduced scalars are valid when `run()` returns —
    /// so a large reduce cadence stays safe for the common "read scalars after
    /// the run" pattern without forcing a per-step collective.
    bool is_last_step = false;

    /// @brief Re-entrant force re-evaluation handle (Optimize drivers only).
    ///
    /// Published by `cpu::OptimizeDriver::run` after `setup_pass`; null on
    /// every other driver. Invoking it re-runs the driver's force-phase
    /// bindings (ZeroForces/ZeroEnergy → evaluator → MPI reduce) on the
    /// frame's *current* geometry, updating `frame.scalars.energy` and
    /// `frame.atoms.f*`. An Optimizer component calls it during the
    /// `line_search` phase to evaluate trial geometries `x + alpha p` without
    /// owning a loop — the Driver still owns the outer loop. `need_grad` is
    /// reserved for a future energy-only fast path; the current force
    /// evaluators always produce energy and forces together. Always check for
    /// non-null before calling (only set under OptimizeDriver).
    std::function<void(bool need_grad)> reeval_forces = nullptr;

    /// @brief True while an Optimizer is evaluating line-search trial
    /// geometries via `reeval_forces`. Migration-type components
    /// (`comm::Migrate`) skip when set: the optimizer's per-row buffers are
    /// keyed to the current owned-row order, so reshuffling rows mid-search
    /// would corrupt it. Halo exchange still runs, keeping trial forces
    /// correct for the bounded step. Cleared before the optimizer commits the
    /// accepted geometry, so the driver's normal force phase migrates as usual.
    bool line_search_trial = false;

    /// @brief Outer-driver-loop status — owned by the category driver
    ///        (Evaluate / Dynamics / Optimize), reset before each iteration,
    ///        consulted at the end of each iteration to decide early stop.
    ///
    /// Only `ConvergenceCheck`-tagged components write to this. SCF
    /// internals manage their own internal convergence and never touch
    /// the outer status. See facade/iteration_status.h for the contract.
    IterationStatus outer_status{};
};

} // namespace cpu
} // namespace ATV
