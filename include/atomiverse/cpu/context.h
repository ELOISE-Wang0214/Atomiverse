/// @file cpu/context.h
/// @brief cpu::Context -- CPU execution context with OpenMP, MPI, and
/// (MD-only) spatial-decomposition awareness.
#pragma once

#include "atomiverse/facade/iteration_status.h"

namespace ATV {
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
