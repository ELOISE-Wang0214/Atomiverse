/// @file cpu/context.h
/// @brief cpu::Context -- CPU execution context with OpenMP and MPI awareness.
#pragma once

namespace ATV {
namespace cpu {

/// @brief CPU execution context holding thread and MPI configuration.
///
/// Thread count: 0 = auto-detect via omp_get_max_threads().
/// MPI fields: rank/size are populated in Backend::init() when MPI is available.
struct Context {
    int n_threads = 0; ///< OpenMP threads (0 = auto-detect).
    int mpi_rank = 0;  ///< MPI rank within communicator (0 if no MPI).
    int mpi_size = 1;  ///< MPI communicator size (1 if no MPI).
};

} // namespace cpu
} // namespace ATV
