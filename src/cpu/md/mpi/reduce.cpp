/// @file src/cpu/md/mpi/reduce.cpp
/// @brief Phase 1 no-op implementation of `ReduceEnergy`.

#include "atomiverse/cpu/md/mpi/reduce.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

void ReduceEnergy::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size == 1) return;
    // Phase 3: batched MPI_Allreduce of {energy, kinetic, virial[6], n_local}.
}

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
