/// @file src/cpu/md/mpi/atom_migration.cpp
/// @brief Phase 1 no-op implementation of `AtomMigration`.

#include "atomiverse/cpu/md/mpi/atom_migration.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

void AtomMigration::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size == 1) return;
    // Phase 3: MPI_Alltoallv of atoms + bonded tuples to new owners;
    // bump dd->neighbor_generation on successful migration.
}

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
