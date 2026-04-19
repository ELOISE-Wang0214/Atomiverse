/// @file src/cpu/md/mpi/reduce.cpp
/// @brief `ReduceEnergy` -- batched end-of-step MPI_Allreduce of the
///        per-rank scalars and stress tensor (Phase 0.4 / P3.3).
///
/// Bound at `end_of_step @ -50`, after `KineticEnergy @ -100` (so KE is
/// the locally-summed value) and before `SnapshotPublisher @ +10` (so any
/// async observer reading from the snapshot sees the globally-reduced
/// scalars). Single-rank: structural no-op.

#include "atomiverse/cpu/md/mpi/reduce.h"

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

#ifdef ATV_HAS_MPI

void ReduceEnergy::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1) return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL) return;

    auto& f = ctx.frame;
    // Pack {energy, kinetic, stress[6]} into a single buffer for one
    // 8-double allreduce. n_local is omitted -- it is authoritative on
    // each rank from `dd->n_local` and does not need to be summed.
    double sendbuf[8];
    sendbuf[0] = f.energy;
    sendbuf[1] = f.kinetic;
    for (int k = 0; k < 6; ++k) sendbuf[2 + k] = f.stress[k];

    double recvbuf[8];
    MPI_Allreduce(sendbuf, recvbuf, 8, MPI_DOUBLE, MPI_SUM, dd->cart_comm);

    f.energy  = recvbuf[0];
    f.kinetic = recvbuf[1];
    for (int k = 0; k < 6; ++k) f.stress[k] = recvbuf[2 + k];
}

#else // !ATV_HAS_MPI

void ReduceEnergy::launch(const cpu::LaunchContext& /*ctx*/) {}

#endif // ATV_HAS_MPI

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
