/// @file src/cpu/md/mpi/reduce.cpp
/// @brief `Reduce` -- batched end-of-step MPI_Allreduce of the
///        per-rank scalars and stress tensor.
///
/// Bound at `end_of_step @ -50`, after `KineticEnergy @ -100` (so KE is
/// the locally-summed value) and before `Publish @ +10` (so any
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

void Reduce::launch(const cpu::LaunchContext& ctx) {
    if (ctx.runtime.mpi_size <= 1)
        return;
    auto* dd = ctx.runtime.dd;
    if (!dd || dd->cart_comm == MPI_COMM_NULL)
        return;

    // Cadence gate (see Reduce doc): with every_ > 1, fire only on step 0, the
    // final step of the run, and every every_-th step. Gated on the
    // all-ranks-identical ctx.step / is_last_step so the collective is still
    // entered in lockstep — no rank-divergent skip, hence no deadlock.
    // every_ == 1 (default) reduces every step, byte-identical to before; the
    // larger cadence is opt-in via MD::set_scalar_reduce_every. Always firing
    // on the last step keeps the globally-reduced scalars valid for the common
    // "read frame.scalars after run()" pattern, so the cadence is safe whenever
    // no *per-step* consumer reads them.
    if (every_ > 1 && !ctx.runtime.is_last_step && ctx.step != 0
        && (ctx.step % every_) != 0)
        return;

    auto& f = ctx.frame;
    // Pack {energy, energy_pair, energy_bond, energy_angle, energy_dihedral,
    // kinetic, stress[6], stress_kinetic[6]} into a single buffer for one
    // 18-double allreduce. The 4 per-term energy fields were added by
    // msu-pe-cuda-components-01 (Phase G.3); they reduce identically to
    // `energy` (each is a per-rank partial sum of one evaluator family).
    // `n_local` is omitted -- it is authoritative on each rank from
    // `dd->n_local` and does not need to be summed. Both stress fields
    // reduce identically (each is a per-rank partial sum of Voigt
    // components).
    double sendbuf[18];
    sendbuf[0] = f.scalars.energy;
    sendbuf[1] = f.scalars.energy_pair;
    sendbuf[2] = f.scalars.energy_bond;
    sendbuf[3] = f.scalars.energy_angle;
    sendbuf[4] = f.scalars.energy_dihedral;
    sendbuf[5] = f.scalars.kinetic;
    for (int k = 0; k < 6; ++k)
        sendbuf[6 + k] = f.scalars.stress[k];
    for (int k = 0; k < 6; ++k)
        sendbuf[12 + k] = f.scalars.stress_kinetic[k];

    double recvbuf[18];
    MPI_Allreduce(sendbuf, recvbuf, 18, MPI_DOUBLE, MPI_SUM, dd->cart_comm);

    f.scalars.energy = recvbuf[0];
    f.scalars.energy_pair = recvbuf[1];
    f.scalars.energy_bond = recvbuf[2];
    f.scalars.energy_angle = recvbuf[3];
    f.scalars.energy_dihedral = recvbuf[4];
    f.scalars.kinetic = recvbuf[5];
    for (int k = 0; k < 6; ++k)
        f.scalars.stress[k] = recvbuf[6 + k];
    for (int k = 0; k < 6; ++k)
        f.scalars.stress_kinetic[k] = recvbuf[12 + k];
}

#else // !ATV_HAS_MPI

void Reduce::launch(const cpu::LaunchContext& /*ctx*/) {}

#endif // ATV_HAS_MPI

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
