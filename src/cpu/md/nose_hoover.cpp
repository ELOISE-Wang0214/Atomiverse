/// @file cpu/md/nose_hoover.cpp
/// @brief CPU Nose-Hoover thermostat implementation.

#include "atomiverse/cpu/md/nose_hoover.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/owned_range.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

NoseHoover::NoseHoover(double T_target, double tau, double dt)
    : T_target_(T_target), tau_(tau), dt_(dt) {}

void NoseHoover::build(const cpu::BuildContext& ctx) {
    // `atoms.size()` at build time is the *global* pre-partition atom count
    // (the storage is allocated once by the user before Backend::set_frame
    // partitions). This gives the correct global `n_dof` for `Q_`.
    int n_global = ctx.frame.atoms.size();
    int n_dof = 3 * n_global - 3;  // 3N - 3 (COM removal)
    Q_ = static_cast<double>(n_dof) * kBoltzmann * T_target_ * tau_ * tau_;
}

void NoseHoover::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    // Reconstruct global DOF from the same factor used in build().
    const int n_global =
        (ctx.runtime.dd && ctx.runtime.mpi_size > 1) ? f.atoms.size()
                                                     : n;
    const int n_dof = 3 * n_global - 3;

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate) ||
        ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {

        // Per-rank kinetic over owned atoms only.
        double ke_local = 0;
        #pragma omp parallel for schedule(static) reduction(+:ke_local)
        for (int i = 0; i < n; ++i) {
            double v2 = f.atoms.vx[i] * f.atoms.vx[i] + f.atoms.vy[i] * f.atoms.vy[i] + f.atoms.vz[i] * f.atoms.vz[i];
            ke_local += f.atoms.mass[i] * v2;
        }
        ke_local *= double(0.5) * kMvv2e;

        // Phase 0.5: under multi-rank reduce KE across the Cartesian comm so
        // the thermostat couples to the *global* bath rather than this
        // rank's slab. Adds two MPI_Allreduces per integration step (one in
        // initial_integrate, one in final_integrate); the cost is a single
        // double per call.
        double ke = ke_local;
#ifdef ATV_HAS_MPI
        if (ctx.runtime.mpi_size > 1 && ctx.runtime.dd
            && ctx.runtime.dd->cart_comm != MPI_COMM_NULL) {
            MPI_Allreduce(&ke_local, &ke, 1, MPI_DOUBLE, MPI_SUM,
                          ctx.runtime.dd->cart_comm);
        }
#endif

        // Update thermostat variable
        xi_ += (dt_ * double(0.5)) * (double(2) * ke - static_cast<double>(n_dof) * kBoltzmann * T_target_) / Q_;

        // Rescale velocities (owners only).
        double scale = std::exp(-dt_ * double(0.5) * xi_);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            f.atoms.vx[i] *= scale;
            f.atoms.vy[i] *= scale;
            f.atoms.vz[i] *= scale;
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
