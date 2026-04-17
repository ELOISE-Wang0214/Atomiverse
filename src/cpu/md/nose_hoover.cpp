/// @file cpu/md/nose_hoover.cpp
/// @brief CPU Nose-Hoover thermostat implementation.

#include "atomiverse/cpu/md/nose_hoover.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

namespace {

inline int owned_count(const cpu::LaunchContext& ctx) {
    const auto* dd = ctx.runtime.dd;
    if (dd && ctx.runtime.mpi_size > 1) return dd->n_local;
    return ctx.frame.atoms.size();
}

} // namespace

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

        // Phase 2 limitation: local KE only. The thermostat is exact in
        // single-rank; under multi-rank it uses per-rank KE which diverges
        // from the true global bath -- a proper NVT multi-rank run requires
        // MPI_Allreduce of `ke` before the xi_ update. Tracked in Phase 3.
        // The Phase 2.7 exit-gate test uses NVE (no NoseHoover) so this
        // limitation is not on the critical path.
        double ke = 0;
        #pragma omp parallel for schedule(static) reduction(+:ke)
        for (int i = 0; i < n; ++i) {
            double v2 = f.atoms.vx[i] * f.atoms.vx[i] + f.atoms.vy[i] * f.atoms.vy[i] + f.atoms.vz[i] * f.atoms.vz[i];
            ke += f.atoms.mass[i] * v2;
        }
        ke *= double(0.5) * kMvv2e;

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
