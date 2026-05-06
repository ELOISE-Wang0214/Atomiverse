/// @file cpu/md/barostat.cpp
/// @brief CPU Berendsen isotropic barostat.

#include "atomiverse/cpu/md/barostat.h"

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/owned_range.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

Barostat::Barostat(double P_target, double tau_p, double dt)
    : P_target_(P_target), tau_p_(tau_p), dt_(dt) {}

void Barostat::launch(const cpu::LaunchContext& ctx) {
    if (!enabled_) return;
    auto& f = ctx.frame;

    // Volume from box (orthorhombic). For multi-rank the box is
    // identical across ranks, so V is the global volume.
    const double V = f.box(0, 0) * f.box(1, 1) * f.box(2, 2);
    if (V <= 0.0) return;

    // Pressure from the trace of the unnormalised virial sum:
    //   P = -(σ_xx + σ_yy + σ_zz) / (3·V)
    const double P =
        -(f.stress[0] + f.stress[1] + f.stress[2]) / (3.0 * V);

    // Berendsen scaling factor μ = (1 - dt/τ_p · (P - P_target))^(1/3).
    // Clamped to [0.5, 2.0] so a single step cannot collapse or explode
    // the box if P diverges (e.g. transient force spike on step 0).
    const double mu_cubed = 1.0 - (dt_ / tau_p_) * (P - P_target_);
    const double mu = std::cbrt(std::max(0.5, std::min(2.0, mu_cubed)));

    // Scale box edges. Single-rank or multi-rank: deterministic, identical
    // computation on every rank (P, V, μ are all global) so the box stays
    // consistent across ranks without any explicit broadcast.
    f.box(0, 0) *= mu;
    f.box(1, 1) *= mu;
    f.box(2, 2) *= mu;

    // Affine-rescale owned atom positions only. Ghost slots are
    // overwritten by the next step's HaloExchangeForward.
    const int n = owned_count(ctx);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.x[i] *= mu;
        f.atoms.y[i] *= mu;
        f.atoms.z[i] *= mu;
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
