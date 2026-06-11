/// @file cpu/md/nve_limit.cpp
/// @brief CPU velocity-Verlet with a per-step displacement (distance)
///        cap. LAMMPS `fix nve/limit` (`lammps/src/fix_nve_limit.cpp`).
///
/// Audit verdict (spec lammps-api-align-02-nve): ALIGNED. The cap below
/// reproduces `fix_nve_limit.cpp` exactly — `vmax = xmax/dt` mirrors
/// `vlimitsq = (xlimit/dtv)^2` (L66), the Euclidean `vsq = v·v` test
/// matches L105/L168, the rescale `s = vmax/sqrt(vsq)` is algebraically
/// identical to LAMMPS `scale = sqrt(vlimitsq/vsq)`, and the cap is
/// applied in BOTH the `initial_integrate` and `final_integrate`
/// branches (matching `FixNVELimit::initial_integrate`/`final_integrate`).

#include "atomiverse/cpu/md/nve_limit.h"

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/owned_range.h"
#include "atomiverse/cpu/md/pbc.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace md {

NVELimit::NVELimit(double dt, double xmax) : dt_(dt), xmax_(xmax) {}

void NVELimit::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const double vmax = xmax_ / dt_;
    const double vmax_sq = vmax * vmax;
    const double dtf = 0.5 * dt_ * kFtm2v;

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            const double dtfm = dtf / f.atoms.mass[i];
            f.atoms.vx[i] += dtfm * f.atoms.fx[i];
            f.atoms.vy[i] += dtfm * f.atoms.fy[i];
            f.atoms.vz[i] += dtfm * f.atoms.fz[i];

            const double vsq = f.atoms.vx[i] * f.atoms.vx[i]
                               + f.atoms.vy[i] * f.atoms.vy[i]
                               + f.atoms.vz[i] * f.atoms.vz[i];
            if (vsq > vmax_sq) {
                const double s = vmax / std::sqrt(vsq);
                f.atoms.vx[i] *= s;
                f.atoms.vy[i] *= s;
                f.atoms.vz[i] *= s;
            }
            f.atoms.x[i] += dt_ * f.atoms.vx[i];
            f.atoms.y[i] += dt_ * f.atoms.vy[i];
            f.atoms.z[i] += dt_ * f.atoms.vz[i];
        }
        wrap_orthorhombic(f, n);

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            const double dtfm = dtf / f.atoms.mass[i];
            f.atoms.vx[i] += dtfm * f.atoms.fx[i];
            f.atoms.vy[i] += dtfm * f.atoms.fy[i];
            f.atoms.vz[i] += dtfm * f.atoms.fz[i];

            const double vsq = f.atoms.vx[i] * f.atoms.vx[i]
                               + f.atoms.vy[i] * f.atoms.vy[i]
                               + f.atoms.vz[i] * f.atoms.vz[i];
            if (vsq > vmax_sq) {
                const double s = vmax / std::sqrt(vsq);
                f.atoms.vx[i] *= s;
                f.atoms.vy[i] *= s;
                f.atoms.vz[i] *= s;
            }
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
