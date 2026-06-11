/// @file cpu/md/velocity_verlet.cpp
/// @brief CPU VelocityVerlet implementation with OpenMP.

#include "atomiverse/cpu/md/velocity_verlet.h"

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/component/owned_range.h"
#include "atomiverse/cpu/md/pbc.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace md {

// ── VelocityVerlet ───────────────────────────────────────────────

VelocityVerlet::VelocityVerlet(double dt) : dt_(dt) {}

void VelocityVerlet::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        const double dtf = double(0.5) * dt_ * kFtm2v;

#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            double dtfm = dtf / f.atoms.mass[i];
            f.atoms.vx[i] += dtfm * f.atoms.fx[i];
            f.atoms.vy[i] += dtfm * f.atoms.fy[i];
            f.atoms.vz[i] += dtfm * f.atoms.fz[i];
            f.atoms.x[i] += dt_ * f.atoms.vx[i];
            f.atoms.y[i] += dt_ * f.atoms.vy[i];
            f.atoms.z[i] += dt_ * f.atoms.vz[i];
        }

        wrap_orthorhombic(f, n);

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
        const double dtf = double(0.5) * dt_ * kFtm2v;

#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            double dtfm = dtf / f.atoms.mass[i];
            f.atoms.vx[i] += dtfm * f.atoms.fx[i];
            f.atoms.vy[i] += dtfm * f.atoms.fy[i];
            f.atoms.vz[i] += dtfm * f.atoms.fz[i];
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
