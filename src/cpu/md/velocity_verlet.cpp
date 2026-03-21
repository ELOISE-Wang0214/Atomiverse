/// @file cpu/md/velocity_verlet.cpp
/// @brief CPU VelocityVerlet implementation with OpenMP.

#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace md {

// ── PBC wrapping (orthorhombic) ──────────────────────────────────

static void wrap_pbc(cpu::Frame& f) {
    Int n = f.n_atoms();
    Real lx = f.box(0, 0), ly = f.box(1, 1), lz = f.box(2, 2);

    #pragma omp parallel for schedule(static)
    for (Int i = 0; i < n; ++i) {
        f.x[i] -= lx * std::floor(f.x[i] / lx);
        f.y[i] -= ly * std::floor(f.y[i] / ly);
        f.z[i] -= lz * std::floor(f.z[i] / lz);
    }
}

// ── VelocityVerlet ───────────────────────────────────────────────

VelocityVerlet::VelocityVerlet(Real dt) : dt_(dt) {}

void VelocityVerlet::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n = f.n_atoms();

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        const Real dtf = Real(0.5) * dt_ * kFtm2v;

        #pragma omp parallel for schedule(static)
        for (Int i = 0; i < n; ++i) {
            Real dtfm = dtf / f.mass[i];
            f.vx[i] += dtfm * f.fx[i];
            f.vy[i] += dtfm * f.fy[i];
            f.vz[i] += dtfm * f.fz[i];
            f.x[i] += dt_ * f.vx[i];
            f.y[i] += dt_ * f.vy[i];
            f.z[i] += dt_ * f.vz[i];
        }

        wrap_pbc(f);

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
        const Real dtf = Real(0.5) * dt_ * kFtm2v;

        #pragma omp parallel for schedule(static)
        for (Int i = 0; i < n; ++i) {
            Real dtfm = dtf / f.mass[i];
            f.vx[i] += dtfm * f.fx[i];
            f.vy[i] += dtfm * f.fy[i];
            f.vz[i] += dtfm * f.fz[i];
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
