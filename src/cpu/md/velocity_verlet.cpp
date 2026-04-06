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
    Int n = f.atoms.size();
    Real lx = f.box(0, 0), ly = f.box(1, 1), lz = f.box(2, 2);

    #pragma omp parallel for schedule(static)
    for (Int i = 0; i < n; ++i) {
        f.atoms.x[i] -= lx * std::floor(f.atoms.x[i] / lx);
        f.atoms.y[i] -= ly * std::floor(f.atoms.y[i] / ly);
        f.atoms.z[i] -= lz * std::floor(f.atoms.z[i] / lz);
    }
}

// ── VelocityVerlet ───────────────────────────────────────────────

VelocityVerlet::VelocityVerlet(Real dt) : dt_(dt) {}

void VelocityVerlet::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n = f.atoms.size();

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        const Real dtf = Real(0.5) * dt_ * kFtm2v;

        #pragma omp parallel for schedule(static)
        for (Int i = 0; i < n; ++i) {
            Real dtfm = dtf / f.atoms.mass[i];
            f.atoms.vx[i] += dtfm * f.atoms.fx[i];
            f.atoms.vy[i] += dtfm * f.atoms.fy[i];
            f.atoms.vz[i] += dtfm * f.atoms.fz[i];
            f.atoms.x[i] += dt_ * f.atoms.vx[i];
            f.atoms.y[i] += dt_ * f.atoms.vy[i];
            f.atoms.z[i] += dt_ * f.atoms.vz[i];
        }

        wrap_pbc(f);

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
        const Real dtf = Real(0.5) * dt_ * kFtm2v;

        #pragma omp parallel for schedule(static)
        for (Int i = 0; i < n; ++i) {
            Real dtfm = dtf / f.atoms.mass[i];
            f.atoms.vx[i] += dtfm * f.atoms.fx[i];
            f.atoms.vy[i] += dtfm * f.atoms.fy[i];
            f.atoms.vz[i] += dtfm * f.atoms.fz[i];
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
