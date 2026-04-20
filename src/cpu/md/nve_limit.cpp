/// @file cpu/md/nve_limit.cpp
/// @brief CPU velocity-Verlet with per-atom velocity cap.

#include "atomiverse/cpu/md/nve_limit.h"

#include "atomiverse/core/units.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/facade/stage.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace md {

namespace {

inline int owned_count(const cpu::LaunchContext& ctx) {
    const auto* dd = ctx.runtime.dd;
    if (dd && ctx.runtime.mpi_size > 1) return dd->n_local;
    return ctx.frame.atoms.size();
}

inline void wrap_pbc(cpu::Frame& f, int n_owned) {
    const double lx = f.box(0, 0), ly = f.box(1, 1), lz = f.box(2, 2);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_owned; ++i) {
        f.atoms.x[i] -= lx * std::floor(f.atoms.x[i] / lx);
        f.atoms.y[i] -= ly * std::floor(f.atoms.y[i] / ly);
        f.atoms.z[i] -= lz * std::floor(f.atoms.z[i] / lz);
    }
}

} // namespace

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
        wrap_pbc(f, n);

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
