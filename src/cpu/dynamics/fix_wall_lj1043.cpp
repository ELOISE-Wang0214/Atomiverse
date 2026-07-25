/// @file src/cpu/dynamics/fix_wall_lj1043.cpp
#include "atomiverse/cpu/dynamics/fix_wall_lj1043.h"

#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/owned_range.h"

namespace ATV {
namespace cpu {
namespace md {

namespace {

constexpr double kTwoPi = 6.28318530717958647692528676655900576;

// Positive integer powers only (matches LAMMPS MathSpecial::powint usage).
double pow_pos(double x, int n) {
    double r = 1.0;
    for (int i = 0; i < n; ++i)
        r *= x;
    return r;
}

} // namespace

FixWallLJ1043::FixWallLJ1043(int axis_, int side_, double wall_pos_,
                             double epsilon_, double sigma_, double cutoff_)
    : axis(axis_), side(side_), wall_pos(wall_pos_), epsilon(epsilon_),
      sigma(sigma_), cutoff(cutoff_) {
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("FixWallLJ1043: axis must be 0, 1, or 2");
    if (!(side == -1 || side == 1))
        throw std::invalid_argument(
            "FixWallLJ1043: side must be -1 (lo) or +1 (hi)");
    if (!(epsilon > 0.0))
        throw std::invalid_argument("FixWallLJ1043: epsilon must be > 0");
    if (!(sigma > 0.0))
        throw std::invalid_argument("FixWallLJ1043: sigma must be > 0");
    if (!(cutoff > 0.0))
        throw std::invalid_argument("FixWallLJ1043: cutoff must be > 0");
    precompute_();
}

void FixWallLJ1043::precompute_() {
    // Match LAMMPS FixWallLJ1043::precompute exactly.
    const double s = sigma;
    const double e = epsilon;
    const double sqrt2 = std::sqrt(2.0);
    coeff1_ = kTwoPi * 2.0 / 5.0 * e * pow_pos(s, 10);
    coeff2_ = kTwoPi * e * pow_pos(s, 4);
    coeff3_ = kTwoPi * sqrt2 / 3.0 * e * pow_pos(s, 3);
    coeff4_ = 0.61 / sqrt2 * s;
    coeff5_ = coeff1_ * 10.0;
    coeff6_ = coeff2_ * 4.0;
    coeff7_ = coeff3_ * 3.0;

    const double rc = cutoff;
    const double rinv = 1.0 / rc;
    const double r2inv = rinv * rinv;
    const double r4inv = r2inv * r2inv;
    // offset = coeff1/rc^10 - coeff2/rc^4 - coeff3/(rc+coeff4)^3
    const double rpc = rc + coeff4_;
    offset_ = coeff1_ * r4inv * r4inv * r2inv - coeff2_ * r4inv
              - coeff3_ / (rpc * rpc * rpc);
}

void FixWallLJ1043::launch(const cpu::LaunchContext& ctx) {
    // Honour public-member edits between launches.
    precompute_();

    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const int dim = axis;
    const int sgn = side; // -1 lo, +1 hi
    const double coord = wall_pos;
    const double rc = cutoff;

    const StorageReal* __restrict x = (dim == 0)   ? f.span(kAtom_x)
                                      : (dim == 1) ? f.span(kAtom_y)
                                                   : f.span(kAtom_z);
    AccReal* __restrict force = (dim == 0)   ? f.span(kAtom_fx)
                                : (dim == 1) ? f.span(kAtom_fy)
                                             : f.span(kAtom_fz);

    double u = 0.0;
    const bool want_force = ctx.forces_requested;

#pragma omp parallel for reduction(+ : u) schedule(static)
    for (int i = 0; i < n; ++i) {
        // delta = distance into the fluid half-space (LAMMPS wall_particle).
        const double xi = static_cast<double>(x[i]);
        double delta;
        if (sgn < 0)
            delta = xi - coord; // lo: atom above wall
        else
            delta = coord - xi; // hi: atom below wall

        if (delta <= 0.0)
            continue; // on or behind wall — skip (foothold: no hard error)
        if (delta >= rc)
            continue;

        const double rinv = 1.0 / delta;
        const double r2inv = rinv * rinv;
        const double r4inv = r2inv * r2inv;
        const double r10inv = r4inv * r4inv * r2inv;
        const double rp = delta + coeff4_;
        const double rp_inv = 1.0 / rp;
        const double rp3_inv = rp_inv * rp_inv * rp_inv;
        const double rp4_inv = rp3_inv * rp_inv;

        // fwall = side * (coeff5/r^11 - coeff6/r^5 - coeff7/(r+c4)^4)
        const double fwall = static_cast<double>(sgn)
                             * (coeff5_ * r10inv * rinv - coeff6_ * r4inv * rinv
                                - coeff7_ * rp4_inv);
        if (want_force)
            force[i] -= fwall;

        // U = coeff1/r^10 - coeff2/r^4 - coeff3/(r+c4)^3 - offset
        u += coeff1_ * r10inv - coeff2_ * r4inv - coeff3_ * rp3_inv - offset_;
    }

    last_u_ = u;
    f.meta.get(kEnergy) += u;
}

} // namespace md
} // namespace cpu
} // namespace ATV
