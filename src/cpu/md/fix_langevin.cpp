/// @file cpu/md/fix_langevin.cpp
/// @brief LAMMPS-equivalent Langevin thermostat (force-additive form).

#include "atomiverse/cpu/md/fix_langevin.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/md/mpi/domain.h"

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

// Box-Muller Gaussian, unit variance, from a splitmix64 stream.
struct GaussStream {
    std::uint64_t s;
    double next_u01() {
        s += 0x9E3779B97F4A7C15ULL;
        std::uint64_t z = s;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z =  z ^ (z >> 31);
        return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
    }
    double next() {
        const double u = next_u01();
        const double v = next_u01();
        const double r = std::sqrt(-2.0 * std::log(u > 1e-300 ? u : 1e-300));
        return r * std::cos(2.0 * 3.14159265358979323846 * v);
    }
};

} // namespace

FixLangevin::FixLangevin(double T_target, double damp,
                         double dt, std::uint32_t seed)
    : T_target_(T_target), damp_(damp), dt_(dt), seed_(seed) {}

void FixLangevin::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    // Per-step seed mixed in. Atom identity = global tag, so the same
    // atom on different rank counts gets the same kick (modulo
    // floating-point reduction noise elsewhere).
    const std::uint64_t base =
        static_cast<std::uint64_t>(seed_) * 0x9E3779B97F4A7C15ULL
        + static_cast<std::uint64_t>(ctx.step) * 0xBF58476D1CE4E5B9ULL;

    // gfactor1 multiplies v -> drag force.
    // gfactor2_per_sqrt_mass * sqrt(m) * gauss -> random force.
    const double inv_damp_ftm2v = 1.0 / (damp_ * kFtm2v);
    const double gfactor2_pre =
        std::sqrt(2.0 * kBoltzmann * T_target_
                  / (damp_ * dt_ * kMvv2e)) / kFtm2v;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        GaussStream g{base + static_cast<std::uint64_t>(f.atoms.tag[i])
                              * 0xD6E8FEB86659FD93ULL};

        const double m = f.atoms.mass[i];
        const double gfactor1 = -m * inv_damp_ftm2v;
        const double gfactor2 = gfactor2_pre * std::sqrt(m);

        const double rx = g.next();
        const double ry = g.next();
        const double rz = g.next();

        f.atoms.fx[i] += gfactor1 * f.atoms.vx[i] + gfactor2 * rx;
        f.atoms.fy[i] += gfactor1 * f.atoms.vy[i] + gfactor2 * ry;
        f.atoms.fz[i] += gfactor1 * f.atoms.vz[i] + gfactor2 * rz;
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
