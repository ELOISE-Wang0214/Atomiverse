/// @file cpu/md/init.cpp
/// @brief CPU MD infrastructure component implementations.

#include "atomiverse/cpu/md/init.h"
#include "atomiverse/core/units.h"

#include <algorithm>
#include <cmath>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

// ── ZeroForces ───────────────────────────────────────────────────

void ZeroForces::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n = f.atoms.size();
    #pragma omp parallel for schedule(static)
    for (Int i = 0; i < n; ++i) {
        f.atoms.fx[i] = 0;
        f.atoms.fy[i] = 0;
        f.atoms.fz[i] = 0;
    }
}

// ── ZeroEnergy ───────────────────────────────────────────────────

void ZeroEnergy::launch(const cpu::LaunchContext& ctx) {
    ctx.frame.energy = 0;
}

// ── ZeroKinetic ──────────────────────────────────────────────────

void ZeroKinetic::launch(const cpu::LaunchContext& ctx) {
    ctx.frame.kinetic = 0;
}

// ── KineticEnergy ────────────────────────────────────────────────

void KineticEnergy::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    Int n = f.atoms.size();
    Real ke = 0;

    #pragma omp parallel for schedule(static) reduction(+:ke)
    for (Int i = 0; i < n; ++i) {
        Real v2 = f.atoms.vx[i] * f.atoms.vx[i] + f.atoms.vy[i] * f.atoms.vy[i] + f.atoms.vz[i] * f.atoms.vz[i];
        ke += f.atoms.mass[i] * v2;
    }

    f.kinetic = Real(0.5) * kMvv2e * ke;
}

// ── InitVelocity ─────────────────────────────────────────────────

InitVelocity::InitVelocity(Real temperature, unsigned seed)
    : temperature_(temperature), seed_(seed) {}

void InitVelocity::build(const cpu::BuildContext& ctx) {
    auto& f = ctx.frame;
    Int n = f.atoms.size();

    std::mt19937 rng(seed_);
    std::normal_distribution<Real> gauss(0, 1);

    for (Int i = 0; i < n; ++i) {
        Real sigma_v = std::sqrt(kBoltzmann * temperature_ / f.atoms.mass[i]) / std::sqrt(kMvv2e);
        f.atoms.vx[i] = gauss(rng) * sigma_v;
        f.atoms.vy[i] = gauss(rng) * sigma_v;
        f.atoms.vz[i] = gauss(rng) * sigma_v;
    }

    // Remove COM velocity
    Real vcm_x = 0, vcm_y = 0, vcm_z = 0, mtot = 0;
    for (Int i = 0; i < n; ++i) {
        vcm_x += f.atoms.mass[i] * f.atoms.vx[i];
        vcm_y += f.atoms.mass[i] * f.atoms.vy[i];
        vcm_z += f.atoms.mass[i] * f.atoms.vz[i];
        mtot += f.atoms.mass[i];
    }
    vcm_x /= mtot;
    vcm_y /= mtot;
    vcm_z /= mtot;

    for (Int i = 0; i < n; ++i) {
        f.atoms.vx[i] -= vcm_x;
        f.atoms.vy[i] -= vcm_y;
        f.atoms.vz[i] -= vcm_z;
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
