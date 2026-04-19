/// @file cpu/md/init.cpp
/// @brief CPU MD infrastructure component implementations.

#include "atomiverse/cpu/md/init.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/md/mpi/domain.h"

#include <algorithm>
#include <cmath>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

namespace {

/// @brief Outer-loop bound for integrators and force-zeroers.
///
/// Multi-rank (`dd` non-null and `mpi_size > 1`): the authoritative count
/// of locally-owned atoms. Slots `[n_local, n_local + n_ghost)` hold ghost
/// copies that must not be touched by integrators or by the local
/// zeroing pass (`ZeroForces` must leave force contributions written
/// into ghost rows undisturbed so `HaloExchangeReverse` can ship them
/// back). Single-rank or DD-less contexts: fall back to `atoms.size()`.
inline int owned_count(const cpu::LaunchContext& ctx) {
    const auto* dd = ctx.runtime.dd;
    if (dd && ctx.runtime.mpi_size > 1) return dd->n_local;
    return ctx.frame.atoms.size();
}

inline int owned_count(const cpu::BuildContext& ctx) {
    const auto* dd = ctx.runtime.dd;
    if (dd && ctx.runtime.mpi_size > 1) return dd->n_local;
    return ctx.frame.atoms.size();
}

} // namespace

// ── ZeroForces ───────────────────────────────────────────────────

void ZeroForces::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    // Zero the full array (owners + every possible ghost slot up to capacity).
    // ZeroForces sits at order -1000 in pre_force, BEFORE HaloExchangeForward
    // populates ghosts for the current step, so we cannot rely on
    // `dd->n_ghost` here -- it still holds the previous step's count.
    // Zeroing [0, atoms.size()) is the cheapest correct option: it covers
    // every slot the next halo could write into, leaves owner forces clean
    // for force kernels, and gives `HaloExchangeReverse` a known-zero
    // baseline on ghost rows.
    const int n = f.atoms.size();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.fx[i] = 0;
        f.atoms.fy[i] = 0;
        f.atoms.fz[i] = 0;
    }
}

// ── ZeroEnergy ───────────────────────────────────────────────────

void ZeroEnergy::launch(const cpu::LaunchContext& ctx) {
    ctx.frame.energy = 0;
}

// ── ZeroStress ───────────────────────────────────────────────────

void ZeroStress::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    for (int k = 0; k < 6; ++k) f.stress[k] = 0;
}

// ── ZeroKinetic ──────────────────────────────────────────────────

void ZeroKinetic::launch(const cpu::LaunchContext& ctx) {
    ctx.frame.kinetic = 0;
}

// ── KineticEnergy ────────────────────────────────────────────────

void KineticEnergy::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);  // owned only -- ghosts duplicate owners' KE
    double ke = 0;

    #pragma omp parallel for schedule(static) reduction(+:ke)
    for (int i = 0; i < n; ++i) {
        double v2 = f.atoms.vx[i] * f.atoms.vx[i] + f.atoms.vy[i] * f.atoms.vy[i] + f.atoms.vz[i] * f.atoms.vz[i];
        ke += f.atoms.mass[i] * v2;
    }

    f.kinetic = double(0.5) * kMvv2e * ke;
}

// ── InitVelocity ─────────────────────────────────────────────────

InitVelocity::InitVelocity(double temperature, unsigned seed)
    : temperature_(temperature), seed_(seed) {}

void InitVelocity::build(const cpu::BuildContext& ctx) {
    // Phase 2: InitVelocity operates on locally-owned atoms only. COM
    // removal is *local* -- for multi-rank usage that requires a
    // globally-zero COM, the caller should set velocities pre-partition
    // on a single rank and let `Backend::set_frame()` partition the
    // already-initialised velocities. A proper multi-rank InitVelocity
    // with MPI_Allreduce-based global COM is tracked in Phase 3.
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    std::mt19937 rng(seed_);
    std::normal_distribution<double> gauss(0, 1);

    for (int i = 0; i < n; ++i) {
        double sigma_v = std::sqrt(kBoltzmann * temperature_ / f.atoms.mass[i]) / std::sqrt(kMvv2e);
        f.atoms.vx[i] = gauss(rng) * sigma_v;
        f.atoms.vy[i] = gauss(rng) * sigma_v;
        f.atoms.vz[i] = gauss(rng) * sigma_v;
    }

    // Local COM velocity removal.
    double vcm_x = 0, vcm_y = 0, vcm_z = 0, mtot = 0;
    for (int i = 0; i < n; ++i) {
        vcm_x += f.atoms.mass[i] * f.atoms.vx[i];
        vcm_y += f.atoms.mass[i] * f.atoms.vy[i];
        vcm_z += f.atoms.mass[i] * f.atoms.vz[i];
        mtot += f.atoms.mass[i];
    }
    if (mtot > 0) {
        vcm_x /= mtot;
        vcm_y /= mtot;
        vcm_z /= mtot;
        for (int i = 0; i < n; ++i) {
            f.atoms.vx[i] -= vcm_x;
            f.atoms.vy[i] -= vcm_y;
            f.atoms.vz[i] -= vcm_z;
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
