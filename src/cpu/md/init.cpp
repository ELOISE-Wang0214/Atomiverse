/// @file cpu/md/init.cpp
/// @brief CPU MD infrastructure component implementations.

#include "atomiverse/cpu/md/init.h"

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/owned_range.h"

#include <algorithm>
#include <cmath>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

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

void ZeroEnergy::launch(const cpu::LaunchContext& ctx) { ctx.frame.energy = 0; }

// ── ZeroStress ───────────────────────────────────────────────────

void ZeroStress::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    for (int k = 0; k < 6; ++k)
        f.stress[k] = 0;
    // `stress_kinetic` is overwritten (not accumulated) by VirialStress
    // each step, so it does not strictly need zeroing -- but zero it
    // here for symmetry and so observers that read it before VirialStress
    // ever ran see a deterministic zero rather than uninitialised data.
    for (int k = 0; k < 6; ++k)
        f.stress_kinetic[k] = 0;
}

// ── ZeroKinetic ──────────────────────────────────────────────────

void ZeroKinetic::launch(const cpu::LaunchContext& ctx) {
    ctx.frame.kinetic = 0;
}

// ── KineticEnergy ────────────────────────────────────────────────

void KineticEnergy::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx); // owned only -- ghosts duplicate owners' KE
    double ke = 0;

#pragma omp parallel for schedule(static) reduction(+ : ke)
    for (int i = 0; i < n; ++i) {
        double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                    + f.atoms.vy[i] * f.atoms.vy[i]
                    + f.atoms.vz[i] * f.atoms.vz[i];
        ke += f.atoms.mass[i] * v2;
    }

    f.kinetic = double(0.5) * kMvv2e * ke;
}

// ── ZeroLinearMomentum ───────────────────────────────────────────

void ZeroLinearMomentum::launch(const cpu::LaunchContext& ctx) {
    if (step_++ % freq_ != 0)
        return;

    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    double p[3] = {0.0, 0.0, 0.0};
    double mtot = 0.0;
    for (int i = 0; i < n; ++i) {
        const double m = f.atoms.mass[i];
        p[0] += m * f.atoms.vx[i];
        p[1] += m * f.atoms.vy[i];
        p[2] += m * f.atoms.vz[i];
        mtot += m;
    }

#ifdef ATV_HAS_MPI
    if (ctx.runtime.mpi_size > 1 && ctx.runtime.dd
        && ctx.runtime.dd->cart_comm != MPI_COMM_NULL) {
        double local[4] = {p[0], p[1], p[2], mtot};
        double global[4];
        MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_SUM,
                      ctx.runtime.dd->cart_comm);
        p[0] = global[0];
        p[1] = global[1];
        p[2] = global[2];
        mtot = global[3];
    }
#endif
    if (mtot <= 0.0)
        return;

    const double vx = p[0] / mtot;
    const double vy = p[1] / mtot;
    const double vz = p[2] / mtot;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.vx[i] -= vx;
        f.atoms.vy[i] -= vy;
        f.atoms.vz[i] -= vz;
    }
}

// ── InitVelocity ─────────────────────────────────────────────────

InitVelocity::InitVelocity(double temperature, unsigned seed, VelocityDist dist)
    : temperature_(temperature), seed_(seed), dist_(dist) {}

void InitVelocity::build(const cpu::BuildContext& ctx) {
    // Mirrors LAMMPS `Velocity::create()` (velocity.cpp:160):
    //   1. draw v from `dist` (uniform / gaussian), scale by 1/sqrt(m);
    //   2. remove COM velocity (`zero_momentum()`, velocity.cpp:738);
    //   3. rescale so the measured temperature equals the requested T
    //      (`rescale(t, t_desired)`, velocity.cpp:368).
    //
    // Phase 2: InitVelocity operates on locally-owned atoms only. COM
    // removal is *local* -- for multi-rank usage that requires a
    // globally-zero COM, the caller should set velocities pre-partition
    // on a single rank and let `Backend::set_frame()` partition the
    // already-initialised velocities. The temperature rescale below is
    // MPI-aware (it reduces the global KE) so the measured *global*
    // temperature is exact even though the COM removal stays local. A
    // proper multi-rank InitVelocity with MPI_Allreduce-based global COM
    // is tracked in Phase 3.
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    // ── 1. draw velocities from the chosen distribution ──────────────
    // LAMMPS scales the raw draw by `1/sqrt(m)` (velocity.cpp:287-291);
    // the absolute scale is then fixed by the temperature rescale in
    // step 3, so the per-component prefactor only needs to be consistent.
    std::mt19937 rng(seed_);
    if (dist_ == VelocityDist::gaussian) {
        std::normal_distribution<double> draw(0.0, 1.0); // N(0,1)
        for (int i = 0; i < n; ++i) {
            const double s = 1.0 / std::sqrt(f.atoms.mass[i]);
            f.atoms.vx[i] = draw(rng) * s;
            f.atoms.vy[i] = draw(rng) * s;
            f.atoms.vz[i] = draw(rng) * s;
        }
    } else { // VelocityDist::uniform
        // LAMMPS draws `random->uniform() - 0.5` → U(-0.5, 0.5).
        std::uniform_real_distribution<double> draw(-0.5, 0.5);
        for (int i = 0; i < n; ++i) {
            const double s = 1.0 / std::sqrt(f.atoms.mass[i]);
            f.atoms.vx[i] = draw(rng) * s;
            f.atoms.vy[i] = draw(rng) * s;
            f.atoms.vz[i] = draw(rng) * s;
        }
    }

    // ── 2. local COM velocity removal (LAMMPS `mom yes` default) ─────
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

    // ── 3. rescale so the measured temperature equals T exactly ──────
    // LAMMPS `rescale(t_old, t_new)`: factor = sqrt(t_new / t_old).
    // KE = 0.5 * kMvv2e * sum(m v^2); T = 2*KE / (Ndof * kB) with
    // Ndof = 3N - 3 (COM removed), matching `compute temp` extra_dof.
    double sum_mv2 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                          + f.atoms.vy[i] * f.atoms.vy[i]
                          + f.atoms.vz[i] * f.atoms.vz[i];
        sum_mv2 += f.atoms.mass[i] * v2;
    }
    long long n_global = n;
#ifdef ATV_HAS_MPI
    if (ctx.runtime.mpi_size > 1 && ctx.runtime.dd
        && ctx.runtime.dd->cart_comm != MPI_COMM_NULL) {
        double local = sum_mv2, global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM,
                      ctx.runtime.dd->cart_comm);
        sum_mv2 = global;
        long long n_local = n, ng = 0;
        MPI_Allreduce(&n_local, &ng, 1, MPI_LONG_LONG, MPI_SUM,
                      ctx.runtime.dd->cart_comm);
        n_global = ng;
    }
#endif
    const long long ndof = 3 * n_global - 3;
    if (ndof > 0 && sum_mv2 > 0.0) {
        const double ke = 0.5 * kMvv2e * sum_mv2;
        const double t_measured =
            2.0 * ke / (static_cast<double>(ndof) * kBoltzmann);
        const double factor = std::sqrt(temperature_ / t_measured);
        for (int i = 0; i < n; ++i) {
            f.atoms.vx[i] *= factor;
            f.atoms.vy[i] *= factor;
            f.atoms.vz[i] *= factor;
        }
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
