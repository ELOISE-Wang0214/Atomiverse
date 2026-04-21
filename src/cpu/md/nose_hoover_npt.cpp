/// @file cpu/md/nose_hoover_npt.cpp
/// @brief Single-chain Martyna-Tobias-Klein isotropic NPT.

#include "atomiverse/cpu/md/nose_hoover_npt.h"

#include "atomiverse/core/units.h"
#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/md/mpi/domain.h"
#include "atomiverse/cpu/md/owned_range.h"
#include "atomiverse/cpu/md/pbc.h"

#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef ATV_HAS_MPI
#include <mpi.h>
#endif

namespace ATV {
namespace cpu {
namespace md {

namespace {

/// @brief Hard clamp on `|ε̇·dt|` used by the box-scaling update.
///
/// `ε̇` is the strain-rate log-volume velocity; the exact update scales
/// the box by `exp(ε̇·dt)` each half-step. A transient large `ε̇`
/// (early equilibration, large pressure error, user-chosen small
/// `τ_P`) can blow the box up or collapse it to zero in a single
/// step. Clamping the stride before `exp` matches the LAMMPS `fix
/// npt iso` guard and keeps the integrator stable without materially
/// affecting the converged NPT ensemble (the clamp only bites during
/// the first few steps of a poorly-initialised run).
constexpr double kEpsDtMaxFraction = 0.1;

// sinh(x)/x with stable Taylor branch near 0.
inline double sinhc(double x) {
    if (std::abs(x) < 1e-6) {
        const double x2 = x * x;
        return 1.0 + x2 / 6.0 + x2 * x2 / 120.0;
    }
    return std::sinh(x) / x;
}

} // namespace

NoseHooverNPT::NoseHooverNPT(Config cfg)
    : T_target_(cfg.T_target), tau_T_(cfg.tau_T),
      P_target_(cfg.P_target), tau_P_(cfg.tau_P), dt_(cfg.dt) {}

void NoseHooverNPT::build(const cpu::BuildContext& ctx) {
    // atoms.size() at build time = global N (set by user before
    // Backend::set_frame partitions).
    const int N = static_cast<int>(ctx.frame.atoms.size());
    n_dof_ = 3 * N - 3;
    kT_    = kBoltzmann * T_target_;

    // Standard MTK masses (Tuckerman 2006 §III).
    Q_t_ = static_cast<double>(n_dof_) * kT_ * tau_T_ * tau_T_;
    Q_b_ = kT_ * tau_T_ * tau_T_;
    W_   = (static_cast<double>(n_dof_) + 3.0) * kT_ * tau_P_ * tau_P_;

    // Initial KE at equipartition — used in step 0's first NHC update
    // before any velocity update has produced a real KE_global_.
    KE_global_ = 0.5 * static_cast<double>(n_dof_) * kT_;
}

double NoseHooverNPT::compute_global_KE(const cpu::LaunchContext& ctx) const {
    const auto& f = ctx.frame;
    const int n   = owned_count(ctx);
    double ke_local = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:ke_local)
    for (int i = 0; i < n; ++i) {
        const double v2 = f.atoms.vx[i] * f.atoms.vx[i]
                        + f.atoms.vy[i] * f.atoms.vy[i]
                        + f.atoms.vz[i] * f.atoms.vz[i];
        ke_local += f.atoms.mass[i] * v2;
    }
    ke_local *= 0.5 * kMvv2e;

    double ke = ke_local;
#ifdef ATV_HAS_MPI
    if (ctx.runtime.mpi_size > 1 && ctx.runtime.dd
        && ctx.runtime.dd->cart_comm != MPI_COMM_NULL) {
        MPI_Allreduce(&ke_local, &ke, 1, MPI_DOUBLE, MPI_SUM,
                      ctx.runtime.dd->cart_comm);
    }
#endif
    return ke;
}

void NoseHooverNPT::nhc_atoms_half(const cpu::LaunchContext& ctx, double half_dt) {
    // Tuckerman 2006 single-chain step.
    const double quarter_dt = 0.5 * half_dt;
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    // --- update eta_t_dot using current KE ---
    double G = (2.0 * KE_global_ - n_dof_ * kT_) / Q_t_;
    eta_t_dot_ += G * quarter_dt;

    // --- scale velocities ---
    const double scale = std::exp(-eta_t_dot_ * half_dt);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.vx[i] *= scale;
        f.atoms.vy[i] *= scale;
        f.atoms.vz[i] *= scale;
    }
    KE_global_ *= scale * scale;

    // --- propagate eta_t (informational, no force on it) ---
    eta_t_ += eta_t_dot_ * half_dt;

    // --- recompute G with new KE, update eta_t_dot again ---
    G = (2.0 * KE_global_ - n_dof_ * kT_) / Q_t_;
    eta_t_dot_ += G * quarter_dt;
}

void NoseHooverNPT::nhc_baro_half(double half_dt) {
    const double quarter_dt = 0.5 * half_dt;

    double G = (W_ * eps_dot_ * eps_dot_ - kT_) / Q_b_;
    eta_b_dot_ += G * quarter_dt;

    const double scale = std::exp(-eta_b_dot_ * half_dt);
    eps_dot_ *= scale;

    eta_b_ += eta_b_dot_ * half_dt;

    G = (W_ * eps_dot_ * eps_dot_ - kT_) / Q_b_;
    eta_b_dot_ += G * quarter_dt;
}

void NoseHooverNPT::update_eps_dot(const cpu::LaunchContext& ctx, double half_dt) {
    const auto& f = ctx.frame;
    const double V = f.box(0,0) * f.box(1,1) * f.box(2,2);
    if (V <= 0.0) return;

    // P_internal in `frame.stress` convention (matches Berendsen barostat
    // and tests/mpi/test_npt_multirank.cpp).  At end_of_step time
    // `frame.stress` is reduced; at initial_integrate time it carries
    // the previous step's reduced value (a "reasonable stale" estimate).
    const double P_int = -(f.stress[0] + f.stress[1] + f.stress[2]) / (3.0 * V);

    // G_eps drives box growth when P_target > P_int (need to expand to
    // hit lower target) AND when KE is high (kinetic-stress contribution).
    // The (P_target - P_int) sign matches the existing Berendsen
    // convention because both compute pressure as `-trace/(3V)`.
    const double G_eps = (3.0 * V * (P_target_ - P_int)
                          + (6.0 / static_cast<double>(n_dof_)) * KE_global_)
                         / W_;
    eps_dot_ += G_eps * half_dt;
}

void NoseHooverNPT::scale_v_by_baro(const cpu::LaunchContext& ctx, double half_dt) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const double factor = -(1.0 + 3.0 / static_cast<double>(n_dof_)) * eps_dot_;
    const double scale  = std::exp(factor * half_dt);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.vx[i] *= scale;
        f.atoms.vy[i] *= scale;
        f.atoms.vz[i] *= scale;
    }
    KE_global_ *= scale * scale;
}

void NoseHooverNPT::kick_v_by_force(const cpu::LaunchContext& ctx, double half_dt) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const double dtf = half_dt * kFtm2v;

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        const double dtfm = dtf / f.atoms.mass[i];
        f.atoms.vx[i] += dtfm * f.atoms.fx[i];
        f.atoms.vy[i] += dtfm * f.atoms.fy[i];
        f.atoms.vz[i] += dtfm * f.atoms.fz[i];
    }
    // KE invalidated; will be recomputed at end_of_step.
}

void NoseHooverNPT::propagate_x_and_box(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);

    // Cap per-step box scaling to ~10% to avoid integer overflow in cell-list
    // grid sizing if eps_dot diverges (NaN/inf) from a transient stress
    // spike. Once clamped, the bound forces the NHC-baro thermostat to take
    // over within a few steps and bring eps_dot back to physical values.
    const double eps_dt = std::clamp(eps_dot_ * dt_,
                                     -kEpsDtMaxFraction, kEpsDtMaxFraction);
    const double a       = std::exp(eps_dt);
    const double sinhc_h = sinhc(eps_dt * 0.5);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.x[i] = f.atoms.x[i] * a + f.atoms.vx[i] * dt_ * sinhc_h;
        f.atoms.y[i] = f.atoms.y[i] * a + f.atoms.vy[i] * dt_ * sinhc_h;
        f.atoms.z[i] = f.atoms.z[i] * a + f.atoms.vz[i] * dt_ * sinhc_h;
    }

    // Update box edges: V_new = V * exp(3 eps_dot dt) ⇒ each edge *= exp(eps_dot dt).
    f.box(0, 0) *= a;
    f.box(1, 1) *= a;
    f.box(2, 2) *= a;

    // Wrap (orthorhombic) owned atoms into the freshly-scaled box.
    wrap_orthorhombic(f, n);

    // The per-rank slab bounds (sublo/subhi) were set by Backend::set_frame
    // from the box at that point. NPT just changed the box -- if we don't
    // refresh the slab, the next AtomMigration / HaloExchange will see
    // owned-atom positions that no longer fit the stale boundary, yielding
    // out-of-range cell-list indices and a crash. compute_sublo_subhi is
    // O(1) (three multiplies) and is called identically on every rank, so
    // the partition stays consistent.
    if (ctx.runtime.dd && ctx.runtime.mpi_size > 1) {
        ctx.runtime.dd->compute_sublo_subhi(f.box);
    }
}

void NoseHooverNPT::launch(const cpu::LaunchContext& ctx) {
    const double half_dt = 0.5 * dt_;

    if (ctx.stage == static_cast<int>(ATV::md::Phase::initial_integrate)) {
        // First half of Trotter (P1 -- P6).
        nhc_atoms_half(ctx, half_dt);
        nhc_baro_half(half_dt);
        update_eps_dot(ctx, half_dt);   // uses prev step's reduced stress
        scale_v_by_baro(ctx, half_dt);
        kick_v_by_force(ctx, half_dt);
        propagate_x_and_box(ctx);

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::final_integrate)) {
        // Force has been recomputed; second half-kick (P7).
        kick_v_by_force(ctx, half_dt);

    } else if (ctx.stage == static_cast<int>(ATV::md::Phase::end_of_step)) {
        // ReduceEnergy has run (frame.stress globally reduced).
        // Recompute KE_global_ from CURRENT v before pressure step.
        KE_global_ = compute_global_KE(ctx);

        scale_v_by_baro(ctx, half_dt);   // P8
        update_eps_dot(ctx, half_dt);    // P9 -- new pressure
        nhc_baro_half(half_dt);          // P10
        nhc_atoms_half(ctx, half_dt);    // P11

        // Rewrite frame.kinetic so observers bound at end_of_step @ > -40
        // see the post-NPT kinetic energy.
        ctx.frame.kinetic = KE_global_;
    }
}

} // namespace md
} // namespace cpu
} // namespace ATV
