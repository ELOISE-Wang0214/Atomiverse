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
// 0.5% per half-step ≈ 1% per step.  The clamp is only meant to fence
// off the first few ill-conditioned equilibration steps of an
// over-compressed start state; a well-initialised run sees `eps_dot·dt`
// bounded by the barostat natural frequency `2π/tau_P ≈ 2×10⁻³` at
// `tau_P=1000 fs`, comfortably below this fence.
constexpr double kEpsDtMaxFraction = 0.005;

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
      P_target_(cfg.P_target), tau_P_(cfg.tau_P), dt_(cfg.dt),
      drag_(cfg.drag), mtchain_(cfg.mtchain), mpchain_(cfg.mpchain) {}

void NoseHooverNPT::build(const cpu::BuildContext& ctx) {
    // atoms.size() at build time = global N (set by user before
    // Backend::set_frame partitions).
    const int N = static_cast<int>(ctx.frame.atoms.size());
    n_dof_ = 3 * N - 3;
    kT_    = kBoltzmann * T_target_;

    // Standard MTK masses (Tuckerman 2006 §III).
    //   Q_t[0] = N_f kT tau_T^2 drives the full KE reservoir;
    //   Q_t[k] = kT tau_T^2  for k >= 1  thermostats the previous link.
    Q_t_.assign(static_cast<std::size_t>(mtchain_), 0.0);
    Q_b_.assign(static_cast<std::size_t>(mpchain_), 0.0);
    eta_t_.assign(static_cast<std::size_t>(mtchain_), 0.0);
    eta_t_dot_.assign(static_cast<std::size_t>(mtchain_), 0.0);
    eta_b_.assign(static_cast<std::size_t>(mpchain_), 0.0);
    eta_b_dot_.assign(static_cast<std::size_t>(mpchain_), 0.0);

    const double tau_T2 = tau_T_ * tau_T_;
    const double tau_P2 = tau_P_ * tau_P_;
    Q_t_[0] = static_cast<double>(n_dof_) * kT_ * tau_T2;
    for (int k = 1; k < mtchain_; ++k) Q_t_[static_cast<std::size_t>(k)] = kT_ * tau_T2;
    Q_b_[0] = kT_ * tau_P2;                    // baro-NHC driven by eps_dot KE
    for (int k = 1; k < mpchain_; ++k) Q_b_[static_cast<std::size_t>(k)] = kT_ * tau_P2;

    // Barostat mass: MTK-consistent `(N_f + d) kT τ_P²` (Tuckerman 2006
    // Eq. 5.6.22). Numerically equivalent to LAMMPS `(N+1) kT τ_P²`
    // once our `3V` G_ε factor is balanced against LAMMPS's `V` factor
    // (both give `d ε̇ /dt ≈ V P / (N kT τ_P²)` in the large-N limit).
    W_ = (static_cast<double>(n_dof_) + 3.0) * kT_ * tau_P2;

    // LAMMPS-style empirical drag on *every* chain link of the backward
    // sweep. Matches `fix_nh.cpp::nhc_temp_integrate` (line ~2063–2076):
    //
    //     tdrag_factor = 1 - drag * update->dt * t_freq / nc_tchain
    //                  = 1 - drag * dt / tau_T   (with nc_tchain=1)
    //
    // and is multiplied into each `eta_dot[ich]` for `ich = mtchain-1 .. 0`
    // after the accelerated-then-exp-scaled update. An early Atomiverse
    // version applied `tdrag_factor` only to the innermost link with
    // `half_dt`, which was ~6x weaker damping (fewer links × half the
    // magnitude) and let chain oscillations ring on stiff-bond systems
    // (amorphous PE Stage 2) until atoms overlapped. LAMMPS's per-link
    // drag kills that mode.
    tdrag_factor_ = 1.0 - drag_ * dt_ / tau_T_;
    pdrag_factor_ = 1.0 - drag_ * dt_ / tau_P_;

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
    // M-chain NHC half-step, mirroring LAMMPS `fix_nh.cpp::nhc_temp_integrate`
    // (release branch). Per-link structure inside the backward sweep:
    //
    //     eta_dot[k] *= expfac           // pre-scale by outer link
    //     eta_dot[k] += G_k * dt/4       // accelerate
    //     eta_dot[k] *= tdrag_factor     // *** drag at every link ***
    //     eta_dot[k] *= expfac           // post-scale by outer link
    //
    // with `expfac = exp(-dt/8 * eta_dot[k+1])`; for the outermost link
    // (k = M-1) the "next" link `eta_dot[M]` is a zero sentinel, so
    // `expfac = 1` and the two scalings are no-ops — drag and G still
    // apply. The forward sweep re-runs the same structure with the
    // post-scale KE, minus drag (LAMMPS does not damp on the way out).
    const double quarter_dt = 0.5 * half_dt;
    const double eighth_dt  = 0.5 * quarter_dt;
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const int M = mtchain_;

    auto G_t = [&](int k) -> double {
        if (k == 0)
            return (2.0 * KE_global_ - n_dof_ * kT_) / Q_t_[0];
        return (Q_t_[static_cast<std::size_t>(k - 1)]
                  * eta_t_dot_[static_cast<std::size_t>(k - 1)]
                  * eta_t_dot_[static_cast<std::size_t>(k - 1)] - kT_)
                / Q_t_[static_cast<std::size_t>(k)];
    };

    // Backward pass: M-1 → 0. Drag applied at every link.
    for (int k = M - 1; k >= 0; --k) {
        const double expfac = (k < M - 1)
            ? std::exp(-eta_t_dot_[static_cast<std::size_t>(k + 1)] * eighth_dt)
            : 1.0;
        eta_t_dot_[static_cast<std::size_t>(k)] *= expfac;
        eta_t_dot_[static_cast<std::size_t>(k)] += G_t(k) * quarter_dt;
        eta_t_dot_[static_cast<std::size_t>(k)] *= tdrag_factor_;
        eta_t_dot_[static_cast<std::size_t>(k)] *= expfac;
    }

    // Scale atom velocities by the innermost thermostat momentum.
    const double scale = std::exp(-eta_t_dot_[0] * half_dt);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.vx[i] *= scale;
        f.atoms.vy[i] *= scale;
        f.atoms.vz[i] *= scale;
    }
    KE_global_ *= scale * scale;

    // Propagate the chain positions.
    for (int k = 0; k < M; ++k)
        eta_t_[static_cast<std::size_t>(k)] +=
            eta_t_dot_[static_cast<std::size_t>(k)] * half_dt;

    // Forward pass: 0 → M-1. No drag (LAMMPS pattern).
    for (int k = 0; k < M; ++k) {
        const double expfac = (k < M - 1)
            ? std::exp(-eta_t_dot_[static_cast<std::size_t>(k + 1)] * eighth_dt)
            : 1.0;
        eta_t_dot_[static_cast<std::size_t>(k)] *= expfac;
        eta_t_dot_[static_cast<std::size_t>(k)] += G_t(k) * quarter_dt;
        eta_t_dot_[static_cast<std::size_t>(k)] *= expfac;
    }
}

void NoseHooverNPT::nhc_baro_half(double half_dt) {
    // Same M-chain pattern as `nhc_atoms_half`, but the innermost link
    // thermostats the barostat momentum `eps_dot` instead of atom KE.
    // Drag (`pdrag_factor`) is applied at every link in the backward
    // sweep and dropped in the forward sweep, matching LAMMPS
    // `fix_nh.cpp::nhc_press_integrate`.
    const double quarter_dt = 0.5 * half_dt;
    const double eighth_dt  = 0.5 * quarter_dt;
    const int M = mpchain_;

    auto G_b = [&](int k) -> double {
        if (k == 0) return (W_ * eps_dot_ * eps_dot_ - kT_) / Q_b_[0];
        return (Q_b_[static_cast<std::size_t>(k - 1)]
                  * eta_b_dot_[static_cast<std::size_t>(k - 1)]
                  * eta_b_dot_[static_cast<std::size_t>(k - 1)] - kT_)
                / Q_b_[static_cast<std::size_t>(k)];
    };

    // Backward pass with per-link drag.
    for (int k = M - 1; k >= 0; --k) {
        const double expfac = (k < M - 1)
            ? std::exp(-eta_b_dot_[static_cast<std::size_t>(k + 1)] * eighth_dt)
            : 1.0;
        eta_b_dot_[static_cast<std::size_t>(k)] *= expfac;
        eta_b_dot_[static_cast<std::size_t>(k)] += G_b(k) * quarter_dt;
        eta_b_dot_[static_cast<std::size_t>(k)] *= pdrag_factor_;
        eta_b_dot_[static_cast<std::size_t>(k)] *= expfac;
    }

    // Scale the box-volume momentum by the innermost baro-thermostat.
    eps_dot_ *= std::exp(-eta_b_dot_[0] * half_dt);

    for (int k = 0; k < M; ++k)
        eta_b_[static_cast<std::size_t>(k)] +=
            eta_b_dot_[static_cast<std::size_t>(k)] * half_dt;

    // Forward pass, no drag.
    for (int k = 0; k < M; ++k) {
        const double expfac = (k < M - 1)
            ? std::exp(-eta_b_dot_[static_cast<std::size_t>(k + 1)] * eighth_dt)
            : 1.0;
        eta_b_dot_[static_cast<std::size_t>(k)] *= expfac;
        eta_b_dot_[static_cast<std::size_t>(k)] += G_b(k) * quarter_dt;
        eta_b_dot_[static_cast<std::size_t>(k)] *= expfac;
    }
}

void NoseHooverNPT::update_eps_dot(const cpu::LaunchContext& ctx, double half_dt) {
    const auto& f = ctx.frame;
    const double V = f.box.volume();
    if (V <= 0.0) return;

    // Standard MTK sign convention (Tuckerman 2006, LAMMPS `fix_nh.cpp`):
    //
    //     P_int = trace(σ) / (3V)        ; positive ⇔ compressive
    //     G_ε   = [3V (P_int − P_ext) + (6/N_f) KE_atoms] / W
    //
    // `frame.stress[0..2]` holds the LAMMPS-style pressure-tensor sum
    // (kinetic Σ m·v² plus pair virial Σ r·F), so trace/(3V) is the
    // standard internal pressure, no negation. Both the pressure term
    // and the KE coupling term then live in the same sign convention,
    // so a system at equilibrium has G_ε ≈ small and the box settles.
    //
    // Historical pitfall (fixed 2026-05-19): the previous formula used
    // `P_int = -trace/(3V)` plus `+(P_target - P_int)` to match the
    // Berendsen barostat's internal convention — but the KE coupling
    // was still in standard convention, so the two terms didn't cancel
    // at equilibrium and the box monotonically inflated even after T
    // equilibrated (caught by tests/cpu/md/test_lammps_npt_parity.cpp).
    const double P_int =
        (f.stress[0] + f.stress[1] + f.stress[2]) / (3.0 * V);

    const double G_eps = (3.0 * V * (P_int - P_target_)
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
    f.box.scale_isotropic(a);

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
