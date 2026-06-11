/// @file cpu/minimizer/lbfgs.cpp
/// @brief L-BFGS step implementation: two-loop recursion + multi-criteria
///        convergence. MPI-aware: owned-row loops + cross-rank reductions.

#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/minimizer/global_reduce.h"
#include "atomiverse/cpu/component/owned_range.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace ATV {
namespace cpu {
namespace minimizer {

namespace {

/// Owned-row inner product (per-rank partial; reduce with `global_sum`).
double dot_local(const std::vector<double>& a, const std::vector<double>& b,
                 int dim) {
    double s = 0;
#pragma omp parallel for schedule(static) reduction(+ : s)
    for (int j = 0; j < dim; ++j)
        s += a[j] * b[j];
    return s;
}

} // namespace

LBFGS::LBFGS(int memory, double alpha0, double ftol, int max_ls,
             double rms_ftol, double max_dtol, double rms_dtol, double etol)
    : memory_(memory), alpha0_(alpha0), ftol_(ftol), max_ls_(max_ls),
      rms_ftol_(rms_ftol), max_dtol_(max_dtol), rms_dtol_(rms_dtol),
      etol_(etol) {}

LBFGS& LBFGS::line_search(double c1, double contract, double min_alpha) {
    line_search_enabled_ = true;
    ls_c1_ = c1;
    ls_contract_ = contract;
    ls_min_alpha_ = min_alpha;
    return *this;
}

void LBFGS::build(const cpu::BuildContext& ctx) {
    // Size every buffer to the frame's allocated capacity: under multi-rank
    // `atoms.size()` is the backing-store capacity (≥ n_local + n_ghost), so
    // the owned-row count can never outgrow it. All per-step work happens in
    // these build-time buffers — no allocation in launch().
    const int n = ctx.frame.atoms.size();
    const int dim = n * 3;

    s_history_.resize(memory_);
    y_history_.resize(memory_);
    rho_history_.assign(memory_, 0);
    for (int i = 0; i < memory_; ++i) {
        s_history_[i].assign(dim, 0);
        y_history_[i].assign(dim, 0);
    }

    x_prev_.assign(dim, 0);
    g_prev_.assign(dim, 0);
    direction_.assign(dim, 0);
    x_ws_.assign(dim, 0);
    g_ws_.assign(dim, 0);
    q_ws_.assign(dim, 0);
    r_ws_.assign(dim, 0);
    alpha_ws_.assign(memory_, 0);

    n_stored_ = 0;
    oldest_ = 0;
    has_prev_ = false;
    x_prev_valid_ = false;
    converged_ = false;
    energy_prev_ = 0;
    step_ = 0;
    ls_had_trial_ = false;
    ls_alpha_ = 0;
    ls_gd_ = 0;
    last_generation_ = ~std::uint64_t{0};
    trace_.clear();

    ATV_INFO(
        "[lbfgs] build memory={} alpha0={:.6f} ftol={:.6f} rms_ftol={:.6f} "
        "max_dtol={:.4f} rms_dtol={:.4f} etol={:.6e} line_search={}",
        memory_, alpha0_, ftol_, rms_ftol_, max_dtol_, rms_dtol_, etol_,
        line_search_enabled_ ? "armijo_backtrack" : "off");
}

void LBFGS::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const int dim = n * 3;

    // Capacity growth (driver-level rebuild_frame after a migration influx)
    // is the one legal re-allocation point.
    if (dim > static_cast<int>(x_ws_.size())) {
        for (auto& h : s_history_)
            h.resize(dim, 0);
        for (auto& h : y_history_)
            h.resize(dim, 0);
        x_prev_.resize(dim, 0);
        g_prev_.resize(dim, 0);
        direction_.resize(dim, 0);
        x_ws_.resize(dim, 0);
        g_ws_.resize(dim, 0);
        q_ws_.resize(dim, 0);
        r_ws_.resize(dim, 0);
    }

    // Atom migration reshuffles local row indexing: curvature history,
    // previous-step state, and the in-flight line-search trial all key off
    // row order and must be discarded. `neighbor_generation` is bumped in
    // lockstep on every rank (Migrate gates the exchange on a global
    // allreduce), so this reset — and everything downstream of it — stays
    // rank-uniform. Single-rank: generation never changes after the first
    // observation.
    {
        const auto* dd = ctx.runtime.dd;
        const std::uint64_t gen =
            (dd && ctx.runtime.mpi_size > 1) ? dd->neighbor_generation : 0;
        if (gen != last_generation_) {
            if (last_generation_ != ~std::uint64_t{0}) {
                n_stored_ = 0;
                oldest_ = 0;
                has_prev_ = false;
                x_prev_valid_ = false;
                ls_had_trial_ = false;
                ATV_DEBUG("[lbfgs] migration detected (gen {}) — history "
                          "reset, restarting as steepest descent",
                          gen);
            }
            last_generation_ = gen;
        }
    }

    const int phase = ctx.stage;

    // ── Convergence check phase ──────────────────────────────────────
    if (phase == static_cast<int>(ATV::minimizer::Phase::convergence_check)) {
        // Per-rank force partials over owned rows.
        double local_max_force = 0;
        double sum_f2 = 0;
#pragma omp parallel for schedule(static) reduction(+ : sum_f2)                \
    reduction(max : local_max_force)
        for (int i = 0; i < n; ++i) {
            const double fx = f.atoms.fx(i), fy = f.atoms.fy(i),
                         fz = f.atoms.fz(i);
            const double f2 = fx * fx + fy * fy + fz * fz;
            sum_f2 += f2;
            local_max_force = std::max(local_max_force, std::sqrt(f2));
        }

        // Displacement partials relative to pre-step positions (x_prev_ is
        // saved in search_direction before the step is applied). Invalid on
        // the first step and right after a migration reshuffle — the
        // criteria are then skipped exactly like step 0.
        const bool disp_valid = step_ > 0 && x_prev_valid_;
        double local_max_disp = 0;
        double sum_d2 = 0;
        if (disp_valid) {
#pragma omp parallel for schedule(static) reduction(+ : sum_d2)                \
    reduction(max : local_max_disp)
            for (int i = 0; i < n; ++i) {
                const double dx = f.atoms.x(i) - x_prev_[i * 3 + 0];
                const double dy = f.atoms.y(i) - x_prev_[i * 3 + 1];
                const double dz = f.atoms.z(i) - x_prev_[i * 3 + 2];
                sum_d2 += dx * dx + dy * dy + dz * dz;
                local_max_disp = std::max(
                    local_max_disp,
                    std::max({std::abs(dx), std::abs(dy), std::abs(dz)}));
            }
        }

        // Two batched collectives: SUM{Σf², Σd², n_owned} + MAX{f, d}.
        // Entered by every rank unconditionally (disp partials are zero
        // when invalid), so the collectives stay in lockstep.
        double sums[3] = {sum_f2, sum_d2, static_cast<double>(n)};
        global_sum(ctx, sums, 3);
        double maxs[2] = {local_max_force, local_max_disp};
        global_max(ctx, maxs, 2);

        const double n_global = std::max(sums[2], 1.0);
        max_force_ = maxs[0];
        const double rms_force = std::sqrt(sums[0] / (3 * n_global));
        const double max_disp = maxs[1];
        const double rms_disp =
            disp_valid ? std::sqrt(sums[1] / (3 * n_global)) : double(0);
        // frame.scalars.energy is the globally-reduced total (mpi::Reduce
        // runs at force @ +20, before this phase).
        const double energy_change = (step_ > 0)
                                         ? std::abs(f.scalars.energy
                                                    - energy_prev_)
                                         : double(1e30);

        // Record trace
        OptRecord rec;
        rec.step = step_;
        rec.energy = f.scalars.energy;
        rec.max_force = max_force_;
        rec.rms_force = rms_force;
        rec.max_displacement = max_disp;
        rec.rms_displacement = rms_disp;
        rec.energy_change = energy_change;
        trace_.push_back(rec);

        // All-criteria convergence (skip step 0 for displacement/energy
        // criteria; skip displacement criteria while x_prev_ is invalid).
        converged_ = (max_force_ < ftol_) && (rms_force < rms_ftol_)
                     && (step_ == 0 || !disp_valid || max_disp < max_dtol_)
                     && (step_ == 0 || !disp_valid || rms_disp < rms_dtol_)
                     && (step_ == 0 || energy_change < etol_);

        ATV_INFO("[lbfgs] step={} energy_eV={:.6f} max_force_eVA={:.6f} "
                 "rms_force_eVA={:.6f} max_disp_A={:.6e} rms_disp_A={:.6e} "
                 "dE_eV={:.6e} converged={}",
                 step_, f.scalars.energy, max_force_, rms_force, max_disp,
                 rms_disp, energy_change, converged_);

        // Report convergence to the outer OptimizeDriver loop so it can
        // break early. All ranks computed identical global metrics, so the
        // verdict — and the driver's early break — is rank-uniform.
        if (converged_) {
            ctx.runtime.outer_status.converged = true;
            ctx.runtime.outer_status.reason =
                "lbfgs multi-criteria convergence";
        }

        energy_prev_ = f.scalars.energy;
        return;
    }

    // ── Search direction phase ───────────────────────────────────────

    // Gradient = -force (we minimize energy, force = -grad E)
    auto& x = x_ws_;
    auto& g = g_ws_;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        x[i * 3 + 0] = f.atoms.x(i);
        x[i * 3 + 1] = f.atoms.y(i);
        x[i * 3 + 2] = f.atoms.z(i);
        g[i * 3 + 0] = -f.atoms.fx(i);
        g[i * 3 + 1] = -f.atoms.fy(i);
        g[i * 3 + 2] = -f.atoms.fz(i);
    }

    // ── Armijo line-search rejection & rollback ──────────────────────
    //
    // On step k > 0 we evaluate the trial taken during step k-1's
    // search_direction: x_k = x_{k-1} + alpha * d_{k-1}, E_k is now in
    // f.scalars.energy (globally reduced by mpi::Reduce). Armijo requires
    // E_k - E_{k-1} <= c1 * alpha * (g_{k-1} . d_{k-1}); `ls_gd_` was
    // globally summed when the trial was taken, so accept/reject is
    // rank-uniform.
    //
    // On rejection we roll back to x_{k-1}, scale alpha by `contract`,
    // re-take the step, and return without updating L-BFGS history. The
    // driver's next force phase will evaluate the shrunk trial, so the
    // full backtracking loop plays out over successive driver steps.
    if (line_search_enabled_ && ls_had_trial_) {
        const double armijo_rhs = ls_c1_ * ls_alpha_ * ls_gd_; // < 0 on descent
        const double dE = f.scalars.energy - ls_origin_energy_;
        const bool armijo_ok = (dE <= armijo_rhs);
        if (!armijo_ok && ls_alpha_ > ls_min_alpha_) {
            ls_alpha_ *= ls_contract_;
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                f.atoms.x(i) =
                    x_prev_[i * 3 + 0] + ls_alpha_ * direction_[i * 3 + 0];
                f.atoms.y(i) =
                    x_prev_[i * 3 + 1] + ls_alpha_ * direction_[i * 3 + 1];
                f.atoms.z(i) =
                    x_prev_[i * 3 + 2] + ls_alpha_ * direction_[i * 3 + 2];
            }
            ATV_DEBUG("[lbfgs] armijo reject step={} dE={:.3e} rhs={:.3e} "
                      "new_alpha={:.3e}",
                      step_, dE, armijo_rhs, ls_alpha_);
            return; // history update + new direction happen on the next accepted trial
        }
    }

    ATV_DEBUG("[lbfgs] search_direction step={} stored_history={}", step_,
              n_stored_);

    // Update L-BFGS history with previous step's data. The candidate (s, y)
    // pair is written straight into its history slot; the slot is only
    // committed (counters advanced) when the globally-summed curvature
    // condition holds, otherwise the next candidate overwrites it.
    if (has_prev_) {
        const int idx = (n_stored_ < memory_) ? n_stored_ : oldest_;
        auto& s = s_history_[idx];
        auto& y = y_history_[idx];
        double ys_local = 0;
#pragma omp parallel for schedule(static) reduction(+ : ys_local)
        for (int j = 0; j < dim; ++j) {
            s[j] = x[j] - x_prev_[j];
            y[j] = g[j] - g_prev_[j];
            ys_local += y[j] * s[j];
        }
        const double ys = global_sum(ctx, ys_local);

        // Only store if curvature condition is satisfied
        if (ys > double(1e-10)) {
            rho_history_[idx] = double(1) / ys;
            if (n_stored_ < memory_) {
                ++n_stored_;
            } else {
                oldest_ = (oldest_ + 1) % memory_;
            }
        }
    }

    // Two-loop recursion to compute search direction. Every inner product
    // is a per-rank partial reduced with global_sum so all ranks share the
    // identical direction for their owned rows.
    auto& q = q_ws_;
    std::copy(g.begin(), g.begin() + dim, q.begin());

    // Backward pass: from newest to oldest
    for (int k = n_stored_ - 1; k >= 0; --k) {
        const int idx = (oldest_ + k) % memory_;
        const double sq = global_sum(ctx, dot_local(s_history_[idx], q, dim));
        const double a = rho_history_[idx] * sq;
        alpha_ws_[k] = a;
        const auto& yk = y_history_[idx];
#pragma omp parallel for schedule(static)
        for (int j = 0; j < dim; ++j)
            q[j] -= a * yk[j];
    }

    // Initial Hessian approximation: H0 = gamma * I
    double gamma = 1;
    if (n_stored_ > 0) {
        const int newest = (oldest_ + n_stored_ - 1) % memory_;
        const double yy = global_sum(
            ctx, dot_local(y_history_[newest], y_history_[newest], dim));
        if (yy > double(1e-30))
            gamma = double(1) / (rho_history_[newest] * yy);
    }

    auto& r = r_ws_;
#pragma omp parallel for schedule(static)
    for (int j = 0; j < dim; ++j)
        r[j] = gamma * q[j];

    // Forward pass: from oldest to newest
    for (int k = 0; k < n_stored_; ++k) {
        const int idx = (oldest_ + k) % memory_;
        const double yr = global_sum(ctx, dot_local(y_history_[idx], r, dim));
        const double b = rho_history_[idx] * yr;
        const auto& sk = s_history_[idx];
        const double coeff = alpha_ws_[k] - b;
#pragma omp parallel for schedule(static)
        for (int j = 0; j < dim; ++j)
            r[j] += sk[j] * coeff;
    }

    // direction = -H * g (descent direction)
#pragma omp parallel for schedule(static)
    for (int j = 0; j < dim; ++j)
        direction_[j] = -r[j];

    // Save current state for next step's history update
    std::copy(x.begin(), x.begin() + dim, x_prev_.begin());
    std::copy(g.begin(), g.begin() + dim, g_prev_.begin());
    has_prev_ = true;
    x_prev_valid_ = true;

    // Apply step: x_new = x + alpha * direction
    double alpha = alpha0_;

    // Scale step by inverse max-force to keep steps reasonable. The cap
    // must be computed from the GLOBAL direction maximum — a per-rank cap
    // would make ranks step with different alphas and tear the geometry.
    double dir_max = 0;
#pragma omp parallel for schedule(static) reduction(max : dir_max)
    for (int j = 0; j < dim; ++j)
        dir_max = std::max(dir_max, std::abs(direction_[j]));
    dir_max = global_max(ctx, dir_max);
    if (dir_max > double(1e-20))
        alpha =
            std::min(alpha, double(0.1) / dir_max); // cap step at 0.1 Angstrom

    // Remember the trial so the next step's search_direction can test Armijo.
    if (line_search_enabled_) {
        ls_alpha_ = alpha;
        ls_gd_ = global_sum(ctx, dot_local(g, direction_, dim));
        ls_origin_energy_ = f.scalars.energy; // energy at x_prev_ (pre-step)
        ls_had_trial_ = true;
    }

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.x(i) += alpha * direction_[i * 3 + 0];
        f.atoms.y(i) += alpha * direction_[i * 3 + 1];
        f.atoms.z(i) += alpha * direction_[i * 3 + 2];
    }

    ATV_DEBUG("[lbfgs] applied step={} alpha={:.6e} dir_max={:.6e}", step_,
              alpha, dir_max);

    ++step_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
