/// @file cpu/minimizer/lbfgs.cpp
/// @brief L-BFGS step implementation: two-loop recursion + in-step line
///        search + multi-criteria convergence. MPI-aware: owned-row loops +
///        cross-rank reductions.

#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/component/owned_range.h"
#include "atomiverse/cpu/minimizer/global_reduce.h"

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

LBFGS& LBFGS::line_search(LineSearchMode mode, double c1, double c2,
                          double max_step, double min_alpha) {
    ls_mode_ = mode;
    ls_c1_ = c1;
    ls_c2_ = c2;
    ls_max_step_ = max_step;
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
    ls_phi0_ = 0;
    ls_gd0_ = 0;
    ls_dir_max_ = 0;
    last_generation_ = ~std::uint64_t{0};
    trace_.clear();

    const char* mode_name = (ls_mode_ == LineSearchMode::StrongWolfe)
                                ? "strong_wolfe"
                                : (ls_mode_ == LineSearchMode::Armijo)
                                      ? "armijo"
                                      : "fixed";
    ATV_INFO(
        "[lbfgs] build memory={} alpha0={:.6f} ftol={:.6f} rms_ftol={:.6f} "
        "max_dtol={:.4f} rms_dtol={:.4f} etol={:.6e} line_search={}",
        memory_, alpha0_, ftol_, rms_ftol_, max_dtol_, rms_dtol_, etol_,
        mode_name);
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

    // Atom migration reshuffles local row indexing: curvature history and
    // previous-step state key off row order and must be discarded.
    // `neighbor_generation` is bumped in lockstep on every rank (Migrate gates
    // the exchange on a global allreduce), so this reset stays rank-uniform.
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
        const double energy_change = (step_ > 0)
                                         ? std::abs(f.scalars.energy
                                                    - energy_prev_)
                                         : double(1e30);

        OptRecord rec;
        rec.step = step_;
        rec.energy = f.scalars.energy;
        rec.max_force = max_force_;
        rec.rms_force = rms_force;
        rec.max_displacement = max_disp;
        rec.rms_displacement = rms_disp;
        rec.energy_change = energy_change;
        trace_.push_back(rec);

        converged_ = (max_force_ < ftol_) && (rms_force < rms_ftol_)
                     && (step_ == 0 || !disp_valid || max_disp < max_dtol_)
                     && (step_ == 0 || !disp_valid || rms_disp < rms_dtol_)
                     && (step_ == 0 || energy_change < etol_);

        ATV_INFO("[lbfgs] step={} energy_eV={:.6f} max_force_eVA={:.6f} "
                 "rms_force_eVA={:.6f} max_disp_A={:.6e} rms_disp_A={:.6e} "
                 "dE_eV={:.6e} converged={}",
                 step_, f.scalars.energy, max_force_, rms_force, max_disp,
                 rms_disp, energy_change, converged_);

        if (converged_) {
            ctx.runtime.outer_status.converged = true;
            ctx.runtime.outer_status.reason =
                "lbfgs multi-criteria convergence";
        }

        energy_prev_ = f.scalars.energy;
        return;
    }

    // ── Line-search phase ────────────────────────────────────────────
    //
    // search_direction (run first this step) has already stored `direction_`,
    // the descent slope `ls_gd0_ = g.direction`, the origin energy `ls_phi0_`,
    // and the global `ls_dir_max_`. Here we choose the step length and write
    // the accepted geometry x_prev_ + alpha*direction_. The driver's force
    // phase (next) recomputes forces/energy at that geometry for the
    // convergence check, so we do not re-evaluate after committing alpha.
    if (phase == static_cast<int>(ATV::minimizer::Phase::line_search)) {
        // Displacement cap: bound the trial step so a wild move cannot blow up
        // an expensive (SCF) evaluator. alpha is dimensionless on direction_.
        const double cap =
            (ls_dir_max_ > 1e-20) ? ls_max_step_ / ls_dir_max_ : 1.0;

        // Fixed mode (or no re-eval handle available): legacy single step.
        if (ls_mode_ == LineSearchMode::Fixed || !ctx.runtime.reeval_forces) {
            double alpha = alpha0_;
            if (ls_dir_max_ > 1e-20)
                alpha = std::min(alpha0_, 0.1 / ls_dir_max_); // legacy 0.1 Å cap
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                f.atoms.x(i) = x_prev_[i * 3 + 0] + alpha * direction_[i * 3 + 0];
                f.atoms.y(i) = x_prev_[i * 3 + 1] + alpha * direction_[i * 3 + 1];
                f.atoms.z(i) = x_prev_[i * 3 + 2] + alpha * direction_[i * 3 + 2];
            }
            return;
        }

        // No descent (step 0, or degenerate): leave the geometry unchanged;
        // the convergence check decides whether we are already at a minimum.
        if (!(ls_gd0_ < 0))
            return;

        // phi(alpha) writes the trial geometry, replays the force phase via
        // the driver's re-eval handle, and returns the globally-reduced
        // energy and directional derivative g(x+alpha p).p.
        auto phi = [&](double alpha, bool need_slope) -> LsSample {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                f.atoms.x(i) = x_prev_[i * 3 + 0] + alpha * direction_[i * 3 + 0];
                f.atoms.y(i) = x_prev_[i * 3 + 1] + alpha * direction_[i * 3 + 1];
                f.atoms.z(i) = x_prev_[i * 3 + 2] + alpha * direction_[i * 3 + 2];
            }
            ctx.runtime.reeval_forces(need_slope);
            LsSample s;
            s.f = f.scalars.energy; // globally reduced by md::mpi::Reduce
            double gd_local = 0;
#pragma omp parallel for schedule(static) reduction(+ : gd_local)
            for (int i = 0; i < n; ++i) {
                gd_local += -f.atoms.fx(i) * direction_[i * 3 + 0]
                            + -f.atoms.fy(i) * direction_[i * 3 + 1]
                            + -f.atoms.fz(i) * direction_[i * 3 + 2];
            }
            s.dphi = global_sum(ctx, gd_local);
            return s;
        };

        const double alpha_max = std::max(cap, ls_min_alpha_);
        const double alpha_init = std::min(1.0, alpha_max);
        LineSearchResult lr;
        // Suppress migration while trials replay the force phase (see
        // Context::line_search_trial); halo exchange still refreshes ghosts.
        ctx.runtime.line_search_trial = true;
        if (ls_mode_ == LineSearchMode::StrongWolfe)
            lr = strong_wolfe(phi, ls_phi0_, ls_gd0_, alpha_init, ls_c1_, ls_c2_,
                              max_ls_, alpha_max);
        else // Armijo
            lr = armijo_backtrack(phi, ls_phi0_, ls_gd0_, alpha_init, ls_c1_,
                                  ls_contract_, max_ls_, ls_min_alpha_);
        ctx.runtime.line_search_trial = false;

        // Commit the accepted geometry (the last trial may have been a
        // different alpha). On a failed search (e.g. machine-flat near the
        // minimum) lr.alpha may be ~0 → effectively no move, which the
        // convergence check then resolves.
        const double alpha = lr.alpha;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            f.atoms.x(i) = x_prev_[i * 3 + 0] + alpha * direction_[i * 3 + 0];
            f.atoms.y(i) = x_prev_[i * 3 + 1] + alpha * direction_[i * 3 + 1];
            f.atoms.z(i) = x_prev_[i * 3 + 2] + alpha * direction_[i * 3 + 2];
        }
        ATV_DEBUG("[lbfgs] line_search step={} mode={} alpha={:.6e} ok={} "
                  "n_eval={}",
                  step_, static_cast<int>(ls_mode_), alpha, lr.ok,
                  lr.n_value_evals);
        return;
    }

    // ── Search-direction phase ───────────────────────────────────────

    // Gradient = -force (we minimize energy, force = -grad E).
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

    ATV_DEBUG("[lbfgs] search_direction step={} stored_history={}", step_,
              n_stored_);

    // Update L-BFGS history with the previous step's accepted move. The
    // strong-Wolfe line search guarantees y^T s > 0 automatically (N&W §3.1);
    // we still assert it before committing the pair so a degenerate step
    // (Armijo/Fixed mode, or a search that bailed) cannot push a non-positive
    // curvature pair that would break the inverse-Hessian update.
    if (has_prev_) {
        const int idx = (n_stored_ < memory_) ? n_stored_ : oldest_;
        auto& s = s_history_[idx];
        auto& y = y_history_[idx];
        double ys_local = 0;
        double ss_local = 0;
        double yy_local = 0;
#pragma omp parallel for schedule(static)                                      \
    reduction(+ : ys_local, ss_local, yy_local)
        for (int j = 0; j < dim; ++j) {
            s[j] = x[j] - x_prev_[j];
            y[j] = g[j] - g_prev_[j];
            ys_local += y[j] * s[j];
            ss_local += s[j] * s[j];
            yy_local += y[j] * y[j];
        }
        double red[3] = {ys_local, ss_local, yy_local};
        global_sum(ctx, red, 3);
        const double ys = red[0];
        const double snorm = std::sqrt(red[1]);
        const double ynorm = std::sqrt(red[2]);

        // Curvature condition y^T s > eps*||s||*||y|| keeps H positive definite.
        if (ys > double(1e-12) * snorm * ynorm && ys > double(0)) {
            rho_history_[idx] = double(1) / ys;
            if (n_stored_ < memory_)
                ++n_stored_;
            else
                oldest_ = (oldest_ + 1) % memory_;
        } else {
            ATV_DEBUG("[lbfgs] skipped curvature pair step={} ys={:.3e} "
                      "(non-positive — keeps H PD)",
                      step_, ys);
        }
    }

    // Two-loop recursion → direction = -H g. Each inner product is a per-rank
    // partial reduced with global_sum so every rank shares the identical
    // direction for its owned rows.
    auto& q = q_ws_;
    std::copy(g.begin(), g.begin() + dim, q.begin());

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

    // direction = -H g (descent direction).
#pragma omp parallel for schedule(static)
    for (int j = 0; j < dim; ++j)
        direction_[j] = -r[j];

    // Save the pre-step state for the next history update.
    std::copy(x.begin(), x.begin() + dim, x_prev_.begin());
    std::copy(g.begin(), g.begin() + dim, g_prev_.begin());
    has_prev_ = true;
    x_prev_valid_ = true;

    // Descent precheck (N&W §3.1 / Alg PR+): the quasi-Newton direction is a
    // descent direction iff g.p < 0. A bad (s,y) history or numerical noise
    // can violate it; reset to steepest descent so the line search is
    // well-posed. The dot products are globally reduced for rank-uniformity.
    double gd0 = global_sum(ctx, dot_local(g, direction_, dim));
    if (!(gd0 < 0)) {
#pragma omp parallel for schedule(static)
        for (int j = 0; j < dim; ++j)
            direction_[j] = -g[j];
        gd0 = -global_sum(ctx, dot_local(g, g, dim));
        ATV_DEBUG("[lbfgs] non-descent direction step={} — reset to -g", step_);
    }

    // Global max|direction| component for the line-search displacement cap.
    double dir_max = 0;
#pragma omp parallel for schedule(static) reduction(max : dir_max)
    for (int j = 0; j < dim; ++j)
        dir_max = std::max(dir_max, std::abs(direction_[j]));
    ls_dir_max_ = global_max(ctx, dir_max);
    ls_gd0_ = gd0;
    ls_phi0_ = f.scalars.energy; // energy at x_prev_ (pre-step, this geometry)

    ++step_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
