/// @file cpu/minimizer/cg.cpp
/// @brief Polak-Ribière+ CG step + multi-criteria convergence (mirrors LBFGS).
///        MPI-aware: owned-row loops + cross-rank reductions.

#include "atomiverse/cpu/minimizer/cg.h"
#include "atomiverse/core/log.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/component/owned_range.h"
#include "atomiverse/cpu/minimizer/global_reduce.h"
#include "atomiverse/cpu/minimizer/line_search.h"

#include <algorithm>
#include <cmath>

namespace ATV {
namespace cpu {
namespace minimizer {

CG::CG(double alpha0, double ftol, double rms_ftol, double max_dtol,
       double rms_dtol, double etol)
    : alpha0_(alpha0), ftol_(ftol), rms_ftol_(rms_ftol), max_dtol_(max_dtol),
      rms_dtol_(rms_dtol), etol_(etol) {}

CG& CG::line_search(LineSearchMode mode, double c1, double c2, double max_step,
                    double min_alpha) {
    ls_mode_ = mode;
    ls_c1_ = c1;
    ls_c2_ = c2;
    ls_max_step_ = max_step;
    ls_min_alpha_ = min_alpha;
    return *this;
}

void CG::build(const cpu::BuildContext& ctx) {
    // Capacity-sized buffers (atoms.size() ≥ n_local under multi-rank) —
    // all per-step work happens here, no allocation in launch().
    const int dim = ctx.frame.atoms.size() * 3;
    x_prev_.assign(dim, 0.0);
    g_prev_.assign(dim, 0.0);
    d_prev_.assign(dim, 0.0);
    x_ws_.assign(dim, 0.0);
    g_ws_.assign(dim, 0.0);
    d_ws_.assign(dim, 0.0);
    has_prev_ = false;
    x_prev_valid_ = false;
    converged_ = false;
    step_ = 0;
    last_generation_ = ~std::uint64_t{0};
    trace_.clear();
}

void CG::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = owned_count(ctx);
    const int dim = n * 3;
    const int phase = ctx.stage;

    // Capacity growth is the one legal re-allocation point.
    if (dim > static_cast<int>(x_ws_.size())) {
        x_prev_.resize(dim, 0.0);
        g_prev_.resize(dim, 0.0);
        d_prev_.resize(dim, 0.0);
        x_ws_.resize(dim, 0.0);
        g_ws_.resize(dim, 0.0);
        d_ws_.resize(dim, 0.0);
    }

    // Migration reshuffles local rows → previous-step vectors are keyed to
    // dead indexing; restart as steepest descent. Rank-uniform: the
    // generation bump is gated on a global allreduce inside Migrate.
    {
        const auto* dd = ctx.runtime.md.dd;
        const std::uint64_t gen =
            (dd && ctx.runtime.mpi_size > 1) ? dd->neighbor_generation : 0;
        if (gen != last_generation_) {
            if (last_generation_ != ~std::uint64_t{0}) {
                has_prev_ = false;
                x_prev_valid_ = false;
                ATV_DEBUG("[cg] migration detected (gen {}) — restarting as "
                          "steepest descent",
                          gen);
            }
            last_generation_ = gen;
        }
    }

    // ── Convergence check ──────────────────────────────────────────────
    if (phase == static_cast<int>(ATV::minimizer::Phase::convergence_check)) {
        double local_max_force = 0;
        double sum_f2 = 0;
#pragma omp parallel for schedule(static) reduction(+ : sum_f2)                \
    reduction(max : local_max_force)
        for (int i = 0; i < n; ++i) {
            const double fx = f.atoms.fx(i);
            const double fy = f.atoms.fy(i);
            const double fz = f.atoms.fz(i);
            const double f2 = fx * fx + fy * fy + fz * fz;
            sum_f2 += f2;
            local_max_force = std::max(local_max_force, std::sqrt(f2));
        }

        const bool disp_valid = step_ > 0 && x_prev_valid_;
        double local_max_disp = 0, sum_d2 = 0;
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

        // Batched collectives, entered unconditionally by every rank.
        double sums[3] = {sum_f2, sum_d2, static_cast<double>(n)};
        global_sum(ctx, sums, 3);
        double maxs[2] = {local_max_force, local_max_disp};
        global_max(ctx, maxs, 2);

        const double n_global = std::max(sums[2], 1.0);
        max_force_ = maxs[0];
        const double rms_force = std::sqrt(sums[0] / (3 * n_global));
        const double max_disp = maxs[1];
        const double rms_disp =
            disp_valid ? std::sqrt(sums[1] / (3 * n_global)) : 0.0;
        const double energy_change =
            (step_ > 0) ? std::abs(f.scalars.energy - energy_prev_) : 1e30;

        OptRecord rec{};
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

        ATV_INFO("[cg] step={} energy_eV={:.6f} max_force_eVA={:.6f} "
                 "rms_force_eVA={:.6f} max_disp_A={:.6e} rms_disp_A={:.6e} "
                 "dE_eV={:.6e} converged={}",
                 step_, f.scalars.energy, max_force_, rms_force, max_disp,
                 rms_disp, energy_change, converged_);

        if (converged_) {
            ctx.runtime.outer_status.converged = true;
            ctx.runtime.outer_status.reason = "cg multi-criteria convergence";
        }

        energy_prev_ = f.scalars.energy;
        return;
    }

    // ── Line-search phase ───────────────────────────────────────────────
    //
    // search_direction (this step) stored the PR+ direction in d_prev_, the
    // origin geometry in x_prev_, and the descent slope/origin energy. Here we
    // pick the step length and write x_prev_ + alpha*d_prev_; the driver's
    // force phase (next) recomputes forces/energy at that geometry.
    if (phase == static_cast<int>(ATV::minimizer::Phase::line_search)) {
        const double cap =
            (ls_dir_max_ > 1e-20) ? ls_max_step_ / ls_dir_max_ : 1.0;

        // Fixed mode (or no re-eval handle): legacy single step x + alpha0*d.
        if (ls_mode_ == LineSearchMode::Fixed
            || !ctx.runtime.opt.reeval_forces) {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                f.atoms.x(i) =
                    x_prev_[i * 3 + 0] + alpha0_ * d_prev_[i * 3 + 0];
                f.atoms.y(i) =
                    x_prev_[i * 3 + 1] + alpha0_ * d_prev_[i * 3 + 1];
                f.atoms.z(i) =
                    x_prev_[i * 3 + 2] + alpha0_ * d_prev_[i * 3 + 2];
            }
            return;
        }
        if (!(ls_gd0_ < 0))
            return; // step 0 / degenerate — leave geometry to the conv. check

        auto phi = [&](double alpha, bool need_slope) -> LsSample {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < n; ++i) {
                f.atoms.x(i) = x_prev_[i * 3 + 0] + alpha * d_prev_[i * 3 + 0];
                f.atoms.y(i) = x_prev_[i * 3 + 1] + alpha * d_prev_[i * 3 + 1];
                f.atoms.z(i) = x_prev_[i * 3 + 2] + alpha * d_prev_[i * 3 + 2];
            }
            ctx.runtime.opt.reeval_forces(need_slope);
            LsSample s;
            s.f = f.scalars.energy;
            double gd_local = 0;
#pragma omp parallel for schedule(static) reduction(+ : gd_local)
            for (int i = 0; i < n; ++i) {
                gd_local += -f.atoms.fx(i) * d_prev_[i * 3 + 0]
                            + -f.atoms.fy(i) * d_prev_[i * 3 + 1]
                            + -f.atoms.fz(i) * d_prev_[i * 3 + 2];
            }
            s.dphi = global_sum(ctx, gd_local);
            return s;
        };

        const double alpha_max = std::max(cap, ls_min_alpha_);
        const double alpha_init = std::min(1.0, alpha_max);
        LineSearchResult lr;
        ctx.runtime.opt.line_search_trial =
            true; // suppress migration during trials
        if (ls_mode_ == LineSearchMode::StrongWolfe)
            lr =
                strong_wolfe(phi, ls_phi0_, ls_gd0_, alpha_init, ls_c1_, ls_c2_,
                             /*max_evals=*/20, alpha_max);
        else
            lr =
                armijo_backtrack(phi, ls_phi0_, ls_gd0_, alpha_init, ls_c1_,
                                 ls_contract_, /*max_evals=*/20, ls_min_alpha_);
        ctx.runtime.opt.line_search_trial = false;

        const double alpha = lr.alpha;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            f.atoms.x(i) = x_prev_[i * 3 + 0] + alpha * d_prev_[i * 3 + 0];
            f.atoms.y(i) = x_prev_[i * 3 + 1] + alpha * d_prev_[i * 3 + 1];
            f.atoms.z(i) = x_prev_[i * 3 + 2] + alpha * d_prev_[i * 3 + 2];
        }
        return;
    }

    // ── Search direction (Polak-Ribière+) ──────────────────────────────
    auto& g = g_ws_;
    auto& x = x_ws_;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        x[i * 3 + 0] = f.atoms.x(i);
        x[i * 3 + 1] = f.atoms.y(i);
        x[i * 3 + 2] = f.atoms.z(i);
        g[i * 3 + 0] = -f.atoms.fx(i);
        g[i * 3 + 1] = -f.atoms.fy(i);
        g[i * 3 + 2] = -f.atoms.fz(i);
    }

    auto& d = d_ws_;
    if (!has_prev_) {
        // Step 0: steepest descent.
#pragma omp parallel for schedule(static)
        for (int k = 0; k < dim; ++k)
            d[k] = -g[k];
    } else {
        // β = max(0, (g · (g - g_prev)) / (g_prev · g_prev)) — both inner
        // products are global sums so every rank applies the identical β.
        double num = 0, den = 0;
#pragma omp parallel for schedule(static) reduction(+ : num, den)
        for (int k = 0; k < dim; ++k) {
            num += g[k] * (g[k] - g_prev_[k]);
            den += g_prev_[k] * g_prev_[k];
        }
        double nd[2] = {num, den};
        global_sum(ctx, nd, 2);
        const double beta = (nd[1] > 0.0) ? std::max(0.0, nd[0] / nd[1]) : 0.0;
#pragma omp parallel for schedule(static)
        for (int k = 0; k < dim; ++k)
            d[k] = -g[k] + beta * d_prev_[k];
    }

    // Descent precheck: PR+ does not always yield a descent direction (N&W
    // §5.2). Reset to steepest descent if g·d ≥ 0. Globally reduced.
    double gd_local = 0;
#pragma omp parallel for schedule(static) reduction(+ : gd_local)
    for (int k = 0; k < dim; ++k)
        gd_local += g[k] * d[k];
    double gd0 = global_sum(ctx, gd_local);
    if (!(gd0 < 0)) {
#pragma omp parallel for schedule(static)
        for (int k = 0; k < dim; ++k)
            d[k] = -g[k];
        double gg_local = 0;
#pragma omp parallel for schedule(static) reduction(+ : gg_local)
        for (int k = 0; k < dim; ++k)
            gg_local += g[k] * g[k];
        gd0 = -global_sum(ctx, gg_local);
    }

    // Global max|d| for the displacement cap.
    double dir_max = 0;
#pragma omp parallel for schedule(static) reduction(max : dir_max)
    for (int k = 0; k < dim; ++k)
        dir_max = std::max(dir_max, std::abs(d[k]));
    ls_dir_max_ = global_max(ctx, dir_max);
    ls_gd0_ = gd0;
    ls_phi0_ = f.scalars.energy;

    // Save state for next iteration (buffer swaps — no copies, no allocs).
    // After the swap the current direction is in d_prev_ and the origin
    // geometry in x_prev_, which the line_search phase reads.
    std::swap(x_prev_, x_ws_);
    std::swap(g_prev_, g_ws_);
    std::swap(d_prev_, d_ws_);
    has_prev_ = true;
    x_prev_valid_ = true;
    ++step_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
