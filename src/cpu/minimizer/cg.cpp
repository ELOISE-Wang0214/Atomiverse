/// @file cpu/minimizer/cg.cpp
/// @brief Polak-Ribière+ CG step + multi-criteria convergence (mirrors LBFGS).
///        MPI-aware: owned-row loops + cross-rank reductions.

#include "atomiverse/cpu/minimizer/cg.h"
#include "atomiverse/core/log.h"
#include "atomiverse/cpu/comm/domain.h"
#include "atomiverse/cpu/minimizer/global_reduce.h"
#include "atomiverse/cpu/owned_range.h"

#include <algorithm>
#include <cmath>

namespace ATV {
namespace cpu {
namespace minimizer {

CG::CG(double alpha0, double ftol, double rms_ftol, double max_dtol,
       double rms_dtol, double etol)
    : alpha0_(alpha0), ftol_(ftol), rms_ftol_(rms_ftol), max_dtol_(max_dtol),
      rms_dtol_(rms_dtol), etol_(etol) {}

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
        const auto* dd = ctx.runtime.dd;
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

    // Take the step: x_{k+1} = x + α · d
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        f.atoms.x(i) = x[i * 3 + 0] + alpha0_ * d[i * 3 + 0];
        f.atoms.y(i) = x[i * 3 + 1] + alpha0_ * d[i * 3 + 1];
        f.atoms.z(i) = x[i * 3 + 2] + alpha0_ * d[i * 3 + 2];
    }

    // Save state for next iteration (buffer swaps — no copies, no allocs).
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
