/// @file cpu/minimizer/lbfgs.cpp
/// @brief L-BFGS step implementation: two-loop recursion + multi-criteria convergence.

#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace ATV {
namespace cpu {
namespace minimizer {

LBFGS::LBFGS(int memory, double alpha0, double ftol, int max_ls,
                      double rms_ftol, double max_dtol, double rms_dtol, double etol)
    : memory_(memory), alpha0_(alpha0), ftol_(ftol), max_ls_(max_ls),
      rms_ftol_(rms_ftol), max_dtol_(max_dtol), rms_dtol_(rms_dtol), etol_(etol) {}

LBFGS& LBFGS::line_search(double c1, double contract, double min_alpha) {
    line_search_enabled_ = true;
    ls_c1_ = c1;
    ls_contract_ = contract;
    ls_min_alpha_ = min_alpha;
    return *this;
}

void LBFGS::build(const cpu::BuildContext& ctx) {
    int n = ctx.frame.atoms.size();
    int dim = n * 3;

    s_history_.resize(memory_);
    y_history_.resize(memory_);
    rho_history_.resize(memory_, 0);
    for (int i = 0; i < memory_; ++i) {
        s_history_[i].resize(dim, 0);
        y_history_[i].resize(dim, 0);
    }

    x_prev_.resize(dim, 0);
    g_prev_.resize(dim, 0);
    direction_.resize(dim, 0);

    n_stored_ = 0;
    oldest_ = 0;
    has_prev_ = false;
    converged_ = false;
    energy_prev_ = 0;
    step_ = 0;
    ls_had_trial_ = false;
    ls_alpha_ = 0;
    ls_gd_ = 0;
    trace_.clear();

    ATV_INFO(
        "[lbfgs] build memory={} alpha0={:.6f} ftol={:.6f} rms_ftol={:.6f} max_dtol={:.4f} rms_dtol={:.4f} etol={:.6e} line_search={}",
        memory_, alpha0_, ftol_, rms_ftol_, max_dtol_, rms_dtol_, etol_,
        line_search_enabled_ ? "armijo_backtrack" : "off");
}

void LBFGS::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    int n = f.atoms.size();
    int dim = n * 3;

    int phase = ctx.stage;

    // ── Convergence check phase ──────────────────────────────────────
    if (phase == static_cast<int>(ATV::minimizer::Phase::convergence_check)) {
        // Compute force metrics
        max_force_ = 0;
        double sum_f2 = 0;
        for (int i = 0; i < n; ++i) {
            double fx = f.atoms.fx(i), fy = f.atoms.fy(i), fz = f.atoms.fz(i);
            double mag = std::sqrt(fx*fx + fy*fy + fz*fz);
            if (mag > max_force_) max_force_ = mag;
            sum_f2 += fx*fx + fy*fy + fz*fz;
        }
        double rms_force = std::sqrt(sum_f2 / (3 * n));

        // Compute displacement metrics relative to pre-step positions (x_prev_
        // is saved in search_direction phase before the step is applied).
        double max_disp = 0;
        double sum_d2 = 0;
        if (step_ > 0) {
            for (int i = 0; i < n; ++i) {
                double dx = f.atoms.x(i) - x_prev_[i*3+0];
                double dy = f.atoms.y(i) - x_prev_[i*3+1];
                double dz = f.atoms.z(i) - x_prev_[i*3+2];
                sum_d2 += dx*dx + dy*dy + dz*dz;
                max_disp = std::max(max_disp, std::max({std::abs(dx), std::abs(dy), std::abs(dz)}));
            }
        }
        double rms_disp = (step_ > 0) ? std::sqrt(sum_d2 / (3 * n)) : double(0);
        double energy_change = (step_ > 0) ? std::abs(f.energy - energy_prev_) : double(1e30);

        // Record trace
        OptRecord rec;
        rec.step = step_;
        rec.energy = f.energy;
        rec.max_force = max_force_;
        rec.rms_force = rms_force;
        rec.max_displacement = max_disp;
        rec.rms_displacement = rms_disp;
        rec.energy_change = energy_change;
        trace_.push_back(rec);

        // All-criteria convergence (skip step 0 for displacement/energy criteria)
        converged_ = (max_force_ < ftol_)
                   && (rms_force < rms_ftol_)
                   && (step_ == 0 || max_disp < max_dtol_)
                   && (step_ == 0 || rms_disp < rms_dtol_)
                   && (step_ == 0 || energy_change < etol_);

        ATV_INFO(
            "[lbfgs] step={} energy_eV={:.6f} max_force_eVA={:.6f} rms_force_eVA={:.6f} max_disp_A={:.6e} rms_disp_A={:.6e} dE_eV={:.6e} converged={}",
            step_, f.energy, max_force_, rms_force, max_disp, rms_disp,
            energy_change, converged_);

        // Report convergence to the outer OptimizeDriver loop so it can
        // break early. Driver-level early-stop is the contract introduced
        // with `cpu::OptimizeDriver` in the driver-taxonomy redesign.
        if (converged_) {
            ctx.runtime.outer_status.converged = true;
            ctx.runtime.outer_status.reason = "lbfgs multi-criteria convergence";
        }

        energy_prev_ = f.energy;
        return;
    }

    // ── Search direction phase ───────────────────────────────────────

    // Gradient = -force (we minimize energy, force = -grad E)
    std::vector<double> g(dim);
    std::vector<double> x(dim);
    for (int i = 0; i < n; ++i) {
        x[i*3+0] = f.atoms.x(i);
        x[i*3+1] = f.atoms.y(i);
        x[i*3+2] = f.atoms.z(i);
        g[i*3+0] = -f.atoms.fx(i);
        g[i*3+1] = -f.atoms.fy(i);
        g[i*3+2] = -f.atoms.fz(i);
    }

    // ── Armijo line-search rejection & rollback ──────────────────────
    //
    // On step k > 0 we evaluate the trial taken during step k-1's
    // search_direction: x_k = x_{k-1} + alpha * d_{k-1}, E_k is now in
    // f.energy. Armijo requires E_k - E_{k-1} <= c1 * alpha * (g_{k-1} . d_{k-1}).
    //
    // On rejection we roll back to x_{k-1}, scale alpha by `contract`,
    // re-take the step, and return without updating L-BFGS history. The
    // driver's next force phase will evaluate the shrunk trial, so the
    // full backtracking loop plays out over successive driver steps.
    if (line_search_enabled_ && ls_had_trial_) {
        const double armijo_rhs = ls_c1_ * ls_alpha_ * ls_gd_;  // < 0 on descent
        const double dE = f.energy - ls_origin_energy_;
        const bool armijo_ok = (dE <= armijo_rhs);
        if (!armijo_ok && ls_alpha_ > ls_min_alpha_) {
            ls_alpha_ *= ls_contract_;
            for (int i = 0; i < n; ++i) {
                f.atoms.x(i) = x_prev_[i*3+0] + ls_alpha_ * direction_[i*3+0];
                f.atoms.y(i) = x_prev_[i*3+1] + ls_alpha_ * direction_[i*3+1];
                f.atoms.z(i) = x_prev_[i*3+2] + ls_alpha_ * direction_[i*3+2];
            }
            ATV_DEBUG(
                "[lbfgs] armijo reject step={} dE={:.3e} rhs={:.3e} new_alpha={:.3e}",
                step_, dE, armijo_rhs, ls_alpha_);
            return;  // history update + new direction happen on the next accepted trial
        }
    }

    ATV_DEBUG(
        "[lbfgs] search_direction step={} stored_history={}",
        step_, n_stored_);

    // Update L-BFGS history with previous step's data
    if (has_prev_) {
        std::vector<double> s(dim), y(dim);
        for (int j = 0; j < dim; ++j) {
            s[j] = x[j] - x_prev_[j];
            y[j] = g[j] - g_prev_[j];
        }
        double ys = 0;
        for (int j = 0; j < dim; ++j) ys += y[j] * s[j];

        // Only store if curvature condition is satisfied
        if (ys > double(1e-10)) {
            int idx = (n_stored_ < memory_) ? n_stored_ : oldest_;
            s_history_[idx] = s;
            y_history_[idx] = y;
            rho_history_[idx] = double(1) / ys;
            if (n_stored_ < memory_) {
                ++n_stored_;
            } else {
                oldest_ = (oldest_ + 1) % memory_;
            }
        }
    }

    // Two-loop recursion to compute search direction
    std::vector<double> q(g);
    std::vector<double> alpha_buf(n_stored_);

    // Backward pass: from newest to oldest
    for (int k = n_stored_ - 1; k >= 0; --k) {
        int idx = (oldest_ + k) % memory_;
        double a = 0;
        for (int j = 0; j < dim; ++j)
            a += rho_history_[idx] * s_history_[idx][j] * q[j];
        alpha_buf[k] = a;
        for (int j = 0; j < dim; ++j)
            q[j] -= a * y_history_[idx][j];
    }

    // Initial Hessian approximation: H0 = gamma * I
    double gamma = 1;
    if (n_stored_ > 0) {
        int newest = (oldest_ + n_stored_ - 1) % memory_;
        double yy = 0;
        for (int j = 0; j < dim; ++j)
            yy += y_history_[newest][j] * y_history_[newest][j];
        if (yy > double(1e-30))
            gamma = double(1) / (rho_history_[newest] * yy);
    }

    std::vector<double> r(dim);
    for (int j = 0; j < dim; ++j)
        r[j] = gamma * q[j];

    // Forward pass: from oldest to newest
    for (int k = 0; k < n_stored_; ++k) {
        int idx = (oldest_ + k) % memory_;
        double b = 0;
        for (int j = 0; j < dim; ++j)
            b += rho_history_[idx] * y_history_[idx][j] * r[j];
        for (int j = 0; j < dim; ++j)
            r[j] += s_history_[idx][j] * (alpha_buf[k] - b);
    }

    // direction = -H * g (descent direction)
    for (int j = 0; j < dim; ++j)
        direction_[j] = -r[j];

    // Save current state for next step's history update
    x_prev_ = x;
    g_prev_ = g;
    has_prev_ = true;

    // Apply step: x_new = x + alpha * direction
    double alpha = alpha0_;

    // Scale step by inverse max-force to keep steps reasonable
    double dir_max = 0;
    for (int j = 0; j < dim; ++j)
        dir_max = std::max(dir_max, std::abs(direction_[j]));
    if (dir_max > double(1e-20))
        alpha = std::min(alpha, double(0.1) / dir_max);  // cap step at 0.1 Angstrom

    // Remember the trial so the next step's search_direction can test Armijo.
    if (line_search_enabled_) {
        ls_alpha_ = alpha;
        ls_gd_ = 0;
        for (int j = 0; j < dim; ++j) ls_gd_ += g[j] * direction_[j];
        ls_origin_energy_ = f.energy;  // energy at x_prev_ (pre-step)
        ls_had_trial_ = true;
    }

    for (int i = 0; i < n; ++i) {
        f.atoms.x(i) += alpha * direction_[i*3+0];
        f.atoms.y(i) += alpha * direction_[i*3+1];
        f.atoms.z(i) += alpha * direction_[i*3+2];
    }

    ATV_DEBUG(
        "[lbfgs] applied step={} alpha={:.6e} dir_max={:.6e}",
        step_, alpha, dir_max);

    ++step_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
