/// @file cpu/minimizer/lbfgs.cpp
/// @brief L-BFGS step implementation: two-loop recursion + multi-criteria convergence.

#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/core/units.h"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace ATV {
namespace cpu {
namespace minimizer {

LBFGSStep::LBFGSStep(int memory, double alpha0, double ftol, int max_ls,
                      double rms_ftol, double max_dtol, double rms_dtol, double etol)
    : memory_(memory), alpha0_(alpha0), ftol_(ftol), max_ls_(max_ls),
      rms_ftol_(rms_ftol), max_dtol_(max_dtol), rms_dtol_(rms_dtol), etol_(etol) {}

void LBFGSStep::build(const cpu::BuildContext& ctx) {
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
    trace_.clear();
}

void LBFGSStep::launch(const cpu::LaunchContext& ctx) {
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

        // Compute displacement metrics (relative to previous step)
        double max_disp = 0;
        double sum_d2 = 0;
        if (has_prev_) {
            for (int i = 0; i < n; ++i) {
                double dx = f.atoms.x(i) - x_prev_[i*3+0] / 1.0; // x_prev_ stored as flat array... wait
                // Actually x_prev_ stores positions from the search_direction phase.
                // But positions may have changed. Need to store the pre-step positions.
                // For now, compute from direction_ * alpha (the step taken).
            }
        }

        // Use displacement from direction_ (the actual step taken)
        if (step_ > 0) {
            for (int i = 0; i < n; ++i) {
                for (int d = 0; d < 3; ++d) {
                    double disp = f.atoms.x(i) - x_prev_[i*3+0]; // compare current to pre-step
                    if (d == 1) disp = f.atoms.y(i) - x_prev_[i*3+1];
                    if (d == 2) disp = f.atoms.z(i) - x_prev_[i*3+2];
                    sum_d2 += disp * disp;
                    double abs_disp = std::abs(disp);
                    if (abs_disp > max_disp) max_disp = abs_disp;
                }
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

    for (int i = 0; i < n; ++i) {
        f.atoms.x(i) += alpha * direction_[i*3+0];
        f.atoms.y(i) += alpha * direction_[i*3+1];
        f.atoms.z(i) += alpha * direction_[i*3+2];
    }

    ++step_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
