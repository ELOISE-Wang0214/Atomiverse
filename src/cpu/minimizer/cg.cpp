/// @file cpu/minimizer/cg.cpp
/// @brief Polak-Ribière+ CG step + multi-criteria convergence (mirrors LBFGS).

#include "atomiverse/cpu/minimizer/cg.h"
#include "atomiverse/core/log.h"

#include <algorithm>
#include <cmath>

namespace ATV {
namespace cpu {
namespace minimizer {

CG::CG(double alpha0, double ftol, double rms_ftol,
       double max_dtol, double rms_dtol, double etol)
    : alpha0_(alpha0), ftol_(ftol), rms_ftol_(rms_ftol),
      max_dtol_(max_dtol), rms_dtol_(rms_dtol), etol_(etol) {}

void CG::build(const cpu::BuildContext& ctx) {
    int dim = ctx.frame.atoms.size() * 3;
    x_prev_.assign(dim, 0.0);
    g_prev_.assign(dim, 0.0);
    d_prev_.assign(dim, 0.0);
    has_prev_ = false;
    converged_ = false;
    step_ = 0;
    trace_.clear();
}

void CG::launch(const cpu::LaunchContext& ctx) {
    auto& f = ctx.frame;
    const int n = static_cast<int>(f.atoms.size());
    const int dim = n * 3;
    const int phase = ctx.stage;

    // ── Convergence check ──────────────────────────────────────────────
    if (phase == static_cast<int>(ATV::minimizer::Phase::convergence_check)) {
        max_force_ = 0;
        double sum_f2 = 0;
        for (int i = 0; i < n; ++i) {
            const double fx = f.atoms.fx(i);
            const double fy = f.atoms.fy(i);
            const double fz = f.atoms.fz(i);
            const double mag = std::sqrt(fx*fx + fy*fy + fz*fz);
            if (mag > max_force_) max_force_ = mag;
            sum_f2 += fx*fx + fy*fy + fz*fz;
        }
        const double rms_force = std::sqrt(sum_f2 / (3 * n));

        double max_disp = 0, sum_d2 = 0;
        if (step_ > 0) {
            for (int i = 0; i < n; ++i) {
                const double dx = f.atoms.x(i) - x_prev_[i*3+0];
                const double dy = f.atoms.y(i) - x_prev_[i*3+1];
                const double dz = f.atoms.z(i) - x_prev_[i*3+2];
                sum_d2 += dx*dx + dy*dy + dz*dz;
                max_disp = std::max(max_disp, std::max({std::abs(dx),
                                                         std::abs(dy),
                                                         std::abs(dz)}));
            }
        }
        const double rms_disp = (step_ > 0) ? std::sqrt(sum_d2 / (3 * n)) : 0.0;
        const double energy_change = (step_ > 0)
                                   ? std::abs(f.energy - energy_prev_)
                                   : 1e30;

        OptRecord rec{};
        rec.step = step_;
        rec.energy = f.energy;
        rec.max_force = max_force_;
        rec.rms_force = rms_force;
        rec.max_displacement = max_disp;
        rec.rms_displacement = rms_disp;
        rec.energy_change = energy_change;
        trace_.push_back(rec);

        converged_ = (max_force_ < ftol_) &&
                     (rms_force < rms_ftol_) &&
                     (step_ == 0 || max_disp     < max_dtol_) &&
                     (step_ == 0 || rms_disp     < rms_dtol_) &&
                     (step_ == 0 || energy_change < etol_);

        ATV_INFO(
            "[cg] step={} energy_eV={:.6f} max_force_eVA={:.6f} "
            "rms_force_eVA={:.6f} max_disp_A={:.6e} rms_disp_A={:.6e} "
            "dE_eV={:.6e} converged={}",
            step_, f.energy, max_force_, rms_force, max_disp, rms_disp,
            energy_change, converged_);

        if (converged_) {
            ctx.runtime.outer_status.converged = true;
            ctx.runtime.outer_status.reason = "cg multi-criteria convergence";
        }

        energy_prev_ = f.energy;
        return;
    }

    // ── Search direction (Polak-Ribière+) ──────────────────────────────
    std::vector<double> g(dim), x(dim);
    for (int i = 0; i < n; ++i) {
        x[i*3+0] = f.atoms.x(i);
        x[i*3+1] = f.atoms.y(i);
        x[i*3+2] = f.atoms.z(i);
        g[i*3+0] = -f.atoms.fx(i);
        g[i*3+1] = -f.atoms.fy(i);
        g[i*3+2] = -f.atoms.fz(i);
    }

    std::vector<double> d(dim);
    if (!has_prev_) {
        // Step 0: steepest descent.
        for (int k = 0; k < dim; ++k) d[k] = -g[k];
    } else {
        // β = max(0, (g · (g - g_prev)) / (g_prev · g_prev))
        double num = 0, den = 0;
        for (int k = 0; k < dim; ++k) {
            num += g[k] * (g[k] - g_prev_[k]);
            den += g_prev_[k] * g_prev_[k];
        }
        double beta = (den > 0.0) ? std::max(0.0, num / den) : 0.0;
        for (int k = 0; k < dim; ++k) d[k] = -g[k] + beta * d_prev_[k];
    }

    // Take the step: x_{k+1} = x + α · d
    for (int i = 0; i < n; ++i) {
        f.atoms.x(i) = x[i*3+0] + alpha0_ * d[i*3+0];
        f.atoms.y(i) = x[i*3+1] + alpha0_ * d[i*3+1];
        f.atoms.z(i) = x[i*3+2] + alpha0_ * d[i*3+2];
    }

    // Save state for next iteration.
    x_prev_ = std::move(x);
    g_prev_ = g;
    d_prev_ = std::move(d);
    has_prev_ = true;
    ++step_;
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
