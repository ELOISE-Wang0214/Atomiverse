/// @file src/cpu/minimize/anisotropic_cell_relax.cpp
#include "atomiverse/cpu/minimize/anisotropic_cell_relax.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {
namespace minimizer {

AnisotropicCellRelax::AnisotropicCellRelax(double P_target_ev_a3, double alpha)
    : P_target(P_target_ev_a3), alpha0(alpha) {
    if (!(alpha > 0.0) || !(max_scale > 0.0) || !(bulk_mod > 0.0))
        throw std::invalid_argument("AnisotropicCellRelax: invalid parameters");
}

void AnisotropicCellRelax::apply_scales_(Frame& f, double sx, double sy,
                                         double sz) const {
    if (!(sx > 0.0) || !(sy > 0.0) || !(sz > 0.0))
        throw std::runtime_error("AnisotropicCellRelax: invalid scale");
    f.simbox.scale_diagonal(sx, sy, sz);
    const int n = f.nlocal();
    auto x = f.span(kAtom_x);
    auto y = f.span(kAtom_y);
    auto z = f.span(kAtom_z);
    for (int i = 0; i < n; ++i) {
        x[i] *= sx;
        y[i] *= sy;
        z[i] *= sz;
    }
}

void AnisotropicCellRelax::build(const BuildContext&) {
    pressure_ = p_residual_ = 0;
    sx_ = sy_ = sz_ = 1;
    energy_prev_ = 0;
    step_ = 0;
    has_prev_ = false;
}

void AnisotropicCellRelax::launch(const LaunchContext& ctx) {
    auto& f = ctx.frame;
    const double V = f.simbox.volume();
    if (!(V > 0.0))
        throw std::runtime_error("AnisotropicCellRelax: non-positive volume");
    const auto& s = f.meta.get(kStress);

    if (ctx.at(ATV::minimizer::Phase::search_direction)) {
        // Per-axis pressure P_d = σ_dd / V (0 K).
        const double Px = s[0] / V;
        const double Py = s[1] / V;
        const double Pz = s[2] / V;
        pressure_ = (Px + Py + Pz) / 3.0;
        p_residual_ = std::max({std::fabs(Px - P_target), std::fabs(Py - P_target),
                                std::fabs(Pz - P_target)});
        auto one_s = [&](double P) {
            double ds = alpha0 * (P - P_target) / bulk_mod;
            ds = std::max(-max_scale, std::min(max_scale, ds));
            return 1.0 + ds;
        };
        sx_ = one_s(Px);
        sy_ = one_s(Py);
        sz_ = one_s(Pz);
        ++step_;
        return;
    }
    if (ctx.at(ATV::minimizer::Phase::line_search)) {
        if (step_ > 0)
            apply_scales_(f, sx_, sy_, sz_);
        return;
    }
    if (ctx.at(ATV::minimizer::Phase::convergence_check)) {
        const double Px = s[0] / V, Py = s[1] / V, Pz = s[2] / V;
        pressure_ = (Px + Py + Pz) / 3.0;
        p_residual_ = std::max({std::fabs(Px - P_target), std::fabs(Py - P_target),
                                std::fabs(Pz - P_target)});
        const double energy = f.meta.get(kEnergy);
        const double dE = has_prev_ ? std::fabs(energy - energy_prev_) : 1e30;
        if (ConvergenceMetrics* m = ctx.runtime.loop.metrics) {
            m->energy = energy;
            m->max_force = p_residual_;
            m->rms_force = p_residual_;
            m->max_displacement =
                std::max({std::fabs(sx_ - 1), std::fabs(sy_ - 1),
                          std::fabs(sz_ - 1)});
            m->rms_displacement = m->max_displacement;
            m->energy_change = dE;
            m->has_displacement = has_prev_;
        }
        energy_prev_ = energy;
        has_prev_ = true;
    }
}

} // namespace minimizer
} // namespace cpu
} // namespace ATV
