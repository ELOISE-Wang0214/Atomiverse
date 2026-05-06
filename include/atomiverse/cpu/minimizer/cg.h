/// @file cpu/minimizer/cg.h
/// @brief CG — Polak-Ribière conjugate-gradient Optimizer component.
///
/// Pluggable alternative to `cpu::minimizer::LBFGS` inside
/// `cpu::minimizer::Minimize`. Same outer contract: binds
/// `search_direction` (computes the step) and `convergence_check`
/// (reports `outer_status.converged` to the driver).
///
/// Algorithm (Polak-Ribière+):
///   d_0 = -g_0
///   For k ≥ 1:
///     β = max(0, (g_k · (g_k − g_{k−1})) / (g_{k−1} · g_{k−1}))
///     d_k = -g_k + β · d_{k−1}
///   x_{k+1} = x_k + α · d_k
///
/// Restart: when β = 0 (numerator non-positive), the direction collapses
/// back to steepest descent. No explicit periodic restart in this
/// implementation — the Polak-Ribière+ projection onto non-negative β
/// already handles ill-conditioned regions.
///
/// Convergence rule mirrors `LBFGS` exactly (max/RMS force tolerances,
/// max/RMS displacement tolerances, energy change). This makes the
/// optimizer hot-swap test (`tests/cpu/test_minimize_driver_swap_optimizer.cpp`)
/// meaningful: both classes report convergence under the same criteria.
///
/// References:
///   Polak & Ribière, ESAIM Math. Modelling Numer. Anal. 3, 35 (1969)
///   Nocedal & Wright, Numerical Optimization, Ch. 5
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"  // for OptRecord
#include "atomiverse/facade/component_role.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Polak-Ribière conjugate-gradient Optimizer.
class CG : public cpu::ComponentBase {
public:
    /// @brief Construct CG optimizer.
    /// @param alpha0    Step size [Å²/eV] (default 0.01).
    /// @param ftol      Max force convergence tolerance [eV/Å] (default 0.05).
    /// @param rms_ftol  RMS force tolerance [eV/Å] (default 0.03).
    /// @param max_dtol  Max displacement tolerance [Å] (default 0.002).
    /// @param rms_dtol  RMS displacement tolerance [Å] (default 0.001).
    /// @param etol      Energy change tolerance [eV] (default 1e-5).
    explicit CG(double alpha0 = 0.01,
                double ftol = 0.05,
                double rms_ftol = 0.03,
                double max_dtol = 0.002,
                double rms_dtol = 0.001,
                double etol = 1e-5);

    std::string name() const override { return "cg"; }
    ComponentRole role() const override { return ComponentRole::Optimizer; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {
            {static_cast<int>(ATV::minimizer::Phase::search_direction), 0},
            {static_cast<int>(ATV::minimizer::Phase::convergence_check), 0}
        };
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    bool converged() const { return converged_; }
    double max_force() const { return max_force_; }
    const std::vector<OptRecord>& trace() const { return trace_; }

private:
    double alpha0_;
    double ftol_, rms_ftol_, max_dtol_, rms_dtol_, etol_;

    // Previous-step state.
    std::vector<double> x_prev_;
    std::vector<double> g_prev_;
    std::vector<double> d_prev_;     ///< Previous CG search direction.
    bool   has_prev_ = false;
    bool   converged_ = false;
    double max_force_ = 0;
    double energy_prev_ = 0;
    int    step_ = 0;

    std::vector<OptRecord> trace_;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
