/// @file cpu/minimizer/lbfgs.h
/// @brief minimizer::LBFGS -- L-BFGS geometry optimizer (CPU).
///
/// Two-loop recursion L-BFGS with optional line search.
/// Compose with SCFEvaluator at force phase for ab initio geometry optimization.
///
/// Usage:
/// @code
///   cpu::minimizer::LBFGS opt;
///   opt.add(cpu::minimizer::LBFGSStep{memory, alpha, ftol, max_ls});
///   opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
///   opt.on(ATV::minimizer::Phase::end_of_step).add(Thermo{1});
///   opt.build();
///   opt.run(frame, 100);
/// @endcode
///
/// References:
///   Nocedal, Math. Comp. 35, 773 (1980) -- L-BFGS algorithm
///   Nocedal & Wright, Numerical Optimization, Ch. 7
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/minimizer/init.h"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Per-step optimization record for diagnostics.
struct OptRecord {
    int step = 0;              ///< Optimization step number.
    double energy = 0;           ///< Total energy [eV].
    double max_force = 0;        ///< Maximum atomic force magnitude [eV/Angstrom].
    double rms_force = 0;        ///< RMS of atomic force magnitudes [eV/Angstrom].
    double max_displacement = 0; ///< Maximum displacement from previous step [Angstrom].
    double rms_displacement = 0; ///< RMS displacement from previous step [Angstrom].
    double energy_change = 0;    ///< Energy change from previous step [eV].
};

/// @brief L-BFGS step component.
///
/// Implements the two-loop recursion to compute search direction, then
/// performs a backtracking line search with Armijo condition.
/// Binds to search_direction + convergence_check phases.
///
/// Convergence requires ALL five criteria to be satisfied:
///   - max force < ftol
///   - RMS force < rms_ftol
///   - max displacement < max_dtol
///   - RMS displacement < rms_dtol
///   - |energy change| < etol
///
/// At step 0: steepest descent (no history). From step 1: L-BFGS direction.
class LBFGSStep : public cpu::ComponentBase {
public:
    /// @brief Construct L-BFGS step.
    /// @param memory    Number of correction pairs to store (default 5).
    /// @param alpha0    Initial step size [Angstrom^2/eV] (default 0.01).
    /// @param ftol      Max force convergence tolerance [eV/Angstrom] (default 0.05).
    /// @param max_ls    Maximum backtracking attempts per step (default 10).
    /// @param rms_ftol  RMS force tolerance [eV/Angstrom] (default 0.03).
    /// @param max_dtol  Max displacement tolerance [Angstrom] (default 0.002).
    /// @param rms_dtol  RMS displacement tolerance [Angstrom] (default 0.001).
    /// @param etol      Energy change tolerance [eV] (default 1e-5).
    LBFGSStep(int memory = 5, double alpha0 = 0.01,
              double ftol = 0.05, int max_ls = 10,
              double rms_ftol = 0.03, double max_dtol = 0.002,
              double rms_dtol = 0.001, double etol = 1e-5);

    /// @brief Enable the backtracking line search (Armijo condition).
    ///
    /// When enabled, each driver step that fails `E_{k+1} - E_k <= c1 * alpha_k
    /// * g_k.d_k` is *rolled back*: the positions return to `x_k`, alpha is
    /// scaled by `contract`, and the trial is re-taken from `x_k` along the
    /// same direction. History (s, y pairs) is only updated on accepted
    /// trials. `min_alpha` caps the shrink to avoid infinite loops.
    ///
    /// Must be called before `opt.build()` for the setting to take effect.
    /// @param c1         Armijo constant in (0, 1) (default 1e-4).
    /// @param contract   Backtracking contraction factor (default 0.5).
    /// @param min_alpha  Minimum alpha before accepting the trial even on
    ///                   Armijo failure (default 1e-6).
    /// @return Reference to self for chaining.
    LBFGSStep& line_search(double c1 = 1e-4,
                           double contract = 0.5,
                           double min_alpha = 1e-6);

    std::string name() const override { return "lbfgs_step"; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {
            {static_cast<int>(ATV::minimizer::Phase::search_direction), 0},
            {static_cast<int>(ATV::minimizer::Phase::convergence_check), 0}
        };
    }

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Check if optimization has converged.
    bool converged() const { return converged_; }

    /// @brief Current max force magnitude.
    double max_force() const { return max_force_; }

    /// @brief Per-step optimization trace.
    const std::vector<OptRecord>& trace() const { return trace_; }

private:
    int memory_;
    double alpha0_;
    double ftol_;
    int max_ls_;
    double rms_ftol_;
    double max_dtol_;
    double rms_dtol_;
    double etol_;

    // L-BFGS history (circular buffer)
    std::vector<std::vector<double>> s_history_;  // position differences
    std::vector<std::vector<double>> y_history_;  // gradient differences
    std::vector<double> rho_history_;             // 1 / (y^T s)
    int n_stored_ = 0;
    int oldest_ = 0;

    // Previous step state
    std::vector<double> x_prev_;
    std::vector<double> g_prev_;
    std::vector<double> direction_;  // search direction
    bool has_prev_ = false;
    bool converged_ = false;
    double max_force_ = 0;
    double energy_prev_ = 0;
    int step_ = 0;

    // Line search (optional; off by default to preserve legacy behaviour).
    bool line_search_enabled_ = false;
    double ls_c1_ = 1e-4;        // Armijo constant.
    double ls_contract_ = 0.5;   // Backtracking factor.
    double ls_min_alpha_ = 1e-6; // Floor on alpha.
    double ls_alpha_ = 0;        // Alpha chosen for the most recent trial.
    double ls_gd_ = 0;           // g_k · d_k at the trial's origin point.
    double ls_origin_energy_ = 0;// E at the trial's origin (x_prev_).
    bool   ls_had_trial_ = false;

    // Optimization trace
    std::vector<OptRecord> trace_;
};

/// @brief L-BFGS minimizer driver.
///
/// Auto-adds ZeroForces and ZeroEnergy at force phase.
class LBFGS : public Driver<ATV::minimizer::Phase, cpu::device_tag> {
public:
    LBFGS() {
        this->add(ZeroForces{});
        this->add(ZeroEnergy{});
    }
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
