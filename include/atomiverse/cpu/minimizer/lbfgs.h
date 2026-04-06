/// @file cpu/minimizer/lbfgs.h
/// @brief minimizer::LBFGS — L-BFGS geometry optimizer (CPU).
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
///   Nocedal, Math. Comp. 35, 773 (1980) — L-BFGS algorithm
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
#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Per-step optimization record for diagnostics.
struct OptRecord {
    int step = 0;              ///< Optimization step number.
    Real energy = 0;           ///< Total energy [eV].
    Real max_force = 0;        ///< Maximum atomic force magnitude [eV/Angstrom].
    Real rms_force = 0;        ///< RMS of atomic force magnitudes [eV/Angstrom].
    Real max_displacement = 0; ///< Maximum displacement from previous step [Angstrom].
    Real rms_displacement = 0; ///< RMS displacement from previous step [Angstrom].
    Real energy_change = 0;    ///< Energy change from previous step [eV].
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
    /// @param max_ls    Maximum line search iterations (default 10).
    /// @param rms_ftol  RMS force tolerance [eV/Angstrom] (default 0.03).
    /// @param max_dtol  Max displacement tolerance [Angstrom] (default 0.002).
    /// @param rms_dtol  RMS displacement tolerance [Angstrom] (default 0.001).
    /// @param etol      Energy change tolerance [eV] (default 1e-5).
    LBFGSStep(int memory = 5, Real alpha0 = Real(0.01),
              Real ftol = Real(0.05), int max_ls = 10,
              Real rms_ftol = Real(0.03), Real max_dtol = Real(0.002),
              Real rms_dtol = Real(0.001), Real etol = Real(1e-5));

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
    Real max_force() const { return max_force_; }

    /// @brief Per-step optimization trace.
    const std::vector<OptRecord>& trace() const { return trace_; }

private:
    int memory_;
    Real alpha0_;
    Real ftol_;
    int max_ls_;
    Real rms_ftol_;
    Real max_dtol_;
    Real rms_dtol_;
    Real etol_;

    // L-BFGS history (circular buffer)
    std::vector<std::vector<Real>> s_history_;  // position differences
    std::vector<std::vector<Real>> y_history_;  // gradient differences
    std::vector<Real> rho_history_;             // 1 / (y^T s)
    int n_stored_ = 0;
    int oldest_ = 0;

    // Previous step state
    std::vector<Real> x_prev_;
    std::vector<Real> g_prev_;
    std::vector<Real> direction_;  // search direction
    bool has_prev_ = false;
    bool converged_ = false;
    Real max_force_ = 0;
    Real energy_prev_ = 0;
    int step_ = 0;

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
