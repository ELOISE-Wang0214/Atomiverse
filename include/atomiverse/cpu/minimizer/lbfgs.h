/// @file cpu/minimizer/lbfgs.h
/// @brief minimizer::LBFGS — L-BFGS Optimizer component for the
///        `cpu::minimizer::Minimize` driver (declared in `minimize.h`).
///
/// Two-loop L-BFGS recursion with optional Armijo line search.
/// `Minimize` is the OptimizeDriver subclass; `LBFGS` is the Optimizer
/// component plugged into it.
///
/// **MPI / domain decomposition.** Under multi-rank the component
/// operates on the owned rows `[0, dd->n_local)` only; every inner
/// product of the two-loop recursion, the step-size cap, the Armijo
/// directional derivative, and the convergence metrics are reduced
/// across ranks (`minimizer/global_reduce.h`), so all ranks take the
/// identical step and reach the identical convergence verdict. Atom
/// migration reshuffles local row indexing, so on every
/// `dd->neighbor_generation` bump the curvature history and the
/// displacement reference are discarded (the next step restarts as
/// steepest descent) — convergence is unaffected, only its rate.
///
/// Usage:
/// @code
///   cpu::minimizer::Minimize opt;
///   opt.add(cpu::minimizer::LBFGS{memory, alpha, ftol, max_ls});
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

#include <cstdint>
#include <string>
#include <vector>

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/driver/backend_traits.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/minimizer/line_search.h"
#include "atomiverse/cpu/minimizer/minimize.h"
#include "atomiverse/cpu/driver/optimize_driver.h"
#include "atomiverse/cpu/async/snapshot_publisher.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace minimizer {

/// @brief Per-step optimization record for diagnostics.
struct OptRecord {
    int step = 0;         ///< Optimization step number.
    double energy = 0;    ///< Total energy [eV].
    double max_force = 0; ///< Maximum atomic force magnitude [eV/Angstrom].
    double rms_force = 0; ///< RMS of atomic force magnitudes [eV/Angstrom].
    double max_displacement =
        0; ///< Maximum displacement from previous step [Angstrom].
    double rms_displacement =
        0; ///< RMS displacement from previous step [Angstrom].
    double energy_change = 0; ///< Energy change from previous step [eV].
};

/// @brief L-BFGS Optimizer component.
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
class LBFGS : public cpu::ComponentBase {
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
    LBFGS(int memory = 5, double alpha0 = 0.01, double ftol = 0.05,
          int max_ls = 10, double rms_ftol = 0.03, double max_dtol = 0.002,
          double rms_dtol = 0.001, double etol = 1e-5);

    /// @brief Select the line-search policy used at the `line_search` phase.
    ///
    /// In `StrongWolfe` (default) and `Armijo` modes, each optimization step
    /// performs an in-step line search along the quasi-Newton direction
    /// `p = -H g`: trial geometries `x_k + alpha p` are evaluated by replaying
    /// the force phase through `Context::reeval_forces`, and `alpha` is the
    /// dimensionless multiplier on `p` (initial trial 1.0). `StrongWolfe`
    /// enforces N&W conditions 3.7a+3.7b (so `yᵀs > 0` is automatic and every
    /// stored curvature pair is valid); `Armijo` enforces only 3.6a (energy
    /// per trial, no gradient). `Fixed` reproduces the legacy step
    /// `x -= min(alpha0, cap)·(-p)` with no trial evaluation.
    ///
    /// Must be called before `opt.build()`.
    /// @param mode      Fixed / Armijo / StrongWolfe (default StrongWolfe).
    /// @param c1        Sufficient-decrease constant in (0,1) (default 1e-4).
    /// @param c2        Wolfe curvature constant, c1<c2<1 (default 0.9).
    /// @param max_step  Max trial displacement [Å]; caps alpha (default 0.2).
    /// @param min_alpha Alpha floor before giving up (default 1e-8).
    /// @return Reference to self for chaining.
    LBFGS& line_search(LineSearchMode mode = LineSearchMode::StrongWolfe,
                       double c1 = 1e-4, double c2 = 0.9, double max_step = 0.2,
                       double min_alpha = 1e-8);

    std::string name() const override { return "lbfgs"; }
    ComponentRole role() const override { return ComponentRole::Optimizer; }

    std::vector<DefaultBinding> default_bindings() const override {
        return {
            {static_cast<int>(ATV::minimizer::Phase::search_direction), 0},
            {static_cast<int>(ATV::minimizer::Phase::line_search), 0},
            {static_cast<int>(ATV::minimizer::Phase::convergence_check), 0}};
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
    std::vector<std::vector<double>> s_history_; // position differences
    std::vector<std::vector<double>> y_history_; // gradient differences
    std::vector<double> rho_history_;            // 1 / (y^T s)
    int n_stored_ = 0;
    int oldest_ = 0;

    // Previous step state
    std::vector<double> x_prev_;
    std::vector<double> g_prev_;
    std::vector<double> direction_; // search direction
    bool has_prev_ = false;
    bool x_prev_valid_ = false; ///< False after migration reshuffles rows.

    // Per-step workspaces, sized once in build() — no launch() allocation.
    std::vector<double> x_ws_;
    std::vector<double> g_ws_;
    std::vector<double> q_ws_;
    std::vector<double> r_ws_;
    std::vector<double> alpha_ws_;

    /// @brief Last observed `dd->neighbor_generation`; a bump means atom
    /// migration reshuffled local rows → history must be discarded.
    std::uint64_t last_generation_ = ~std::uint64_t{0};
    bool converged_ = false;
    double max_force_ = 0;
    double energy_prev_ = 0;
    int step_ = 0;

    // Line search. Default StrongWolfe; `line_search()` overrides.
    LineSearchMode ls_mode_ = LineSearchMode::StrongWolfe;
    double ls_c1_ = 1e-4;        // Sufficient-decrease constant (3.6a/3.7a).
    double ls_c2_ = 0.9;         // Wolfe curvature constant (3.7b); 0.9 for QN.
    double ls_contract_ = 0.5;   // Armijo backtracking factor.
    double ls_max_step_ = 0.2;   // Max trial displacement [Å] (caps alpha).
    double ls_min_alpha_ = 1e-8; // Alpha floor.
    // Carry-over from the search_direction phase to the line_search phase
    // within a single optimization step (all globally reduced under MPI):
    double ls_phi0_ = 0;   // phi(0) = energy at x_prev_ [eV].
    double ls_gd0_ = 0;    // phi'(0) = g . direction (< 0 on descent).
    double ls_dir_max_ = 0; // global max|direction| component, for the cap.

    // Optimization trace
    std::vector<OptRecord> trace_;
};

} // namespace minimizer
} // namespace cpu
} // namespace ATV
