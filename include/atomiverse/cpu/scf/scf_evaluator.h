/// @file cpu/scf/scf_evaluator.h
/// @brief SCF — CPU force-phase component for ab initio MD.
///
/// Single component for all mean-field methods (RHF, UHF, RKS, UKS).
/// Binds to md::Phase::force. Each launch(): rebuild integrals, run SCF,
/// compute gradient, write forces to frame.
///
/// For AIMD, compose with md::MD driver + VelocityVerlet:
/// @code
///   md::MD sim;
///   sim.add(md::VelocityVerlet{dt});
///   sim.on(md::Phase::force).add(scf::SCF{config, charge});
///   sim.build();
///   sim.run(frame, n_steps);
/// @endcode
#pragma once

#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "atomiverse/cpu/component/hf.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"
#include "atomiverse/cpu/scf/scf_types.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Shared view into the latest SCF state.
///
/// Bundles the three shared-pointer handles published by `SCF`.
/// Observer components capture one of these during construction and then
/// read from it each step; the pointers stay valid across the driver's
/// value-copy of the evaluator.
struct SCFResultHandle {
    std::shared_ptr<const SCFResult> result; ///< Latest SCF result.
    std::shared_ptr<const IntegralPack>
        integrals;                     ///< Latest integrals (shells + atoms).
    std::shared_ptr<const bool> ready; ///< True once a launch has completed.
    std::shared_ptr<const SCFConfig>
        config; ///< SCF configuration (spin, method, grid).
};

/// @brief Ab initio force evaluator component for MD drivers.
///
/// Rebuilds integrals + runs SCF + computes analytic gradient each step.
/// Writes forces to `kAtom_fx/fy/fz` and energy to `kEnergy`.
/// Optionally uses ASPC predictor-corrector for density matrix warm-start.
class SCF : public cpu::HF {
  public:
    /// @brief Construct from SCF configuration.
    /// @param config  SCF parameters (method, spin, XC, convergence, grid).
    /// @param basis   Basis set to resolve at each step from frame positions.
    /// @param charge  Molecular charge (0 = neutral).
    SCF(SCFConfig config, BasisSet basis = BasisSet::STO3G, int charge = 0);

    /// @brief Enable ASPC predictor-corrector for density warm-start.
    /// @param order           Extrapolation order K (default 3).
    /// @param corrector_steps Max SCF corrector iterations (default 2).
    /// @param corrector_tol   Corrector convergence tolerance (default 1e-4).
    /// @param fallback        Fall back to full SCF if corrector fails.
    /// @return Reference to this for chaining (preserves value category:
    ///         lvalue chain → `SCF&`, rvalue chain →
    ///         `SCF&&` so `driver.add(...)` deduces an rvalue
    ///         and moves the temporary into ownership).
    SCF& predictor(int order = 3, int corrector_steps = 2,
                   double corrector_tol = 1e-4, bool fallback = true) &;
    SCF&& predictor(int order = 3, int corrector_steps = 2,
                    double corrector_tol = 1e-4, bool fallback = true) &&;

    /// @brief Set a per-step callback for diagnostics.
    /// @param cb Callback invoked after each force evaluation.
    /// @return Reference to this for chaining (value-category preserving).
    SCF& on_evaluate(std::function<void(const SCFResult&)> cb) &;
    SCF&& on_evaluate(std::function<void(const SCFResult&)> cb) &&;

    /// @brief Toggle gradient computation.
    ///
    /// `with_gradient(true)`  (default) — `launch()` runs SCF + analytic
    ///                                    nuclear gradient and writes
    ///                                    forces into the force columns.
    ///                                    `capabilities()` includes
    ///                                    `Capability::Force`.
    /// `with_gradient(false)` — `launch()` runs SCF only; forces are not
    ///                          touched. `capabilities()` drops
    ///                          `Capability::Force`. Use for energy-only
    ///                          single-point work in `cpu::Evaluate`.
    ///
    /// Must be called before `driver.build()` (capabilities are inspected
    /// during validation).
    /// @return Reference to self for chaining (value-category preserving).
    SCF& with_gradient(bool enabled) & noexcept {
        compute_force_ = enabled;
        return *this;
    }
    SCF&& with_gradient(bool enabled) && noexcept {
        compute_force_ = enabled;
        return std::move(*this);
    }

    std::string name() const override { return "scf_evaluator"; }
    /* role via Potential/kind base */
    /// @brief Capabilities reflect what `launch()` writes to the frame.
    ///        Always: Energy, DensityMatrix, Orbitals.
    ///        Force is added when gradient computation is enabled (the
    ///        default; toggle via `with_gradient(false)`).

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

    /// @brief Whether a full SCF result is available from a completed step.
    bool has_last_result() const { return *has_last_result_; }

    /// @brief Most recent SCFResult produced by launch().
    /// @pre has_last_result() == true.
    const SCFResult& last_result() const { return *last_result_; }

    /// @brief Most recent IntegralPack rebuilt for the current geometry.
    /// @pre has_last_result() == true.
    const IntegralPack& last_integrals() const { return *last_integrals_; }

    /// @brief SCF configuration (method, spin, grid, etc.).
    const SCFConfig& config() const { return config_; }

    /// @brief Shared handle to the latest SCFResult.
    ///
    /// Returns a shared_ptr that stays valid across driver ownership moves
    /// (the driver copies the SCF by value when adding components).
    /// Observers such as `OrbitalInfoObserver` use this to read state written
    /// by the driver-owned evaluator.
    ///
    /// @note `*has_last_result_handle()` is initially false; check it before
    ///       dereferencing `last_result_handle()`.
    std::shared_ptr<const SCFResult> last_result_handle() const {
        return last_result_;
    }

    /// @brief Shared handle to the latest IntegralPack. @see last_result_handle.
    std::shared_ptr<const IntegralPack> last_integrals_handle() const {
        return last_integrals_;
    }

    /// @brief Shared readiness flag (true after the first converged launch()).
    std::shared_ptr<const bool> has_last_result_handle() const {
        return has_last_result_;
    }

    /// @brief Bundled handle for observer components.
    SCFResultHandle result_handle() const {
        return {last_result_, last_integrals_, has_last_result_,
                config_shared_};
    }

  private:
    /// @brief Install the ASPC predictor — shared body for the two
    ///        ref-qualified `predictor()` overloads.
    void set_predictor_(int order, int corrector_steps, double corrector_tol,
                        bool fallback);

    /// @brief Install the per-step callback — shared body for the two
    ///        ref-qualified `on_evaluate()` overloads.
    void set_on_evaluate_(std::function<void(const SCFResult&)> cb);

    SCFConfig config_;
    BasisSet basis_;
    int charge_;

    bool compute_force_ = true;

    bool use_predictor_ = false;
    int corrector_steps_ = 2;
    double corrector_tol_ = 1e-4;
    bool fallback_full_scf_ = true;
    std::unique_ptr<ASPCPredictor> aspc_alpha_;
    std::unique_ptr<ASPCPredictor> aspc_beta_;
    bool has_last_density_ = false;
    xt::xtensor<double, 2> last_P_alpha_;
    xt::xtensor<double, 2> last_P_beta_;

    /// Pre-built inner SCF solver. Constructed in build() with a
    /// representative IntegralPack: the n_basis-sized eigen/density
    /// workspaces are allocated once; geometry-dependent structures
    /// (X, shell pairs, KS grid) are rebuilt inside each run() because
    /// atoms move between launches.
    std::unique_ptr<SCFSolver> solver_;

    /// Pre-built corrector solver for the ASPC path (tighter tol, capped
    /// iterations — a different SCFConfig from solver_). Constructed in
    /// build() when the predictor is enabled so launch() never allocates.
    std::unique_ptr<SCFSolver> corrector_solver_;

    // Shared-pointer storage so observers can capture handles that stay valid
    // after the driver copies this evaluator into its owned-component list.
    std::shared_ptr<bool> has_last_result_{std::make_shared<bool>(false)};
    std::shared_ptr<SCFResult> last_result_{std::make_shared<SCFResult>()};
    std::shared_ptr<IntegralPack> last_integrals_{
        std::make_shared<IntegralPack>()};
    std::shared_ptr<SCFConfig> config_shared_; // initialized in ctor

    std::function<void(const SCFResult&)> on_evaluate_;
};

/* SCF */;
} // namespace scf
} // namespace cpu
} // namespace ATV
