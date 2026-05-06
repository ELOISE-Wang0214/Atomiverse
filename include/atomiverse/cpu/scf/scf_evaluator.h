/// @file cpu/scf/scf_evaluator.h
/// @brief SCFEvaluator — CPU force-phase component for ab initio MD.
///
/// Single component for all mean-field methods (RHF, UHF, RKS, UKS).
/// Binds to md::Phase::force. Each launch(): rebuild integrals, run SCF,
/// compute gradient, write forces to frame.
///
/// For AIMD, compose with md::MD driver + VelocityVerlet:
/// @code
///   md::MD sim;
///   sim.add(md::VelocityVerlet{dt});
///   sim.on(md::Phase::force).add(scf::SCFEvaluator{config, charge});
///   sim.build();
///   sim.run(frame, n_steps);
/// @endcode
#pragma once

#include <memory>
#include <vector>
#include <functional>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Shared view into the latest SCFEvaluator state.
///
/// Bundles the three shared-pointer handles published by `SCFEvaluator`.
/// Observer components capture one of these during construction and then
/// read from it each step; the pointers stay valid across the driver's
/// value-copy of the evaluator.
struct SCFResultHandle {
    std::shared_ptr<const SCFResult> result;       ///< Latest SCF result.
    std::shared_ptr<const IntegralPack> integrals; ///< Latest integrals (shells + atoms).
    std::shared_ptr<const bool> ready;             ///< True once a launch has completed.
    std::shared_ptr<const SCFConfig> config;       ///< SCF configuration (spin, method, grid).
};

/// @brief Ab initio force evaluator component for MD drivers.
///
/// Rebuilds integrals + runs SCF + computes analytic gradient each step.
/// Writes forces to frame.atoms.fx/fy/fz and energy to frame.energy.
/// Optionally uses ASPC predictor-corrector for density matrix warm-start.
class SCFEvaluator : public cpu::ComponentBase {
public:
    /// @brief Construct from SCF configuration.
    /// @param config  SCF parameters (method, spin, XC, convergence, grid).
    /// @param basis   Basis set to resolve at each step from frame positions.
    /// @param charge  Molecular charge (0 = neutral).
    SCFEvaluator(SCFConfig config,
                 BasisSet basis = BasisSet::STO3G,
                 int charge = 0);

    /// @brief Enable ASPC predictor-corrector for density warm-start.
    /// @param order           Extrapolation order K (default 3).
    /// @param corrector_steps Max SCF corrector iterations (default 2).
    /// @param corrector_tol   Corrector convergence tolerance (default 1e-4).
    /// @param fallback        Fall back to full SCF if corrector fails.
    /// @return Reference to this for chaining.
    SCFEvaluator& predictor(int order = 3, int corrector_steps = 2,
                            double corrector_tol = 1e-4, bool fallback = true);

    /// @brief Set a per-step callback for diagnostics.
    /// @param cb Callback invoked after each force evaluation.
    SCFEvaluator& on_evaluate(std::function<void(const SCFResult&)> cb);

    std::string name() const override { return "scf_evaluator"; }
    ComponentRole role() const override { return ComponentRole::Evaluator; }
    /// @brief Capabilities reflect what SCFEvaluator writes each launch.
    ///        Always Energy + DensityMatrix + Orbitals; Force is included
    ///        because gradient computation is currently unconditional.
    ///        Phase 4 introduces a `with_gradient(false)` mode that drops
    ///        the Force bit for energy-only single points.
    Capability capabilities() const override {
        return Capability::Energy | Capability::Force |
               Capability::DensityMatrix | Capability::Orbitals;
    }

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
    /// (the driver copies the SCFEvaluator by value when adding components).
    /// Observers such as `OrbitalInfoObserver` use this to read state written
    /// by the driver-owned evaluator.
    ///
    /// @note `*has_last_result_handle()` is initially false; check it before
    ///       dereferencing `last_result_handle()`.
    std::shared_ptr<const SCFResult> last_result_handle() const { return last_result_; }

    /// @brief Shared handle to the latest IntegralPack. @see last_result_handle.
    std::shared_ptr<const IntegralPack> last_integrals_handle() const { return last_integrals_; }

    /// @brief Shared readiness flag (true after the first converged launch()).
    std::shared_ptr<const bool> has_last_result_handle() const { return has_last_result_; }

    /// @brief Bundled handle for observer components.
    SCFResultHandle result_handle() const {
        return {last_result_, last_integrals_, has_last_result_, config_shared_};
    }

private:
    SCFConfig config_;
    BasisSet basis_;
    int charge_;

    bool use_predictor_ = false;
    int corrector_steps_ = 2;
    double corrector_tol_ = 1e-4;
    bool fallback_full_scf_ = true;
    std::unique_ptr<ASPCPredictor> aspc_alpha_;
    std::unique_ptr<ASPCPredictor> aspc_beta_;
    bool has_last_density_ = false;
    xt::xtensor<double, 2> last_P_alpha_;
    xt::xtensor<double, 2> last_P_beta_;

    // Shared-pointer storage so observers can capture handles that stay valid
    // after the driver copies this evaluator into its owned-component list.
    std::shared_ptr<bool> has_last_result_{std::make_shared<bool>(false)};
    std::shared_ptr<SCFResult> last_result_{std::make_shared<SCFResult>()};
    std::shared_ptr<IntegralPack> last_integrals_{std::make_shared<IntegralPack>()};
    std::shared_ptr<SCFConfig> config_shared_;  // initialized in ctor

    std::function<void(const SCFResult&)> on_evaluate_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
