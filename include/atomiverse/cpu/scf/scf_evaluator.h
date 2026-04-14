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
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

namespace ATV {
namespace cpu {
namespace scf {

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

    void build(const cpu::BuildContext& ctx) override;
    void launch(const cpu::LaunchContext& ctx) override;

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

    std::function<void(const SCFResult&)> on_evaluate_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
