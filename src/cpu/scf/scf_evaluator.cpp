/// @file cpu/scf/scf_evaluator.cpp
/// @brief SCFEvaluator implementation — unified AIMD force component.

#include "atomiverse/cpu/scf/scf_evaluator.h"

namespace ATV {
namespace cpu {
namespace scf {

SCFEvaluator::SCFEvaluator(SCFConfig config, int charge)
    : config_(std::move(config)), charge_(charge) {}

SCFEvaluator& SCFEvaluator::predictor(int order, int corrector_steps,
                                       double corrector_tol, bool fallback) {
    use_predictor_ = true;
    corrector_steps_ = corrector_steps;
    corrector_tol_ = corrector_tol;
    fallback_full_scf_ = fallback;
    aspc_alpha_ = std::make_unique<ASPCPredictor>(order);
    aspc_beta_ = std::make_unique<ASPCPredictor>(order);
    return *this;
}

SCFEvaluator& SCFEvaluator::on_evaluate(std::function<void(const SCFResult&)> cb) {
    on_evaluate_ = std::move(cb);
    return *this;
}

void SCFEvaluator::build(const cpu::BuildContext& /*ctx*/) {
    // ASPC predictors already initialized in predictor() if enabled.
    // No additional build-time work needed.
}

void SCFEvaluator::launch(const cpu::LaunchContext& ctx) {
    auto& frame = ctx.frame;
    int natm = static_cast<int>(frame.atoms.size());

    auto integrals = build_integrals_sto3g(frame, charge_);

    SCFResult result;

    if (use_predictor_ && aspc_alpha_ && aspc_beta_) {
        xt::xtensor<double, 2> Pa_guess, Pb_guess;
        bool have_alpha = aspc_alpha_->predict(Pa_guess);
        bool have_beta = aspc_beta_->predict(Pb_guess);

        if (have_alpha && have_beta) {
            auto corrector_config = config_;
            corrector_config.tol = corrector_tol_;
            corrector_config.max_iter = corrector_steps_;
            result = scf_solve(integrals, corrector_config);

            if (!result.converged && fallback_full_scf_)
                result = scf_solve(integrals, config_);

            xt::xtensor<double, 2> Pa_corr, Pb_corr;
            aspc_alpha_->correct(result.P_alpha, Pa_guess, Pa_corr);
            aspc_beta_->correct(result.P_beta, Pb_guess, Pb_corr);
            aspc_alpha_->record(Pa_corr);
            aspc_beta_->record(Pb_corr);
        } else {
            result = scf_solve(integrals, config_);
            aspc_alpha_->record(result.P_alpha);
            aspc_beta_->record(result.P_beta);
        }
    } else {
        result = scf_solve(integrals, config_);
    }

    auto grad = scf_gradient(frame, integrals, result, config_);

    for (int A = 0; A < natm; ++A) {
        frame.atoms.fx(A) = grad.forces(A, 0);
        frame.atoms.fy(A) = grad.forces(A, 1);
        frame.atoms.fz(A) = grad.forces(A, 2);
    }
    frame.energy = grad.energy;

    if (on_evaluate_)
        on_evaluate_(result);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
