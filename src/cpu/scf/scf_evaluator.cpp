/// @file cpu/scf/scf_evaluator.cpp
/// @brief SCFEvaluator implementation — unified AIMD force component.

#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/core/log.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace scf {

SCFEvaluator::SCFEvaluator(SCFConfig config, BasisSet basis, int charge)
    : config_(std::move(config)), basis_(basis), charge_(charge),
      config_shared_(std::make_shared<SCFConfig>(config_)) {}

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
    ATV_INFO(
        "[scf_eval] build charge={} predictor={}",
        charge_, use_predictor_);
}

void SCFEvaluator::launch(const cpu::LaunchContext& ctx) {
    const auto t_total = LogClock::now();
    auto& frame = ctx.frame;
    int natm = static_cast<int>(frame.atoms.size());
    ATV_INFO(
        "[scf_eval] step={} start natm={} charge={}",
        ctx.step, natm, charge_);

    const auto t_integrals = LogClock::now();
    auto integrals = build_integrals(frame, basis_, charge_);
    ATV_DEBUG(
        "[scf_eval] step={} integrals_ms={:.1f} basis={} n_basis={} electrons={}",
        ctx.step, elapsed_ms(t_integrals),
        basis_set_name(basis_), integrals.n_basis, integrals.n_electrons);

    SCFResult result;

    const auto t_scf = LogClock::now();
    SCFInitialGuess density_guess;
    if (use_predictor_ && aspc_alpha_ && aspc_beta_) {
        xt::xtensor<double, 2> Pa_guess, Pb_guess;
        bool have_alpha = aspc_alpha_->predict(Pa_guess);
        bool have_beta = aspc_beta_->predict(Pb_guess);

        if (have_alpha && have_beta) {
            density_guess.P_alpha = std::move(Pa_guess);
            density_guess.P_beta = std::move(Pb_guess);
            density_guess.has_alpha = true;
            density_guess.has_beta = true;
            auto corrector_config = config_;
            corrector_config.tol = corrector_tol_;
            corrector_config.max_iter = corrector_steps_;
            ATV_INFO("[scf_eval] step={} scf guess source=aspc_predictor", ctx.step);
            result = scf_solve(integrals, corrector_config, &density_guess);

            if (!result.converged && fallback_full_scf_)
                result = scf_solve(integrals, config_, &density_guess);

            xt::xtensor<double, 2> Pa_corr, Pb_corr;
            aspc_alpha_->correct(result.P_alpha, density_guess.P_alpha, Pa_corr);
            aspc_beta_->correct(result.P_beta, density_guess.P_beta, Pb_corr);
            aspc_alpha_->record(Pa_corr);
            aspc_beta_->record(Pb_corr);
            has_last_density_ = true;
            last_P_alpha_ = result.P_alpha;
            last_P_beta_ = result.P_beta;
        } else {
            if (has_last_density_) {
                density_guess.P_alpha = last_P_alpha_;
                density_guess.P_beta = last_P_beta_;
                density_guess.has_alpha = true;
                density_guess.has_beta = true;
                ATV_INFO("[scf_eval] step={} scf guess source=previous_density", ctx.step);
                result = scf_solve(integrals, config_, &density_guess);
            } else {
                result = scf_solve(integrals, config_);
            }
            aspc_alpha_->record(result.P_alpha);
            aspc_beta_->record(result.P_beta);
            has_last_density_ = true;
            last_P_alpha_ = result.P_alpha;
            last_P_beta_ = result.P_beta;
        }
    } else if (has_last_density_) {
        density_guess.P_alpha = last_P_alpha_;
        density_guess.P_beta = last_P_beta_;
        density_guess.has_alpha = true;
        density_guess.has_beta = true;
        ATV_INFO("[scf_eval] step={} scf guess source=previous_density", ctx.step);
        result = scf_solve(integrals, config_, &density_guess);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    } else {
        result = scf_solve(integrals, config_);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_ = result.P_beta;
    }
    ATV_DEBUG(
        "[scf_eval] step={} scf_ms={:.1f} converged={} iterations={} energy_ha={:.10f}",
        ctx.step, elapsed_ms(t_scf), result.converged, result.iterations, result.energy);

    const auto t_gradient = LogClock::now();
    auto grad = scf_gradient(frame, integrals, result, config_);
    ATV_DEBUG(
        "[scf_eval] step={} gradient_ms={:.1f}",
        ctx.step, elapsed_ms(t_gradient));

    double max_force = 0.0;
    for (int A = 0; A < natm; ++A) {
        frame.atoms.fx(A) = grad.forces(A, 0);
        frame.atoms.fy(A) = grad.forces(A, 1);
        frame.atoms.fz(A) = grad.forces(A, 2);
        const double fx = grad.forces(A, 0);
        const double fy = grad.forces(A, 1);
        const double fz = grad.forces(A, 2);
        const double mag = std::sqrt(fx * fx + fy * fy + fz * fz);
        max_force = std::max(max_force, mag);
    }
    frame.energy = grad.energy;

    ATV_INFO(
        "[scf_eval] step={} done total_ms={:.1f} energy_eV={:.6f} max_force_eVA={:.6f} scf_converged={} scf_iterations={}",
        ctx.step, elapsed_ms(t_total), frame.energy, max_force,
        result.converged, result.iterations);

    if (on_evaluate_)
        on_evaluate_(result);

    // Cache SCF state for downstream observers (OrbitalInfoObserver, cube
    // writers). Shell pointers inside `integrals.shells` index into
    // `integrals.shell_data`; moving both together keeps those pointers valid
    // because std::vector move preserves the underlying buffer addresses.
    *last_integrals_ = std::move(integrals);
    *last_result_ = std::move(result);
    *has_last_result_ = true;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
