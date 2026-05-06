/// @file cpu/scf/scf_evaluator.cpp
/// @brief SCFEvaluator implementation — unified AIMD force component.
///
/// As of driver-taxonomy phase 4, the SCF cycle runs through a
/// pre-built `cpu::scf::SCFSolver` member: the orthogonalization matrix,
/// shell-pair lists, and KS grid are allocated once in `build()` and
/// reused across every `launch()`. The ASPC corrector path uses a
/// transient solver because it needs a different (tighter `tol`,
/// smaller `max_iter`) configuration; that's the rare warm-start path.

#include "atomiverse/cpu/scf/scf_evaluator.h"

#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"

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

void SCFEvaluator::build(const cpu::BuildContext& ctx) {
    // Build a representative IntegralPack for the current geometry — used
    // both as the working pack for launch #1 and to size the SCFSolver's
    // persistent workspaces (X, shell-pair lists, KS grid). Subsequent
    // launches reuse those workspaces.
    *last_integrals_ = build_integrals(ctx.frame, basis_, charge_);

    solver_ = std::make_unique<SCFSolver>(config_);
    solver_->build(*last_integrals_);

    ATV_INFO(
        "[scf_eval] build charge={} predictor={} compute_force={} basis={} n_basis={}",
        charge_, use_predictor_, compute_force_,
        basis_set_name(basis_), last_integrals_->n_basis);
}

void SCFEvaluator::launch(const cpu::LaunchContext& ctx) {
    const auto t_total = LogClock::now();
    auto& frame = ctx.frame;
    const int natm = static_cast<int>(frame.atoms.size());
    ATV_INFO(
        "[scf_eval] step={} start natm={} charge={}",
        ctx.step, natm, charge_);

    const auto t_integrals = LogClock::now();
    *last_integrals_ = build_integrals(frame, basis_, charge_);
    auto& integrals = *last_integrals_;
    ATV_DEBUG(
        "[scf_eval] step={} integrals_ms={:.1f} basis={} n_basis={} electrons={}",
        ctx.step, elapsed_ms(t_integrals),
        basis_set_name(basis_), integrals.n_basis, integrals.n_electrons);

    SCFResult result;

    const auto t_scf = LogClock::now();
    SCFInitialGuess density_guess;
    if (use_predictor_ && aspc_alpha_ && aspc_beta_) {
        xt::xtensor<double, 2> Pa_guess, Pb_guess;
        const bool have_alpha = aspc_alpha_->predict(Pa_guess);
        const bool have_beta  = aspc_beta_->predict(Pb_guess);

        if (have_alpha && have_beta) {
            density_guess.P_alpha = std::move(Pa_guess);
            density_guess.P_beta  = std::move(Pb_guess);
            density_guess.has_alpha = true;
            density_guess.has_beta  = true;

            // ASPC corrector run: tighter tol, fewer iterations, transient
            // solver because the SCFConfig differs from the main one.
            SCFConfig corrector_config = config_;
            corrector_config.tol      = corrector_tol_;
            corrector_config.max_iter = corrector_steps_;
            ATV_INFO("[scf_eval] step={} scf guess source=aspc_predictor", ctx.step);
            {
                SCFSolver corrector{corrector_config};
                corrector.build(integrals);
                result = corrector.run(integrals, &density_guess);
            }

            if (!result.converged && fallback_full_scf_)
                result = solver_->run(integrals, &density_guess);

            xt::xtensor<double, 2> Pa_corr, Pb_corr;
            aspc_alpha_->correct(result.P_alpha, density_guess.P_alpha, Pa_corr);
            aspc_beta_->correct(result.P_beta,  density_guess.P_beta,  Pb_corr);
            aspc_alpha_->record(Pa_corr);
            aspc_beta_->record(Pb_corr);
            has_last_density_ = true;
            last_P_alpha_ = result.P_alpha;
            last_P_beta_  = result.P_beta;
        } else {
            if (has_last_density_) {
                density_guess.P_alpha = last_P_alpha_;
                density_guess.P_beta  = last_P_beta_;
                density_guess.has_alpha = true;
                density_guess.has_beta  = true;
                ATV_INFO("[scf_eval] step={} scf guess source=previous_density", ctx.step);
                result = solver_->run(integrals, &density_guess);
            } else {
                result = solver_->run(integrals);
            }
            aspc_alpha_->record(result.P_alpha);
            aspc_beta_->record(result.P_beta);
            has_last_density_ = true;
            last_P_alpha_ = result.P_alpha;
            last_P_beta_  = result.P_beta;
        }
    } else if (has_last_density_) {
        density_guess.P_alpha = last_P_alpha_;
        density_guess.P_beta  = last_P_beta_;
        density_guess.has_alpha = true;
        density_guess.has_beta  = true;
        ATV_INFO("[scf_eval] step={} scf guess source=previous_density", ctx.step);
        result = solver_->run(integrals, &density_guess);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_  = result.P_beta;
    } else {
        result = solver_->run(integrals);
        has_last_density_ = true;
        last_P_alpha_ = result.P_alpha;
        last_P_beta_  = result.P_beta;
    }
    ATV_DEBUG(
        "[scf_eval] step={} scf_ms={:.1f} converged={} iterations={} energy_ha={:.10f}",
        ctx.step, elapsed_ms(t_scf), result.converged, result.iterations, result.energy);

    // Forces — only when gradient is enabled. Energy-only single-point
    // mode (with_gradient(false)) skips this entirely.
    if (compute_force_) {
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
    } else {
        // Energy-only path. Convert result.energy (Ha) to frame.energy (eV).
        frame.energy = result.energy * ATV::kHartree2eV;
        ATV_INFO(
            "[scf_eval] step={} done total_ms={:.1f} energy_eV={:.6f} (energy-only) scf_converged={} scf_iterations={}",
            ctx.step, elapsed_ms(t_total), frame.energy,
            result.converged, result.iterations);
    }

    if (on_evaluate_)
        on_evaluate_(result);

    *last_result_ = std::move(result);
    *has_last_result_ = true;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
