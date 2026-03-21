/// @file cpu/scf/uhf_force.cpp
/// @brief UHFForce implementation — open-shell AIMD runner.

#include "atomiverse/cpu/scf/uhf_force.h"
#include "atomiverse/core/units.h"

#include <cmath>

namespace ATV {
namespace cpu {
namespace scf {

UHFForce::UHFForce(int n_alpha, int n_beta, Real scf_tol, int max_scf_iter, int charge)
    : n_alpha_(n_alpha), n_beta_(n_beta),
      scf_tol_(scf_tol), max_scf_iter_(max_scf_iter), charge_(charge) {}

UHFForce& UHFForce::predictor(int order, int corrector_steps,
                               Real corrector_tol, bool fallback) {
    use_predictor_ = true;
    corrector_steps_ = corrector_steps;
    corrector_tol_ = corrector_tol;
    fallback_full_scf_ = fallback;
    aspc_alpha_ = std::make_unique<ASPCPredictor>(order);
    aspc_beta_ = std::make_unique<ASPCPredictor>(order);
    return *this;
}

void UHFForce::compute_forces(Frame& frame) {
    int natm = static_cast<int>(frame.n_atoms());

    auto integrals = build_integrals_sto3g(frame, charge_);

    UHFResult result;

    if (use_predictor_ && aspc_alpha_ && aspc_beta_) {
        xt::xtensor<Real, 2> Pa_guess, Pb_guess;
        bool have_alpha = aspc_alpha_->predict(Pa_guess);
        bool have_beta = aspc_beta_->predict(Pb_guess);

        if (have_alpha && have_beta) {
            result = uhf_solve(integrals, n_alpha_, n_beta_, corrector_tol_, corrector_steps_);

            if (!result.converged && fallback_full_scf_)
                result = uhf_solve(integrals, n_alpha_, n_beta_, scf_tol_, max_scf_iter_);

            xt::xtensor<Real, 2> Pa_corr, Pb_corr;
            aspc_alpha_->correct(result.P_alpha, Pa_guess, Pa_corr);
            aspc_beta_->correct(result.P_beta, Pb_guess, Pb_corr);
            aspc_alpha_->record(Pa_corr);
            aspc_beta_->record(Pb_corr);
        } else {
            result = uhf_solve(integrals, n_alpha_, n_beta_, scf_tol_, max_scf_iter_);
            aspc_alpha_->record(result.P_alpha);
            aspc_beta_->record(result.P_beta);
        }
    } else {
        result = uhf_solve(integrals, n_alpha_, n_beta_, scf_tol_, max_scf_iter_);
    }

    last_result_ = result;
    has_last_result_ = true;

    auto grad = uhf_gradient(frame, integrals, result);

    for (int A = 0; A < natm; ++A) {
        frame.fx(A) = grad.forces(A, 0);
        frame.fy(A) = grad.forces(A, 1);
        frame.fz(A) = grad.forces(A, 2);
    }
    frame.energy = grad.energy;
}

void UHFForce::run(Frame& frame, int n_steps, Real dt) {
    int natm = static_cast<int>(frame.n_atoms());
    trace_.clear();
    trace_.reserve(n_steps);

    compute_forces(frame);

    for (int step = 0; step < n_steps; ++step) {
        for (int i = 0; i < natm; ++i) {
            Real ax = frame.fx(i) * kFtm2v / frame.mass(i);
            Real ay = frame.fy(i) * kFtm2v / frame.mass(i);
            Real az = frame.fz(i) * kFtm2v / frame.mass(i);
            frame.vx(i) += 0.5 * dt * ax;
            frame.vy(i) += 0.5 * dt * ay;
            frame.vz(i) += 0.5 * dt * az;
        }

        for (int i = 0; i < natm; ++i) {
            frame.x(i) += dt * frame.vx(i);
            frame.y(i) += dt * frame.vy(i);
            frame.z(i) += dt * frame.vz(i);
        }

        compute_forces(frame);

        for (int i = 0; i < natm; ++i) {
            Real ax = frame.fx(i) * kFtm2v / frame.mass(i);
            Real ay = frame.fy(i) * kFtm2v / frame.mass(i);
            Real az = frame.fz(i) * kFtm2v / frame.mass(i);
            frame.vx(i) += 0.5 * dt * ax;
            frame.vy(i) += 0.5 * dt * ay;
            frame.vz(i) += 0.5 * dt * az;
        }

        Real ke = 0;
        for (int i = 0; i < natm; ++i)
            ke += 0.5 * kMvv2e * frame.mass(i)
               * (frame.vx(i)*frame.vx(i) + frame.vy(i)*frame.vy(i) + frame.vz(i)*frame.vz(i));
        frame.kinetic = ke;

        Real temp = 2.0 * ke / (3.0 * natm * kBoltzmann);

        AIMDStepRecord rec;
        rec.step = step;
        rec.potential_energy = frame.energy;
        rec.kinetic_energy = ke;
        rec.total_energy = ke + frame.energy;
        rec.temperature = temp;
        rec.scf_iterations = has_last_result_ ? last_result_.iterations : 0;
        rec.scf_residual = has_last_result_ ? last_result_.residual : 0;
        rec.used_predictor = use_predictor_ && aspc_alpha_ && aspc_alpha_->history_size() > aspc_alpha_->order();
        trace_.push_back(rec);

        if (on_step_)
            on_step_(rec);
    }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
