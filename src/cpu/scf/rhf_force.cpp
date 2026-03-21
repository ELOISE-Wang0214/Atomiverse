/// @file cpu/scf/rhf_force.cpp
/// @brief RHFForce implementation — BOMD and second-gen CPMD runner.

#include "atomiverse/cpu/scf/rhf_force.h"
#include "atomiverse/core/units.h"

#include <cmath>
#include <stdexcept>

namespace ATV {
namespace cpu {
namespace scf {

RHFForce::RHFForce(Real scf_tol, int max_scf_iter, int charge)
    : scf_tol_(scf_tol), max_scf_iter_(max_scf_iter), charge_(charge) {}

RHFForce& RHFForce::predictor(int order, int corrector_steps,
                               Real corrector_tol, bool fallback) {
    use_predictor_ = true;
    corrector_steps_ = corrector_steps;
    corrector_tol_ = corrector_tol;
    fallback_full_scf_ = fallback;
    aspc_ = std::make_unique<ASPCPredictor>(order);
    return *this;
}

void RHFForce::compute_forces(Frame& frame) {
    int natm = static_cast<int>(frame.n_atoms());

    // Rebuild integrals at current geometry
    auto integrals = build_integrals_sto3g(frame, charge_);

    RHFResult result;
    bool used_predictor = false;

    if (use_predictor_ && aspc_) {
        xt::xtensor<Real, 2> P_guess;
        if (aspc_->predict(P_guess)) {
            // Predictor succeeded -- run limited SCF with P_guess as initial density
            // For now, use the standard solver (it starts from Hcore guess anyway).
            // TODO: Pass P_guess as initial density to rhf_solve when supported.
            result = rhf_solve(integrals, corrector_tol_, corrector_steps_);
            used_predictor = true;

            if (!result.converged && fallback_full_scf_) {
                result = rhf_solve(integrals, scf_tol_, max_scf_iter_);
                used_predictor = false;
            }

            // Apply corrector blending
            xt::xtensor<Real, 2> P_corrected;
            aspc_->correct(result.P, P_guess, P_corrected);
            aspc_->record(P_corrected);
        } else {
            // Not enough history -- full SCF
            result = rhf_solve(integrals, scf_tol_, max_scf_iter_);
            aspc_->record(result.P);
        }
    } else {
        // BOMD mode: full SCF every step
        result = rhf_solve(integrals, scf_tol_, max_scf_iter_);
        if (use_predictor_ && aspc_)
            aspc_->record(result.P);
    }

    last_result_ = result;
    has_last_result_ = true;

    // Compute analytic gradient
    auto grad = rhf_gradient(frame, integrals, result);

    // Write forces and energy to frame
    for (int A = 0; A < natm; ++A) {
        frame.fx(A) = grad.forces(A, 0);
        frame.fy(A) = grad.forces(A, 1);
        frame.fz(A) = grad.forces(A, 2);
    }
    frame.energy = grad.energy;
}

void RHFForce::run(Frame& frame, int n_steps, Real dt) {
    int natm = static_cast<int>(frame.n_atoms());
    trace_.clear();
    trace_.reserve(n_steps);

    // Initial forces
    compute_forces(frame);

    for (int step = 0; step < n_steps; ++step) {
        // ── Velocity Verlet: half-step velocity ──
        for (int i = 0; i < natm; ++i) {
            Real ax = frame.fx(i) * kFtm2v / frame.mass(i);
            Real ay = frame.fy(i) * kFtm2v / frame.mass(i);
            Real az = frame.fz(i) * kFtm2v / frame.mass(i);
            frame.vx(i) += 0.5 * dt * ax;
            frame.vy(i) += 0.5 * dt * ay;
            frame.vz(i) += 0.5 * dt * az;
        }

        // ── Full-step position ──
        for (int i = 0; i < natm; ++i) {
            frame.x(i) += dt * frame.vx(i);
            frame.y(i) += dt * frame.vy(i);
            frame.z(i) += dt * frame.vz(i);
        }

        // ── Recompute forces at new positions ──
        compute_forces(frame);

        // ── Half-step velocity (second half) ──
        for (int i = 0; i < natm; ++i) {
            Real ax = frame.fx(i) * kFtm2v / frame.mass(i);
            Real ay = frame.fy(i) * kFtm2v / frame.mass(i);
            Real az = frame.fz(i) * kFtm2v / frame.mass(i);
            frame.vx(i) += 0.5 * dt * ax;
            frame.vy(i) += 0.5 * dt * ay;
            frame.vz(i) += 0.5 * dt * az;
        }

        // ── Compute kinetic energy and temperature ──
        Real ke = 0;
        for (int i = 0; i < natm; ++i)
            ke += 0.5 * kMvv2e * frame.mass(i)
               * (frame.vx(i)*frame.vx(i) + frame.vy(i)*frame.vy(i) + frame.vz(i)*frame.vz(i));
        frame.kinetic = ke;

        Real temp = 2.0 * ke / (3.0 * natm * kBoltzmann);

        // ── Record diagnostics ──
        AIMDStepRecord rec;
        rec.step = step;
        rec.potential_energy = frame.energy;
        rec.kinetic_energy = ke;
        rec.total_energy = ke + frame.energy;
        rec.temperature = temp;
        rec.scf_iterations = has_last_result_ ? last_result_.iterations : 0;
        rec.scf_residual = has_last_result_ ? last_result_.residual : 0;
        rec.used_predictor = use_predictor_ && aspc_ && aspc_->history_size() > aspc_->order();
        trace_.push_back(rec);

        if (on_step_)
            on_step_(rec);
    }
}

} // namespace scf
} // namespace cpu
} // namespace ATV
