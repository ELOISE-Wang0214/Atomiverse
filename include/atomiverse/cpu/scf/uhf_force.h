/// @file cpu/scf/uhf_force.h
/// @brief UHFForce — UHF force evaluator for ab initio molecular dynamics.
///
/// Open-shell AIMD force evaluator. Same interface as RHFForce but with
/// spin-resolved SCF (separate alpha/beta channels).
///
/// CPU-only implementation. Analogous to cpu::scf::RHFForce for closed-shell.
#pragma once

#include <functional>
#include <string>
#include <vector>
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"
#include "atomiverse/cpu/scf/rhf_force.h"  // AIMDStepRecord

namespace ATV {
namespace cpu {
namespace scf {

/// @brief UHF force evaluator for open-shell AIMD.
class UHFForce {
public:
    /// @brief Construct with SCF parameters.
    /// @param n_alpha  Number of alpha electrons.
    /// @param n_beta   Number of beta electrons.
    /// @param scf_tol  SCF convergence tolerance [Hartree].
    /// @param max_scf_iter Maximum SCF iterations.
    /// @param charge   Molecular charge (0 = neutral).
    UHFForce(int n_alpha, int n_beta,
             Real scf_tol = Real(1e-8), int max_scf_iter = 100, int charge = 0);

    /// @brief Enable ASPC predictor-corrector mode.
    /// @param order           Extrapolation order K (default 3).
    /// @param corrector_steps Max SCF corrector iterations (default 2).
    /// @param corrector_tol   Corrector convergence tolerance (default 1e-4).
    /// @param fallback        Fall back to full SCF if corrector fails.
    /// @return Reference to this for chaining.
    UHFForce& predictor(int order = 3, int corrector_steps = 2,
                        Real corrector_tol = Real(1e-4), bool fallback = true);

    /// @brief Run AIMD with velocity Verlet integration.
    ///
    /// Executes the full BOMD/CPMD2 loop with UHF forces:
    ///   1. Initial SCF + gradient
    ///   2. For each step: half-step v, full-step x, SCF+gradient, half-step v
    ///
    /// @param frame   CPU Frame (modified in-place with final positions/velocities/forces).
    /// @param n_steps Number of MD steps.
    /// @param dt      Time step [fs].
    void run(Frame& frame, int n_steps, Real dt);

    /// @brief Compute forces at the current geometry (single-point).
    ///
    /// Writes forces to frame.fx/fy/fz and energy to frame.energy.
    /// @param frame CPU Frame.
    void compute_forces(Frame& frame);

    /// @brief Access the step-by-step diagnostic trace.
    /// @return Vector of per-step diagnostic records.
    const std::vector<AIMDStepRecord>& trace() const { return trace_; }

    /// @brief Set a per-step callback for live monitoring.
    /// @param cb Callback invoked after each MD step with that step's diagnostics.
    void on_step(std::function<void(const AIMDStepRecord&)> cb) { on_step_ = std::move(cb); }

private:
    int n_alpha_, n_beta_;
    Real scf_tol_;
    int max_scf_iter_;
    int charge_;

    bool use_predictor_ = false;
    int corrector_steps_ = 2;
    Real corrector_tol_ = Real(1e-4);
    bool fallback_full_scf_ = true;
    std::unique_ptr<ASPCPredictor> aspc_alpha_;
    std::unique_ptr<ASPCPredictor> aspc_beta_;

    UHFResult last_result_;
    bool has_last_result_ = false;

    std::vector<AIMDStepRecord> trace_;
    std::function<void(const AIMDStepRecord&)> on_step_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
