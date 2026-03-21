/// @file cpu/scf/rhf_force.h
/// @brief RHFForce — RHF force evaluator for ab initio molecular dynamics.
///
/// Predefined force component that computes RHF energy and analytic nuclear
/// gradients at each MD step. Supports full SCF (BOMD mode) and ASPC
/// predictor-corrector (second-generation CPMD mode).
///
/// CPU-only implementation. Analogous to md::LJ for classical force fields.
///
/// Usage (BOMD):
/// @code{.cpp}
/// cpu::scf::RHFForce rhf{1e-8, 100};
/// rhf.run(frame, 100, 0.5);  // 100 steps, 0.5 fs
/// @endcode
///
/// Usage (CPMD2):
/// @code{.cpp}
/// cpu::scf::RHFForce rhf{1e-8, 100};
/// rhf.predictor(3, 2, 1e-4);
/// rhf.run(frame, 1000, 0.5);
/// @endcode
#pragma once

#include <functional>
#include <string>
#include <vector>
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/gradient.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Per-step AIMD diagnostic record.
struct AIMDStepRecord {
    int step;               ///< MD step index (0-based).
    Real total_energy;      ///< KE + PE [eV].
    Real potential_energy;  ///< SCF energy [eV].
    Real kinetic_energy;    ///< Nuclear KE [eV].
    Real temperature;       ///< Instantaneous temperature [K].
    int scf_iterations;     ///< SCF iterations this step.
    Real scf_residual;      ///< Final SCF residual.
    bool used_predictor;    ///< Whether ASPC predictor was used.
};

/// @brief RHF force evaluator for AIMD.
///
/// Encapsulates: integral rebuild -> SCF solve -> gradient computation.
/// Manages density matrix history for ASPC predictor-corrector.
class RHFForce {
public:
    /// @brief Construct with SCF parameters.
    /// @param scf_tol      SCF convergence tolerance [Hartree].
    /// @param max_scf_iter Maximum SCF iterations.
    /// @param charge       Molecular charge (0 = neutral).
    RHFForce(Real scf_tol = Real(1e-8), int max_scf_iter = 100, int charge = 0);

    /// @brief Enable ASPC predictor-corrector mode.
    /// @param order           Extrapolation order K (default 3).
    /// @param corrector_steps Max SCF corrector iterations (default 2).
    /// @param corrector_tol   Corrector convergence tolerance (default 1e-4).
    /// @param fallback        Fall back to full SCF if corrector fails.
    /// @return Reference to this for chaining.
    RHFForce& predictor(int order = 3, int corrector_steps = 2,
                        Real corrector_tol = Real(1e-4), bool fallback = true);

    /// @brief Run AIMD with velocity Verlet integration.
    ///
    /// Executes the full BOMD/CPMD2 loop:
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
    Real scf_tol_;
    int max_scf_iter_;
    int charge_;

    bool use_predictor_ = false;
    int corrector_steps_ = 2;
    Real corrector_tol_ = Real(1e-4);
    bool fallback_full_scf_ = true;
    std::unique_ptr<ASPCPredictor> aspc_;

    // Last converged SCF result (for predictor warm-start)
    RHFResult last_result_;
    bool has_last_result_ = false;

    std::vector<AIMDStepRecord> trace_;
    std::function<void(const AIMDStepRecord&)> on_step_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
