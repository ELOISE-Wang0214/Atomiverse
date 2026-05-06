/// @file cpu/scf/scf_solver.h
/// @brief cpu::scf::SCFSolver — public low-level inner SCF solver.
///
/// **NOT a `Driver<...>` subclass** — `SCFSolver` is the *inner* solver
/// (the SCF self-consistency loop) that lives behind force-providing
/// evaluator components. The user-facing way to compute single-point
/// electronic structure is `cpu::Evaluate + cpu::scf::SCFEvaluator`;
/// `SCFSolver` is the internal machinery `SCFEvaluator` calls each step.
///
/// `SCFSolver` is also a public API for low-level use (debugging,
/// testing, density-matrix dumps, DIIS comparisons, initial-guess
/// experiments). It owns five `internal::*` step helpers, each with
/// build/step lifecycle and pre-allocated workspaces. The helpers are
/// **not** registered into any outer `Driver` Plan.
///
/// Usage (low-level direct SCF):
/// @code
///   cpu::scf::SCFSolver solver{cpu::scf::default_config(integrals)};
///   solver.build(integrals);
///   const auto& result = solver.run(integrals);
///   // result.energy, result.P_alpha, result.eigenvalues_alpha, ...
/// @endcode
///
/// Usage (single-point through the outer facade):
/// @code
///   cpu::Evaluate eval;
///   eval.add(cpu::scf::SCFEvaluator{config});  // owns an SCFSolver internally
///   eval.build();
///   eval.run(frame);
/// @endcode
#pragma once

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/internal/density_builder.h"
#include "atomiverse/cpu/scf/internal/diis_mixer.h"
#include "atomiverse/cpu/scf/internal/dual_convergence_check.h"
#include "atomiverse/cpu/scf/internal/eigen_solver.h"
#include "atomiverse/cpu/scf/internal/hamiltonian_builder.h"
#include "atomiverse/cpu/scf/scf_solve.h"   // SCFConfig, SCFInitialGuess, SCFResult

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Inner SCF self-consistency loop. Public low-level surface.
///
/// Construction → build(reference) → run(integrals[, guess]).
/// `last_result()` is valid after each `run()` call.
class SCFSolver {
public:
    /// @brief Construct from an SCF configuration.
    /// @param config  Method, spin, convergence, grid parameters.
    explicit SCFSolver(SCFConfig config);

    /// @brief One-time setup: builds shell-pair lists, KS grid (if KS),
    /// and the orthogonalization matrix X. Sized to `reference.n_basis`.
    /// Subsequent `run()` calls assume the same n_basis (geometry may
    /// change in-place; n_basis is fixed by the basis set).
    /// @param reference IntegralPack used to size internal workspaces.
    void build(const IntegralPack& reference);

    /// @brief Execute the SCF cycle on `integrals`.
    /// @param integrals  Molecular integrals at the current geometry.
    /// @param guess      Optional density-matrix initial guess.
    /// @return Reference to the owned `SCFResult` (valid until the next
    ///         `run()` or `SCFSolver` destruction).
    SCFResult& run(const IntegralPack& integrals,
                   const SCFInitialGuess* guess = nullptr);

    /// @brief Access the most recent SCF result.
    const SCFResult& last_result() const { return result_; }

    /// @brief Access the current SCF configuration (for inspection).
    const SCFConfig& config() const { return config_; }

private:
    SCFConfig config_;
    bool      built_ = false;

    // Per-build state.
    RealMatrix X_;  ///< Orthogonalization matrix (S^{-1/2}). Owned.

    // Internal SCF step helpers.
    internal::HamiltonianBuilder hamiltonian_;
    internal::EigenSolver        eigen_;
    internal::DensityBuilder     density_;
    internal::DIISMixer          diis_;
    internal::DualConvergenceCheck convergence_;

    // Owned result buffer (returned by reference from run()).
    SCFResult result_;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
