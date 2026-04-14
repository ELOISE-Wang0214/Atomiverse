/// @file cpu/scf/rhf_solve.h
/// @brief CPU RHF solver using xtensor-blas -- reference implementation.
///
/// Full Restricted Hartree-Fock for closed-shell systems.
/// Uses xt::linalg for eigendecomposition and matrix operations.
///
/// @see scf::SCF for the CUDA Driver-based implementation
/// @see cpu::scf::IntegralPack for the input integrals
#pragma once

#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Per-iteration data for convergence analysis.
struct IterationRecord {
    int iter = 0;              ///< Iteration number (1-based).
    double energy = 0;         ///< Total energy at this iteration [Hartree].
    double residual = 0;       ///< ||P_new - P_old||_F.
    bool diis_active = false;  ///< Whether DIIS was used this iteration.
};

/// @brief RHF solver result.
struct RHFResult {
    bool converged = false;              ///< Whether SCF converged within max_iter.
    int iterations = 0;                  ///< Number of iterations performed.
    double energy = 0;                   ///< Total energy (electronic + nuclear) [Hartree].
    double electronic_energy = 0;        ///< Electronic energy only [Hartree].
    double residual = 0;                 ///< Final ||P_new - P_old||_F.
    xt::xtensor<double, 1> eigenvalues;  ///< Orbital energies [Hartree].
    xt::xtensor<double, 2> P;           ///< Final density matrix.
    xt::xtensor<double, 2> C;           ///< Final MO coefficient matrix.
    std::vector<IterationRecord> trace;  ///< Per-iteration convergence trace.
};

/// @brief Run RHF on CPU using xtensor-blas.
///
/// Algorithm:
///   1. X = S^{-1/2} (canonical orthogonalization)
///   2. Iterate:
///      - F = Hcore + G(P) (Fock matrix)
///      - F' = X^T F X (orthogonal basis)
///      - Diagonalize F' -> (epsilon, C')
///      - C = X C' (back-transform)
///      - P = 2 C_occ C_occ^T (closed-shell density)
///      - E = 0.5 Tr[(Hcore + F) P] + V_nuc
///      - Check ||P - P_old||_F < tol
///
/// @param integrals Molecular integrals (S, Hcore, ERI, V_nuc).
/// @param tol Convergence tolerance.
/// @param max_iter Maximum SCF iterations.
/// @return RHFResult with energy, eigenvalues, density, MO coefficients.
RHFResult rhf_solve(const IntegralPack& integrals, double tol, int max_iter);

} // namespace scf
} // namespace cpu
} // namespace ATV
