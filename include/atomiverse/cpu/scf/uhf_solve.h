/// @file cpu/scf/uhf_solve.h
/// @brief CPU UHF solver — unrestricted Hartree-Fock reference implementation.
///
/// Full Unrestricted Hartree-Fock for open-shell and closed-shell systems.
/// Maintains separate Fock matrices, densities, and DIIS histories for
/// alpha and beta spin channels.
///
/// @see cpu::scf::rhf_solve for the closed-shell restricted solver
/// @see cpu::scf::IntegralPack for the input integrals
#pragma once

#include <vector>
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief UHF solver result with spin-resolved quantities.
struct UHFResult {
    bool converged = false;          ///< Whether SCF converged within max_iter.
    int iterations = 0;              ///< Number of iterations performed.
    Real energy = 0;                 ///< Total energy (electronic + nuclear) [Hartree].
    Real electronic_energy = 0;      ///< Electronic energy only [Hartree].
    Real residual = 0;               ///< max(||dP_alpha||_F, ||dP_beta||_F).
    int n_alpha = 0;                 ///< Number of alpha electrons.
    int n_beta = 0;                  ///< Number of beta electrons.
    Real spin_contamination = 0;     ///< <S^2> - S*(S+1): 0 for pure spin states.

    xt::xtensor<Real, 1> eigenvalues_alpha; ///< Alpha orbital energies [Hartree].
    xt::xtensor<Real, 1> eigenvalues_beta;  ///< Beta orbital energies [Hartree].
    xt::xtensor<Real, 2> P_alpha;  ///< Final alpha density matrix.
    xt::xtensor<Real, 2> P_beta;   ///< Final beta density matrix.
    xt::xtensor<Real, 2> C_alpha;  ///< Final alpha MO coefficient matrix.
    xt::xtensor<Real, 2> C_beta;   ///< Final beta MO coefficient matrix.
    std::vector<IterationRecord> trace; ///< Per-iteration convergence trace.
};

/// @brief Run UHF on CPU using xtensor-blas with separate alpha/beta DIIS.
///
/// Algorithm:
///   1. X = S^{-1/2} (canonical orthogonalization, shared by both channels)
///   2. Initial guess: diagonalize Hcore for both channels
///   3. Iterate:
///      - P_total = P_alpha + P_beta
///      - G_alpha_ij = sum_kl P_total_kl*(ij|kl) - P_alpha_kl*(ik|jl)
///      - G_beta_ij  = sum_kl P_total_kl*(ij|kl) - P_beta_kl*(ik|jl)
///      - F_alpha = Hcore + G_alpha, F_beta = Hcore + G_beta
///      - DIIS extrapolation (independent for alpha/beta)
///      - Diagonalize each channel in orthogonal basis
///      - Rebuild P_alpha = C_alpha_occ * C_alpha_occ^T
///      - Rebuild P_beta  = C_beta_occ  * C_beta_occ^T
///      - E = 0.5 * Tr[(Hcore+F_alpha)*P_alpha + (Hcore+F_beta)*P_beta] + V_nuc
///      - Check max(||dP_alpha||_F, ||dP_beta||_F) < tol
///   4. Compute spin contamination <S^2> - S*(S+1)
///
/// @param integrals Molecular integrals (S, Hcore, ERI, V_nuc).
/// @param n_alpha Number of alpha electrons.
/// @param n_beta  Number of beta electrons (0 for pure alpha systems).
/// @param tol     Convergence tolerance on density matrix Frobenius norm.
/// @param max_iter Maximum SCF iterations.
/// @return UHFResult with energies, spin diagnostics, densities, MO coefficients.
UHFResult uhf_solve(const IntegralPack& integrals, int n_alpha, int n_beta,
                    Real tol = Real(1e-8), int max_iter = 100);

} // namespace scf
} // namespace cpu
} // namespace ATV
