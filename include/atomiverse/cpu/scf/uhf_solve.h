/// @file cpu/scf/uhf_solve.h
/// @brief CPU UHF solver -- unrestricted Hartree-Fock reference implementation.
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

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief UHF solver result with spin-resolved quantities.
struct UHFResult {
    bool converged = false;          ///< Whether SCF converged within max_iter.
    int iterations = 0;              ///< Number of iterations performed.
    double energy = 0;               ///< Total energy (electronic + nuclear) [Hartree].
    double electronic_energy = 0;    ///< Electronic energy only [Hartree].
    double residual = 0;             ///< max(||dP_alpha||_F, ||dP_beta||_F).
    int n_alpha = 0;                 ///< Number of alpha electrons.
    int n_beta = 0;                  ///< Number of beta electrons.
    double spin_contamination = 0;   ///< <S^2> - S*(S+1): 0 for pure spin states.

    xt::xtensor<double, 1> eigenvalues_alpha; ///< Alpha orbital energies [Hartree].
    xt::xtensor<double, 1> eigenvalues_beta;  ///< Beta orbital energies [Hartree].
    xt::xtensor<double, 2> P_alpha;  ///< Final alpha density matrix.
    xt::xtensor<double, 2> P_beta;   ///< Final beta density matrix.
    xt::xtensor<double, 2> C_alpha;  ///< Final alpha MO coefficient matrix.
    xt::xtensor<double, 2> C_beta;   ///< Final beta MO coefficient matrix.
    std::vector<IterationRecord> trace; ///< Per-iteration convergence trace.
};

/// @brief Run UHF on CPU using xtensor-blas with separate alpha/beta DIIS.
///
/// @param integrals Molecular integrals (S, Hcore, ERI, V_nuc).
/// @param n_alpha Number of alpha electrons.
/// @param n_beta  Number of beta electrons (0 for pure alpha systems).
/// @param tol     Convergence tolerance on density matrix Frobenius norm.
/// @param max_iter Maximum SCF iterations.
/// @return UHFResult with energies, spin diagnostics, densities, MO coefficients.
UHFResult uhf_solve(const IntegralPack& integrals, int n_alpha, int n_beta,
                    double tol = 1e-8, int max_iter = 100);

} // namespace scf
} // namespace cpu
} // namespace ATV
