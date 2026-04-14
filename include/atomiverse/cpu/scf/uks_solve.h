/// @file cpu/scf/uks_solve.h
/// @brief CPU UKS solver -- unrestricted Kohn-Sham DFT reference implementation.
///
/// Phase 1: CPU-only, supporting LSDA (Slater + VWN5) and PBE GGA functionals.
/// Reuses IntegralPack (S, T, V, ERI) and canonical orthogonalization from RHF/UHF.
/// Replaces HF exchange with XC functionals evaluated on atom-centered numerical grids.
///
/// @see cpu::scf::uhf_solve for the UHF solver (same spin structure, HF exchange)
/// @see cpu::scf::xc_functional.h for XC functional evaluator
/// @see cpu::scf::grid.h for molecular integration grid
#pragma once

#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/xc_functional.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief UKS solver result with spin-resolved quantities and XC decomposition.
struct UKSResult {
    bool converged = false;          ///< Whether SCF converged within max_iter.
    int iterations = 0;              ///< Number of iterations performed.
    double energy = 0;               ///< Total energy E_KS [Hartree].
    double electronic_energy = 0;    ///< Electronic energy only [Hartree].
    double one_electron_energy = 0;  ///< sum_sigma Tr[Hcore * P^sigma] [Hartree].
    double coulomb_energy = 0;       ///< 0.5 Tr[J * P] [Hartree].
    double xc_energy = 0;           ///< E_xc from grid integration [Hartree].
    double residual = 0;             ///< max(||dP^alpha||_F, ||dP^beta||_F).
    int n_alpha = 0;                 ///< Number of alpha electrons.
    int n_beta = 0;                  ///< Number of beta electrons.
    double integrated_n_alpha = 0;   ///< integral rho_alpha dtau (grid quality check).
    double integrated_n_beta = 0;    ///< integral rho_beta dtau.
    double spin_contamination = 0;   ///< <S^2> - S(S+1) (informational).

    xt::xtensor<double, 1> eigenvalues_alpha; ///< Alpha orbital energies [Hartree].
    xt::xtensor<double, 1> eigenvalues_beta;  ///< Beta orbital energies [Hartree].
    xt::xtensor<double, 2> P_alpha;  ///< Final alpha density matrix.
    xt::xtensor<double, 2> P_beta;   ///< Final beta density matrix.
    xt::xtensor<double, 2> C_alpha;  ///< Final alpha MO coefficient matrix.
    xt::xtensor<double, 2> C_beta;   ///< Final beta MO coefficient matrix.
    std::vector<IterationRecord> trace; ///< Per-iteration convergence trace.
};

/// @brief Run UKS on CPU with specified XC functional.
///
/// @param integrals  Molecular integrals (S, Hcore, ERI, V_nuc).
/// @param n_alpha    Number of alpha electrons.
/// @param n_beta     Number of beta electrons.
/// @param xc_type    XC functional (LSDA or PBE).
/// @param n_radial   Radial grid points per atom (default 50).
/// @param n_angular  Angular grid points (default 110).
/// @param tol        Convergence tolerance on density residual.
/// @param max_iter   Maximum SCF iterations.
/// @return UKSResult with energies, XC decomposition, densities, MO coefficients.
UKSResult uks_solve(const IntegralPack& integrals,
                    int n_alpha, int n_beta,
                    XCType xc_type = XCType::LSDA,
                    int n_radial = 50, int n_angular = 110,
                    double tol = 1e-8, int max_iter = 100);

} // namespace scf
} // namespace cpu
} // namespace ATV
