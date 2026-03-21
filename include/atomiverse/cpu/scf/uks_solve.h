/// @file cpu/scf/uks_solve.h
/// @brief CPU UKS solver — unrestricted Kohn-Sham DFT reference implementation.
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
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
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
    Real energy = 0;                 ///< Total energy E_KS [Hartree].
    Real electronic_energy = 0;      ///< Electronic energy only [Hartree].
    Real one_electron_energy = 0;    ///< Σ_σ Tr[Hcore · P^σ] [Hartree].
    Real coulomb_energy = 0;         ///< ½ Tr[J · P] [Hartree].
    Real xc_energy = 0;             ///< E_xc from grid integration [Hartree].
    Real residual = 0;               ///< max(||dP^α||_F, ||dP^β||_F).
    int n_alpha = 0;                 ///< Number of alpha electrons.
    int n_beta = 0;                  ///< Number of beta electrons.
    Real integrated_n_alpha = 0;     ///< ∫ ρ_α dτ (grid quality check).
    Real integrated_n_beta = 0;      ///< ∫ ρ_β dτ.
    Real spin_contamination = 0;     ///< ⟨S²⟩ - S(S+1) (informational).

    xt::xtensor<Real, 1> eigenvalues_alpha; ///< Alpha orbital energies [Hartree].
    xt::xtensor<Real, 1> eigenvalues_beta;  ///< Beta orbital energies [Hartree].
    xt::xtensor<Real, 2> P_alpha;  ///< Final alpha density matrix.
    xt::xtensor<Real, 2> P_beta;   ///< Final beta density matrix.
    xt::xtensor<Real, 2> C_alpha;  ///< Final alpha MO coefficient matrix.
    xt::xtensor<Real, 2> C_beta;   ///< Final beta MO coefficient matrix.
    std::vector<IterationRecord> trace; ///< Per-iteration convergence trace.
};

/// @brief Run UKS on CPU with specified XC functional.
///
/// Algorithm:
///   1. Build molecular grid (Becke partitioning + Lebedev angular)
///   2. X = S^{-1/2} (canonical orthogonalization)
///   3. Initial guess: diagonalize Hcore (+ HOMO-LUMO mixing for open-shell)
///   4. Damped SCF loop:
///      - Build J(P) from ERI (Coulomb only, no HF exchange)
///      - Evaluate XC on grid → E_xc, V^α_xc, V^β_xc
///      - F^σ = Hcore + J + V^σ_xc
///      - Diagonalize, rebuild densities, check convergence
///   5. Final energy: E = Σ_σ Tr[Hcore · P^σ] + ½ Tr[J · P] + E_xc + V_nuc
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
                    Real tol = Real(1e-8), int max_iter = 100);

} // namespace scf
} // namespace cpu
} // namespace ATV
