/// @file cpu/scf/gradient.h
/// @brief Analytic nuclear gradients for RHF and UHF.
///
/// Computes forces F_I = -dE/dR_I from converged SCF results and integral
/// derivatives. Output in eV/Angstrom (converted from Hartree/Bohr at the
/// boundary).
///
/// References:
///   Pulay, Mol. Phys. 17, 197 (1969) — Pulay forces
///   Szabo & Ostlund, Appendix C — RHF gradient derivation
///   Helgaker, Jorgensen, Olsen, Ch. 9 — molecular gradients
#pragma once

#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/uhf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Result of a gradient computation.
struct GradientResult {
    xt::xtensor<Real, 2> forces;  ///< Forces [n_atoms, 3] in eV/Angstrom.
    Real energy;                   ///< Total energy in eV.
};

/// @brief Compute RHF analytic nuclear gradient.
///
/// @param frame     System geometry (positions in Angstrom, type_id = Z).
/// @param integrals Precomputed integrals at current geometry (Hartree/Bohr).
/// @param result    Converged RHF result (P, C, eigenvalues).
/// @return Forces on all atoms (eV/Angstrom) and total energy (eV).
GradientResult rhf_gradient(
    const Frame& frame,
    const IntegralPack& integrals,
    const RHFResult& result
);

/// @brief Compute UHF analytic nuclear gradient.
///
/// @param frame     System geometry (positions in Angstrom, type_id = Z).
/// @param integrals Precomputed integrals at current geometry (Hartree/Bohr).
/// @param result    Converged UHF result (P_alpha, P_beta, C_alpha, C_beta, eigenvalues).
/// @return Forces on all atoms (eV/Angstrom) and total energy (eV).
GradientResult uhf_gradient(
    const Frame& frame,
    const IntegralPack& integrals,
    const UHFResult& result
);

/// @brief Finite-difference gradient for validation.
///
/// Computes numerical gradient via central differences. Cost: 6*n_atoms + 1
/// SCF evaluations. Useful only for testing analytic gradients.
///
/// @param frame   System geometry (positions in Angstrom, type_id = Z).
/// @param charge  Molecular charge (0 = neutral).
/// @param scf_tol SCF convergence tolerance (tight for accurate FD).
/// @param max_iter Maximum SCF iterations.
/// @param h       Displacement step in Angstrom (default 5e-5).
/// @return Forces on all atoms (eV/Angstrom) and total energy (eV).
GradientResult rhf_finite_difference_gradient(
    const Frame& frame,
    int charge = 0,
    Real scf_tol = Real(1e-10),
    int max_iter = 200,
    Real h = Real(5e-5)
);

/// @brief UHF finite-difference gradient for validation.
///
/// Computes numerical gradient via central differences with UHF energy.
/// Cost: 6*n_atoms + 1 SCF evaluations.
///
/// @param frame   System geometry (positions in Angstrom, type_id = Z).
/// @param n_alpha Number of alpha electrons.
/// @param n_beta  Number of beta electrons.
/// @param charge  Molecular charge (0 = neutral).
/// @param scf_tol SCF convergence tolerance (tight for accurate FD).
/// @param max_iter Maximum SCF iterations.
/// @param h       Displacement step [Angstrom] (default 5e-5).
/// @return Forces on all atoms [eV/Angstrom] and total energy [eV].
GradientResult uhf_finite_difference_gradient(
    const Frame& frame,
    int n_alpha,
    int n_beta,
    int charge = 0,
    Real scf_tol = Real(1e-10),
    int max_iter = 200,
    Real h = Real(5e-5)
);

} // namespace scf
} // namespace cpu
} // namespace ATV
