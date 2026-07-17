/// @file cpu/scf/gradient.h
/// @brief Analytic nuclear gradients for mean-field methods (RHF, UHF, RKS, UKS).
///
/// Primary entry point: `scf_gradient(frame, integrals, SCFResult, SCFConfig)`.
/// Dispatches to HF (Coulomb + exchange) or KS (Coulomb + XC grid) based on config.
///
/// References:
///   Pulay, Mol. Phys. 17, 197 (1969) -- Pulay forces
///   Szabo & Ostlund, Appendix C -- RHF gradient derivation
///   Helgaker, Jorgensen, Olsen, Ch. 9 -- molecular gradients
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/internal/integrals.h"
#include "atomiverse/cpu/scf/internal/scf_types.h"
#include "atomiverse/cpu/scf/scf_result.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Unified analytic nuclear gradient for all mean-field methods.
///
/// Dispatches based on SCFConfig:
///   - HF (Restricted or Unrestricted): Coulomb + exchange two-electron terms
///   - KS (Restricted or Unrestricted): Coulomb only + XC grid contribution
///
/// The UHF two-electron formula naturally reduces to the RHF fused form
/// for restricted cases (P_alpha = P_beta = P/2).
///
/// @param frame     System geometry (positions in Angstrom, type_id = Z).
/// @param integrals Precomputed integrals at current geometry (Hartree/Bohr).
/// @param result    Converged SCF result from SCFSolver::run().
/// @param config    SCF configuration (method, spin, XC type, grid parameters).
/// @return Forces on all atoms (eV/Angstrom) and total energy (eV).
GradientResult scf_gradient(const Frame& frame, const IntegralPack& integrals,
                            const SCFResult& result, const SCFConfig& config);

/// @brief Unified finite-difference gradient for validation.
///
/// Central differences using SCFSolver::run() at displaced geometries.
/// Cost: 6*n_atoms + 1 SCF evaluations.
///
/// @param frame  System geometry (positions in Angstrom, type_id = Z).
/// @param config SCF configuration (method, spin, XC, grid, electron counts).
/// @param charge Molecular charge (0 = neutral).
/// @param h      Displacement step in Angstrom (default 5e-5).
/// @param basis  Which basis to rebuild integrals in (default STO-3G).
/// @return Forces on all atoms (eV/Angstrom) and total energy (eV).
GradientResult scf_finite_difference_gradient(const Frame& frame,
                                              const SCFConfig& config,
                                              int charge = 0, double h = 5e-5,
                                              BasisSet basis = BasisSet::STO3G);

} // namespace scf
} // namespace cpu
} // namespace ATV
