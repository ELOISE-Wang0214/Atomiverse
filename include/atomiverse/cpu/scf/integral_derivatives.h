/// @file cpu/scf/integral_derivatives.h
/// @brief Analytic derivatives of STO-3G molecular integrals w.r.t. nuclear positions.
///
/// Computes dS/dR, dT/dR, dV/dR, d(ERI)/dR, and dV_nuc/dR for s+p type Gaussians.
/// Supports elements H through O (Z=1..8) with STO-3G basis (s and sp shells).
/// Uses the angular momentum raising formula for derivatives:
///   d/dA_d <a|O|b> = 2*alpha_a * <a+1_d|O|b> - a_d * <a-1_d|O|b>
///
/// All internal computation in atomic units (Bohr, Hartree).
/// Output in Hartree/Bohr.
///
/// @see cpu::scf::native for the undifferentiated integrals
/// @see cpu::scf::scf_gradient for the gradient assembly
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solve.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Derivatives of all molecular integrals w.r.t. nuclear positions.
///
/// Each derivative tensor is indexed as [atom_index, cart_component, mu, nu]
/// where cart_component = 0,1,2 for x,y,z. All in Hartree/Bohr.
struct IntegralDerivPack {
    xt::xtensor<double, 4> dS;     ///< dS/dR [n_atoms, 3, n_basis, n_basis] [1/Bohr].
    xt::xtensor<double, 4> dT;     ///< dT/dR [n_atoms, 3, n_basis, n_basis] [Hartree/Bohr].
    xt::xtensor<double, 4> dV;     ///< dV/dR [n_atoms, 3, n_basis, n_basis] [Hartree/Bohr].
    xt::xtensor<double, 2> dV_nuc; ///< dV_nuc/dR [n_atoms, 3] [Hartree/Bohr].
    int n_atoms = 0;               ///< Number of atoms.
    int n_basis = 0;               ///< Number of basis functions.
};

/// @brief Compute all integral derivatives for STO-3G at current geometry.
///
/// Reads atom positions from frame.atoms.x/y/z (Angstrom) and nuclear charges
/// from frame.atoms.type. Converts to Bohr internally. Output in Hartree/Bohr.
///
/// @param frame CPU Frame with atom positions and types (atoms.type = Z).
/// @param charge Molecular charge (0 = neutral).
/// @return IntegralDerivPack with all derivative tensors.
IntegralDerivPack build_integral_derivatives_sto3g(const Frame& frame, int charge = 0);

/// @brief Contract two-electron derivative integrals directly into the gradient.
///
/// Avoids materializing the full dERI tensor, which scales as O(n_atoms * n_basis^4)
/// in memory. The contraction is performed on the fly using the same analytic
/// primitive derivative formulas as build_integral_derivatives_sto3g().
///
/// @param[out] grad      Gradient accumulator [n_atoms, 3] in Hartree/Bohr.
/// @param integrals      Integral pack for the current geometry.
/// @param P_total        Total density matrix (restricted: P_alpha, unrestricted: Pa+Pb).
/// @param result         Converged SCF result with spin densities.
/// @param config         SCF configuration selecting HF/KS and spin treatment.
void accumulate_two_electron_gradient_sto3g(
    xt::xtensor<double, 2>& grad,
    const Frame& frame,
    const IntegralPack& integrals,
    const xt::xtensor<double, 2>& P_total,
    const SCFResult& result,
    const SCFConfig& config);

} // namespace scf
} // namespace cpu
} // namespace ATV
