/// @file cpu/scf/integral_derivatives.h
/// @brief Analytic derivatives of STO-3G molecular integrals w.r.t. nuclear positions.
///
/// Computes dS/dR, dT/dR, dV/dR, d(ERI)/dR, and dV_nuc/dR for s-type Gaussians.
/// All internal computation in atomic units (Bohr, Hartree).
/// Output in Hartree/Bohr.
///
/// @see cpu::scf::native for the undifferentiated integrals
/// @see cpu::scf::rhf_gradient for the gradient assembly
#pragma once

#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/frame.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Derivatives of all molecular integrals w.r.t. nuclear positions.
///
/// Each derivative tensor is indexed as [atom_index, cart_component, mu, nu]
/// where cart_component = 0,1,2 for x,y,z. All in Hartree/Bohr.
struct IntegralDerivPack {
    xt::xtensor<Real, 4> dS;     ///< dS/dR [n_atoms, 3, n_basis, n_basis] [1/Bohr].
    xt::xtensor<Real, 4> dT;     ///< dT/dR [n_atoms, 3, n_basis, n_basis] [Hartree/Bohr].
    xt::xtensor<Real, 4> dV;     ///< dV/dR [n_atoms, 3, n_basis, n_basis] [Hartree/Bohr].
    xt::xtensor<Real, 6> dERI;   ///< d(ERI)/dR [n_atoms, 3, n_basis, n_basis, n_basis, n_basis] [Hartree/Bohr].
    xt::xtensor<Real, 2> dV_nuc; ///< dV_nuc/dR [n_atoms, 3] [Hartree/Bohr].
    Int n_atoms = 0;             ///< Number of atoms.
    Int n_basis = 0;             ///< Number of basis functions.
};

/// @brief Compute all integral derivatives for STO-3G at current geometry.
///
/// Reads atom positions from frame.x/y/z (Angstrom) and nuclear charges
/// from frame.type_id. Converts to Bohr internally. Output in Hartree/Bohr.
///
/// @param frame CPU Frame with atom positions and types (type_id = Z).
/// @param charge Molecular charge (0 = neutral).
/// @return IntegralDerivPack with all derivative tensors.
IntegralDerivPack build_integral_derivatives_sto3g(const Frame& frame, int charge = 0);

} // namespace scf
} // namespace cpu
} // namespace ATV
