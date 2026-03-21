/// @file cpu/scf/integrals.h
/// @brief IntegralPack + build_integrals_sto3g() — molecular integrals via native C++.
///
/// All arrays are xtensor. Computes overlap (S), kinetic (T), nuclear attraction (V),
/// core Hamiltonian (Hcore = T + V), two-electron repulsion integrals (ERI),
/// and nuclear repulsion energy (V_nuc) for a given molecular geometry.
///
/// Supports elements H (Z=1), He (Z=2), Li (Z=3) with STO-3G basis.
///
/// @see cpu::scf::rhf_solve for the CPU RHF solver that consumes IntegralPack
#pragma once

#include <vector>
#include <array>
#include <xtensor/xtensor.hpp>

#include "atomiverse/core/real.h"
#include "atomiverse/core/int.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Molecular integrals for a Hartree-Fock / Kohn-Sham calculation.
struct IntegralPack {
    xt::xtensor<Real, 2> S;       ///< Overlap matrix [n_basis x n_basis].
    xt::xtensor<Real, 2> T;       ///< Kinetic energy matrix [n_basis x n_basis].
    xt::xtensor<Real, 2> V;       ///< Nuclear attraction matrix [n_basis x n_basis].
    xt::xtensor<Real, 2> Hcore;   ///< Core Hamiltonian T + V [n_basis x n_basis].
    xt::xtensor<Real, 4> ERI;     ///< Two-electron integrals (ij|kl) [n^4].
    Real V_nuc = 0;               ///< Nuclear repulsion energy [Hartree].
    Int n_basis = 0;              ///< Number of basis functions.
    Int n_electrons = 0;          ///< Total number of electrons.

    /// @brief Atom descriptors (Z, position in Bohr). Populated by build_integrals_sto3g.
    std::vector<native::Atom> atoms;

    /// @brief Shell descriptors with owned coefficient storage.
    /// Each shell has n_prim, alpha (pointer into coeff_storage), coeff (pointer into coeff_storage), center.
    std::vector<native::Shell> shells;

    /// @brief Backing storage for shell exponents and normalized coefficients.
    /// Shells point into this storage. Layout: [alpha_0..alpha_{n-1}, coeff_0..coeff_{n-1}] per shell.
    std::vector<std::array<double, 6>> shell_data;

    /// @brief Number of doubly-occupied orbitals (RHF closed-shell).
    /// @return `n_electrons / 2`.
    Int n_occ() const { return n_electrons / 2; }
};

/// @brief Compute STO-3G integrals for a molecular geometry via native C++.
///
/// Reads atom positions from frame.x/y/z (Angstrom) and nuclear charges
/// from frame.type_id. Converts to Bohr internally.
///
/// @param frame CPU Frame with atom positions and types (type_id = Z).
/// @param charge Molecular charge (0 = neutral, +1 = cation, etc.).
/// @return IntegralPack with all integrals computed.
IntegralPack build_integrals_sto3g(const Frame& frame, int charge = 0);

} // namespace scf
} // namespace cpu
} // namespace ATV
