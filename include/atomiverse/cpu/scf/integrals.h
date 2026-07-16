/// @file cpu/scf/integrals.h
/// @brief IntegralPack + build_integrals() -- molecular integrals via native C++.
///
/// All arrays are xtensor. Computes overlap (S), kinetic (T), nuclear attraction (V),
/// core Hamiltonian (Hcore = T + V), two-electron repulsion integrals (ERI),
/// and nuclear repulsion energy (V_nuc) for a given molecular geometry.
///
/// Supports elements H through O (Z=1..8) with STO-3G basis (s and sp shells).
///
/// @see cpu::scf::SCFSolver for the CPU RHF solver that consumes IntegralPack
#pragma once

#include <array>
#include <vector>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Molecular integrals for a Hartree-Fock / Kohn-Sham calculation.
struct IntegralPack {
    xt::xtensor<double, 2> S; ///< Overlap matrix [n_basis x n_basis].
    xt::xtensor<double, 2> T; ///< Kinetic energy matrix [n_basis x n_basis].
    xt::xtensor<double, 2>
        V; ///< Nuclear attraction matrix [n_basis x n_basis].
    xt::xtensor<double, 2>
        Hcore;                  ///< Core Hamiltonian T + V [n_basis x n_basis].
    xt::xtensor<double, 4> ERI; ///< Two-electron integrals (ij|kl) [n^4].
    double V_nuc = 0;           ///< Nuclear repulsion energy [Hartree].
    int n_basis = 0;            ///< Number of basis functions.
    int n_electrons = 0;        ///< Total number of electrons.

    /// @brief Atom descriptors (Z, position in Bohr). Populated by build_integrals.
    std::vector<native::Atom> atoms;

    /// @brief Shell descriptors with owned coefficient storage.
    /// Each shell has n_prim, alpha (pointer into coeff_storage), coeff (pointer into coeff_storage), center.
    std::vector<native::Shell> shells;

    /// @brief Backing storage for shell exponents and normalized coefficients.
    /// One inner vector per shell, length `2 * n_prim`, laid out as
    /// `[alpha_0 ... alpha_{n_prim-1}, coeff_0 ... coeff_{n_prim-1}]`.
    /// `pack.shells[i].alpha` / `.coeff` point into this storage; any code
    /// that copies `shell_data` must also rebind those pointers.
    std::vector<std::vector<double>> shell_data;

    /// @brief Number of doubly-occupied orbitals (RHF closed-shell).
    /// @return `n_electrons / 2`.
    int n_occ() const { return n_electrons / 2; }
};

/// @brief Compute STO-3G integrals for a molecular geometry via native C++.
///
/// Reads atom positions from the `kAtom_x/y/z` columns (Angstrom) and nuclear
/// charges from `kAtom_type`. Converts to Bohr internally.
///
/// By default the two-electron ERI tensor is NOT computed (compute_eri=false).
/// The CPU SCF solver uses direct SCF via build_shell_pair_list + build_g_rhf_direct,
/// so storing the full O(n^4) ERI tensor is unnecessary for CPU calculations.
/// Pass compute_eri=true when the caller explicitly needs integrals.ERI:
///   - CUDA RHFFock (uploads ERI to device)
///   - Tests that verify ERI values directly
///
/// @param frame       CPU Frame with atom positions and types (atoms.type = Z).
/// @param basis       Which basis set to resolve (STO3G, 6-31G, ...).
/// @param charge      Molecular charge (0 = neutral, +1 = cation, etc.).
/// @param compute_eri If true, compute and store the full ERI tensor [default false].
/// @return IntegralPack with S, T, V, Hcore, V_nuc, shells, atoms (ERI only if requested).
IntegralPack build_integrals(const Frame& frame, BasisSet basis, int charge = 0,
                             bool compute_eri = false);

} // namespace scf
} // namespace cpu
} // namespace ATV
