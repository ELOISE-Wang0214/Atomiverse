/// @file cpu/scf/one_electron.h
/// @brief One-electron integral builders for Hartree-Fock / DFT.
///
/// Each free function has a single responsibility. Composition is the
/// caller's job:
///
/// @code
///   auto atoms  = AtomList{frame, charge};
///   auto shells = ShellList{atoms, BasisSet::STO3G};
///   auto S = build_overlap(shells);
///   auto T = build_kinetic(shells);
///   auto V = build_nuclear_attraction(shells, atoms);
///   auto Hcore = T + V;  // xtensor element-wise expression
///   double vnn = nuclear_repulsion_energy(atoms);
/// @endcode
///
/// These builders use the native integral primitives in native_integrals.h.
/// The libcint reference path stays under test utilities and is not used
/// in production code.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/internal/basis_workspace.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Overlap matrix S_ij = <φ_i|φ_j>.
/// @param shells Resolved shell list.
/// @return (n_basis × n_basis) symmetric matrix.
[[nodiscard]] xt::xtensor<double, 2> build_overlap(const ShellList& shells);

/// @brief Kinetic-energy matrix T_ij = <φ_i|-½∇²|φ_j>.
/// @param shells Resolved shell list.
/// @return (n_basis × n_basis) symmetric matrix.
[[nodiscard]] xt::xtensor<double, 2> build_kinetic(const ShellList& shells);

/// @brief Nuclear-attraction matrix V_ij = -Σ_A Z_A <φ_i|1/|r - R_A||φ_j>.
/// @param shells Resolved shell list.
/// @param atoms  Atoms supplying nuclear charges and positions.
/// @return (n_basis × n_basis) symmetric matrix.
[[nodiscard]] xt::xtensor<double, 2>
build_nuclear_attraction(const ShellList& shells, const AtomList& atoms);

/// @brief Classical nucleus-nucleus repulsion energy
///        V_nuc = Σ_{A<B} Z_A Z_B / |R_A - R_B|, in Hartree.
/// @param atoms AtomList with positions in Bohr.
/// @return Nuclear repulsion energy [Hartree].
[[nodiscard]] double nuclear_repulsion_energy(const AtomList& atoms);

} // namespace scf
} // namespace cpu
} // namespace ATV
