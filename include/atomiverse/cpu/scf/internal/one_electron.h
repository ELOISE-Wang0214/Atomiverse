/// @file cpu/scf/one_electron.h
/// @brief Nuclear repulsion energy for Hartree-Fock / DFT.
///
/// The one-electron integral kernels (overlap / kinetic / nuclear attraction)
/// are private static members of `IntegralPack`
/// (`IntegralPack::overlap/kinetic/nuclear_attraction`, defined in
/// one_electron.cpp) — construct an `IntegralPack{frame, basis}` to obtain
/// S / T / V / Hcore:
///
/// @code
///   auto atoms  = AtomList{frame, charge};
///   double vnn = nuclear_repulsion_energy(atoms);
/// @endcode
///
/// These kernels use the native integral primitives in native_integrals.h.
/// The libcint reference path stays under test utilities and is not used
/// in production code.
#pragma once

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/internal/basis_workspace.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Classical nucleus-nucleus repulsion energy
///        V_nuc = Σ_{A<B} Z_A Z_B / |R_A - R_B|, in Hartree.
/// @param atoms AtomList with positions in Bohr.
/// @return Nuclear repulsion energy [Hartree].
[[nodiscard]] double nuclear_repulsion_energy(const AtomList& atoms);

} // namespace scf
} // namespace cpu
} // namespace ATV
