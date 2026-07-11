/// @file cpu/gpw/basis_data/gpw_data.h
/// @brief In-tree GTH pseudopotential / Gaussian basis lookup (gpw-edl-aimd M0).
///
/// GPW's pseudopotential and basis data live in the binary, not in a third-party
/// data directory. This is what lets the GPW suite build and run real molecules
/// with no external package present. See `gpw_data.cpp` for the provenance of
/// every number (Krack 2005; VandeVondele-Hutter 2007; GTH/HGH).
#pragma once

#include <string>

#include "atomiverse/cpu/gpw/internal/cp2k_basis.h"
#include "atomiverse/cpu/gpw/internal/cp2k_io.h"

namespace ATV {
namespace cpu {
namespace gpw {

/// @brief Look up a GTH pseudopotential from the in-tree table.
///
/// @param element Chemical symbol (e.g. "Zn", "O", "H").
/// @param name    Potential name/alias (e.g. "GTH-PBE-q12").
/// @return Parsed GthPotential (radii [Bohr], couplings h [Hartree]).
/// @throws std::runtime_error if the element/name is not in the table (no fallback).
GthPotential builtin_gth_potential(const std::string& element, const std::string& name);

/// @brief Look up a contracted Gaussian basis from the in-tree table.
///
/// @param element Chemical symbol (e.g. "Zn", "O", "H").
/// @param name    Basis name/alias (e.g. "DZVP-MOLOPT-SR-GTH").
/// @return Parsed GaussianBasis (exponents [Bohr^-2]).
/// @throws std::runtime_error if the element/name is not in the table (no fallback).
GaussianBasis builtin_gth_basis(const std::string& element, const std::string& name);

} // namespace gpw
} // namespace cpu
} // namespace ATV
