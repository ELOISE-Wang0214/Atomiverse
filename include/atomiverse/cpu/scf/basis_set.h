/// @file cpu/scf/basis_set.h
/// @brief Public basis-set identifier for Gaussian methods.
///
/// Atom lists, shell lists, and integral packs are **not** public types —
/// they are private workspace inside concrete methods (`RHF`, `UHF`, …).
#pragma once

#include "atomiverse/cpu/scf/basis_data/basis_table.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Identifier for a Gaussian basis set.
enum class BasisSet {
    STO3G,       ///< Minimal STO-3G (s/sp). Production baseline.
    Pople631G,   ///< 6-31G split-valence (no polarization).
    Pople631Gs,  ///< 6-31G* (d polarization on heavy atoms).
    Pople631pGs, ///< 6-31+G* (heavy-atom diffuse sp + d polarization).
    Def2SVP,     ///< def2-SVP (Ahlrichs split-valence polarized).
    Def2TZVP,    ///< def2-TZVP (triple-zeta).
    AugCcPVDZ,   ///< aug-cc-pVDZ (Dunning DZ + diffuse augmentation).
};

/// @brief Human-readable name for a BasisSet (logs / reports).
[[nodiscard]] const char* basis_set_name(BasisSet basis);

/// @brief Angular convention the basis set is defined with.
///
/// STO-3G / 6-31G / 6-31G* / 6-31+G* are Cartesian (6d by definition);
/// def2-SVP / def2-TZVP / aug-cc-pVDZ are spherical (5d/7f). Throws on an
/// unknown enum value.
/// @param basis Basis-set identifier.
/// @return The defining convention (no fallback).
[[nodiscard]] basis_data::AngularConvention angular_convention(BasisSet basis);

} // namespace scf
} // namespace cpu
} // namespace ATV
