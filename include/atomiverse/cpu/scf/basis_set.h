/// @file cpu/scf/basis_set.h
/// @brief Public basis-set identifier for Gaussian methods.
///
/// Atom lists, shell lists, and integral packs are **not** public types —
/// they are private workspace inside concrete methods (`RHF`, `UHF`, …).
#pragma once

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Identifier for a Gaussian basis set.
enum class BasisSet {
    STO3G,      ///< Minimal STO-3G (s/sp). Production baseline.
    Pople631G,  ///< 6-31G split-valence (no polarization).
    Pople631Gs, ///< 6-31G* (d polarization on heavy atoms).
    Def2SVP,    ///< def2-SVP (Ahlrichs split-valence polarized).
    Def2TZVP,   ///< def2-TZVP (triple-zeta).
};

/// @brief Human-readable name for a BasisSet (logs / reports).
[[nodiscard]] const char* basis_set_name(BasisSet basis);

} // namespace scf
} // namespace cpu
} // namespace ATV
