/// @file cpu/scf/basis_set.h
/// @brief Basis-set identifier + molecule-description types resolved from a
///        Frame (atoms) and a BasisSet (shells).
///
/// These are the layer-1 types used by every SCF-adjacent builder:
///   cpu::scf::make_atom_list     — Frame + charge → AtomList
///   cpu::scf::build_shell_list   — AtomList + BasisSet → ShellList
///
/// AtomList and ShellList are value types. They own their storage and can
/// be moved freely; ShellList fixes up its internal Shell pointers after
/// move/copy so handing a ShellList across function boundaries is safe.
///
/// Together they replace the "one god builder" that previously resolved
/// atoms + shells + integrals in a single call.
#pragma once

#include <array>
#include <vector>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Identifier for a Gaussian basis set.
///
/// Each entry pairs with a data table under
/// `include/atomiverse/cpu/scf/basis_data/`. Stage B1 ships STO-3G only;
/// Stage B2 adds the remaining entries.
enum class BasisSet {
    STO3G,      ///< Minimal STO-3G (s/sp, H–O). Production baseline.
    Pople631G,  ///< 6-31G split-valence (no polarization). Stage B2.
    Pople631Gs, ///< 6-31G* (d polarization on heavy atoms). Stage B2.
    Def2SVP,    ///< def2-SVP (Ahlrichs split-valence polarized). Stage B2.
    Def2TZVP,   ///< def2-TZVP (triple-zeta). Stage B2.
};

/// @brief Human-readable name for a BasisSet, suitable for logs and reports.
[[nodiscard]] const char* basis_set_name(BasisSet basis);

/// @brief Immutable atom list derived from a Frame.
///
/// Stores atom descriptors in Bohr (the unit every downstream integral
/// routine expects) plus the net electron count given the molecular
/// charge. The Frame itself stores positions in Angstrom.
struct AtomList {
    std::vector<native::Atom> atoms; ///< Z + position in Bohr.
    int n_electrons = 0;             ///< Sum(Z) - charge.
};

/// @brief Shell list resolved from an atom list and a basis set.
///
/// `shells` holds flat AO descriptors whose `alpha` / `coeff` point into
/// `shell_storage`. Each inner vector holds `2*n_prim` doubles:
/// `[alpha_0 ... alpha_{n_prim-1}, coeff_0 ... coeff_{n_prim-1}]`.
/// Contraction length is per-shell and can differ across shells — this
/// is essential for any basis beyond STO-3G (e.g. 6-31G has mixed 6/3/1,
/// def2-SVP goes up to 7).
///
/// Copy / move rebinds the pointers automatically so the caller does not
/// need to worry about relocation.
struct ShellList {
    std::vector<native::Shell> shells; ///< Shell descriptors (AO-level).
    std::vector<std::vector<double>>
        shell_storage; ///< [α_0..α_{np-1}, c_0..c_{np-1}] per shell.
    int n_basis = 0;   ///< Total AOs (sum of n_cart over shells).
    BasisSet basis = BasisSet::STO3G; ///< Which basis set produced this list.

    ShellList() = default;
    ShellList(const ShellList& o) { copy_from(o); }
    ShellList(ShellList&& o) noexcept { move_from(std::move(o)); }
    ShellList& operator=(const ShellList& o) {
        if (this != &o)
            copy_from(o);
        return *this;
    }
    ShellList& operator=(ShellList&& o) noexcept {
        if (this != &o)
            move_from(std::move(o));
        return *this;
    }

  private:
    void copy_from(const ShellList& o);
    void move_from(ShellList&& o) noexcept;
    void rebind_pointers();
};

/// @brief Build an AtomList from a Frame + molecular charge.
///
/// Reads `frame.atoms.type` as Z, `frame.atoms.x/y/z` as Angstrom, and
/// converts positions to Bohr. Electron count is `Σ Z − charge`.
[[nodiscard]] AtomList make_atom_list(const Frame& frame, int charge = 0);

/// @brief Resolve a BasisSet over an AtomList into concrete shells + AO count.
///
/// Throws `std::runtime_error` if the basis set has no entry for one of
/// the atom types, or if the basis set is not yet implemented in this
/// Stage (Stage B1 supports STO-3G only).
[[nodiscard]] ShellList build_shell_list(const AtomList& atoms, BasisSet basis);

} // namespace scf
} // namespace cpu
} // namespace ATV
