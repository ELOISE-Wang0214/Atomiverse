/// @file cpu/scf/internal/basis_workspace.h
/// @brief **PRIVATE** AtomList / ShellList workspace — not a user type.
///
/// Frame is the public molecule. These types are Bohr-unit scratch for
/// integral kernels only.
#pragma once

#include <vector>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Private: atoms in Bohr + electron count (from Frame + charge).
struct AtomList {
    std::vector<native::Atom> atoms;
    int n_electrons = 0;

    AtomList() = default;
    AtomList(const Frame& frame, int charge = 0);
};

/// @brief Private: shells + owned coefficient storage for one basis.
struct ShellList {
    std::vector<native::Shell> shells;
    std::vector<std::vector<double>> shell_storage;
    int n_basis = 0;
    BasisSet basis = BasisSet::STO3G;

    ShellList() = default;
    ShellList(const AtomList& atoms, BasisSet basis);

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

} // namespace scf
} // namespace cpu
} // namespace ATV
