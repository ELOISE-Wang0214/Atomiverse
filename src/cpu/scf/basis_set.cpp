/// @file cpu/scf/basis_set.cpp
/// @brief Implementation of AtomList and ShellList builders.

#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/internal/basis_workspace.h"

#include <stdexcept>
#include <string>
#include <utility>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/basis_data/def2_svp.h"
#include "atomiverse/cpu/scf/basis_data/def2_tzvp.h"
#include "atomiverse/cpu/scf/basis_data/pople_631g.h"
#include "atomiverse/cpu/scf/basis_data/pople_631gs.h"
#include "atomiverse/cpu/scf/basis_data/sto3g.h"
#include "atomiverse/cpu/scf/cart.h"

namespace ATV {
namespace cpu {
namespace scf {

const char* basis_set_name(BasisSet basis) {
    switch (basis) {
    case BasisSet::STO3G:
        return "STO-3G";
    case BasisSet::Pople631G:
        return "6-31G";
    case BasisSet::Pople631Gs:
        return "6-31G*";
    case BasisSet::Def2SVP:
        return "def2-SVP";
    case BasisSet::Def2TZVP:
        return "def2-TZVP";
    }
    return "unknown";
}

basis_data::AngularConvention angular_convention(BasisSet basis) {
    switch (basis) {
    case BasisSet::STO3G:
    case BasisSet::Pople631G:
    case BasisSet::Pople631Gs:
        return basis_data::AngularConvention::Cartesian;
    case BasisSet::Def2SVP:
    case BasisSet::Def2TZVP:
        return basis_data::AngularConvention::Spherical;
    }
    throw std::runtime_error("angular_convention: unknown BasisSet");
}

// ── ShellList copy / move rebind ─────────────────────────────────────

void ShellList::rebind_pointers() {
    const int n = n_basis;
    for (int i = 0; i < n; ++i) {
        const int np = shells[i].n_prim;
        shells[i].alpha = shell_storage[i].data();
        shells[i].coeff = shell_storage[i].data() + np;
    }
}

void ShellList::copy_from(const ShellList& o) {
    shells = o.shells;
    shell_storage = o.shell_storage;
    n_basis = o.n_basis;
    basis = o.basis;
    rebind_pointers();
}

void ShellList::move_from(ShellList&& o) noexcept {
    shells = std::move(o.shells);
    shell_storage = std::move(o.shell_storage);
    n_basis = o.n_basis;
    basis = o.basis;
    rebind_pointers();
}

// ── AtomList constructor ─────────────────────────────────────────────

AtomList::AtomList(const Frame& frame, int charge) {
    const int natm = static_cast<int>(frame.nlocal());

    atoms.resize(natm);
    int nelec = -charge;
    const int* __restrict sp_type = frame.span(kAtom_type);
    const StorageReal* __restrict sp_x = frame.span(kAtom_x);
    const StorageReal* __restrict sp_y = frame.span(kAtom_y);
    const StorageReal* __restrict sp_z = frame.span(kAtom_z);
    for (int i = 0; i < natm; ++i) {
        const int Z = sp_type[i];
        atoms[i].Z = Z;
        // Frame positions are in Angstrom; SCF integrals live in Bohr.
        using units::convert;
        using units::Length;
        atoms[i].pos[0] = convert(sp_x[i], Length::Angstrom, Length::Bohr);
        atoms[i].pos[1] = convert(sp_y[i], Length::Angstrom, Length::Bohr);
        atoms[i].pos[2] = convert(sp_z[i], Length::Angstrom, Length::Bohr);
        nelec += Z;
    }
    n_electrons = nelec;
}

// ── Generic builder driven by a checked-in basis_data table ──────────
namespace {
//
// Walks a `std::vector<basis_data::AtomBasis>` indexed by Z and emits
// one native::Shell per Cartesian AO, applying the physical GTO
// normalization (phys_norm_s / phys_norm_p) to the stored contraction
// coefficients. Supports l ∈ {0, 1}; d-shells (l = 2) will land in
// Stage B2b alongside the per-Cartesian-component normalization
// plumbing.

ShellList
build_shell_list_from_table(const AtomList& atoms,
                            const std::vector<basis_data::AtomBasis>& table,
                            BasisSet basis_id) {
    const int natm = static_cast<int>(atoms.atoms.size());

    // First pass: count AOs and validate element coverage.
    int n_basis = 0;
    for (int a = 0; a < natm; ++a) {
        const int Z = atoms.atoms[a].Z;
        if (Z < 1 || Z >= static_cast<int>(table.size())
            || table[static_cast<std::size_t>(Z)].shells.empty()) {
            throw std::runtime_error(
                std::string("cpu::scf::ShellList: '") + basis_set_name(basis_id)
                + "' has no data for Z=" + std::to_string(Z));
        }
        for (const auto& sh : table[static_cast<std::size_t>(Z)].shells) {
            if (sh.l < 0 || sh.l > kMaxL) {
                throw std::runtime_error("cpu::scf::ShellList: l="
                                         + std::to_string(sh.l)
                                         + " not yet supported.");
            }
            n_basis += n_cart_of(sh.l);
        }
    }

    ShellList out;
    out.n_basis = n_basis;
    out.basis = basis_id;
    out.shells.resize(static_cast<std::size_t>(n_basis));
    out.shell_storage.resize(static_cast<std::size_t>(n_basis));

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        const int Z = atoms.atoms[a].Z;
        const double cx = atoms.atoms[a].pos[0];
        const double cy = atoms.atoms[a].pos[1];
        const double cz = atoms.atoms[a].pos[2];

        for (const auto& sh : table[static_cast<std::size_t>(Z)].shells) {
            const int np = static_cast<int>(sh.alpha.size());

            auto emit_shell = [&](int l, int cart, double norm_at(double)) {
                auto& store =
                    out.shell_storage[static_cast<std::size_t>(shell_idx)];
                store.resize(static_cast<std::size_t>(2 * np));
                for (int p = 0; p < np; ++p) {
                    store[static_cast<std::size_t>(p)] =
                        sh.alpha[static_cast<std::size_t>(p)];
                    store[static_cast<std::size_t>(np + p)] =
                        sh.coeff[static_cast<std::size_t>(p)]
                        * norm_at(sh.alpha[static_cast<std::size_t>(p)]);
                }
                auto& s = out.shells[static_cast<std::size_t>(shell_idx)];
                s.n_prim = np;
                s.alpha = store.data();
                s.coeff = store.data() + np;
                s.center[0] = cx;
                s.center[1] = cy;
                s.center[2] = cz;
                s.l = l;
                s.cart = cart;
                ++shell_idx;
            };

            // Emit one native::Shell per Cartesian component. All
            // components of a given l share the same orientation-independent
            // normalization (libcint / PySCF cart convention).
            const int n_cart = n_cart_of(sh.l);
            double (*norm_at)(double) = (sh.l == 0)   ? native::phys_norm
                                        : (sh.l == 1) ? native::phys_norm_p
                                        : (sh.l == 2) ? native::phys_norm_d
                                                      : native::phys_norm_f;
            for (int cart = 0; cart < n_cart; ++cart)
                emit_shell(sh.l, cart, norm_at);
        }
    }
    return out;
}

} // namespace

ShellList::ShellList(const AtomList& atoms, BasisSet basis_id) {
    ShellList tmp;
    switch (basis_id) {
    case BasisSet::STO3G:
        tmp = build_shell_list_from_table(atoms, basis_data::sto3g(), basis_id);
        break;
    case BasisSet::Pople631G:
        tmp = build_shell_list_from_table(atoms, basis_data::pople_631g(),
                                          basis_id);
        break;
    case BasisSet::Pople631Gs:
        tmp = build_shell_list_from_table(atoms, basis_data::pople_631gs(),
                                          basis_id);
        break;
    case BasisSet::Def2SVP:
        tmp = build_shell_list_from_table(atoms, basis_data::def2_svp(),
                                          basis_id);
        break;
    case BasisSet::Def2TZVP:
        tmp = build_shell_list_from_table(atoms, basis_data::def2_tzvp(),
                                          basis_id);
        break;
    default:
        throw std::runtime_error("cpu::scf::ShellList: unknown BasisSet");
    }
    *this = std::move(tmp);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
