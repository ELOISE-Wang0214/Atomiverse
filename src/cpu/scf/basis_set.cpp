/// @file cpu/scf/basis_set.cpp
/// @brief Implementation of AtomList and ShellList builders.

#include "atomiverse/cpu/scf/basis_set.h"

#include <stdexcept>
#include <string>
#include <utility>

#include "atomiverse/core/constants.h"
#include "atomiverse/cpu/scf/basis_data/pople_631g.h"
#include "atomiverse/cpu/scf/sto3g_data.h"

namespace ATV {
namespace cpu {
namespace scf {

const char* basis_set_name(BasisSet basis) {
    switch (basis) {
        case BasisSet::STO3G:      return "STO-3G";
        case BasisSet::Pople631G:  return "6-31G";
        case BasisSet::Pople631Gs: return "6-31G*";
        case BasisSet::Def2SVP:    return "def2-SVP";
        case BasisSet::Def2TZVP:   return "def2-TZVP";
    }
    return "unknown";
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

// ── AtomList builder ─────────────────────────────────────────────────

AtomList make_atom_list(const Frame& frame, int charge) {
    const int natm = static_cast<int>(frame.atoms.size());

    AtomList out;
    out.atoms.resize(natm);
    int n_electrons = -charge;
    for (int i = 0; i < natm; ++i) {
        const int Z = static_cast<int>(frame.atoms.type[i]);
        out.atoms[i].Z = Z;
        out.atoms[i].pos[0] = static_cast<double>(frame.atoms.x[i]) * kAngstromToBohr;
        out.atoms[i].pos[1] = static_cast<double>(frame.atoms.y[i]) * kAngstromToBohr;
        out.atoms[i].pos[2] = static_cast<double>(frame.atoms.z[i]) * kAngstromToBohr;
        n_electrons += Z;
    }
    out.n_electrons = n_electrons;
    return out;
}

// ── ShellList builder: STO-3G branch ─────────────────────────────────

namespace {

ShellList build_shell_list_sto3g(const AtomList& atoms) {
    const int natm = static_cast<int>(atoms.atoms.size());

    int n_basis = 0;
    for (int a = 0; a < natm; ++a)
        n_basis += sto3g_n_basis(atoms.atoms[a].Z);

    ShellList out;
    out.n_basis = n_basis;
    out.basis = BasisSet::STO3G;
    out.shells.resize(n_basis);
    out.shell_storage.resize(n_basis);

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        const int Z = atoms.atoms[a].Z;
        if (Z < 1 || Z > kMaxZ) {
            throw std::runtime_error(
                "cpu::scf::build_shell_list: STO-3G has no data for Z=" +
                std::to_string(Z));
        }
        const double cx = atoms.atoms[a].pos[0];
        const double cy = atoms.atoms[a].pos[1];
        const double cz = atoms.atoms[a].pos[2];

        const auto& basis = sto3g_data[Z];
        for (int s = 0; s < basis.n_raw_shells; ++s) {
            const auto& raw = basis.raw[s];

            auto emit_shell = [&](int l, int cart, const double coeffs[3], double norm_factor_at(double)) {
                const int np = 3;  // STO-3G always 3 primitives per contraction.
                auto& store = out.shell_storage[shell_idx];
                store.resize(static_cast<std::size_t>(2 * np));
                for (int p = 0; p < np; ++p) {
                    store[p] = raw.alpha[p];
                    store[np + p] = coeffs[p] * norm_factor_at(raw.alpha[p]);
                }
                auto& sh = out.shells[shell_idx];
                sh.n_prim = np;
                sh.alpha = store.data();
                sh.coeff = store.data() + np;
                sh.center[0] = cx;
                sh.center[1] = cy;
                sh.center[2] = cz;
                sh.l = l;
                sh.cart = cart;
                ++shell_idx;
            };

            if (raw.type == 0) {
                // Pure s shell → 1 AO.
                emit_shell(/*l=*/0, /*cart=*/0, raw.coeff_s, native::phys_norm);
            } else {
                // SP shell → 1 s AO + 3 p AOs (px, py, pz), shared exponents.
                emit_shell(/*l=*/0, /*cart=*/0, raw.coeff_s, native::phys_norm);
                emit_shell(/*l=*/1, /*cart=*/0, raw.coeff_p, native::phys_norm_p);
                emit_shell(/*l=*/1, /*cart=*/1, raw.coeff_p, native::phys_norm_p);
                emit_shell(/*l=*/1, /*cart=*/2, raw.coeff_p, native::phys_norm_p);
            }
        }
    }
    return out;
}

// ── Generic builder driven by a checked-in basis_data table ──────────
//
// Walks a `std::vector<basis_data::AtomBasis>` indexed by Z and emits
// one native::Shell per Cartesian AO, applying the physical GTO
// normalization (phys_norm_s / phys_norm_p) to the stored contraction
// coefficients. Supports l ∈ {0, 1}; d-shells (l = 2) will land in
// Stage B2b alongside the per-Cartesian-component normalization
// plumbing.

ShellList build_shell_list_from_table(
    const AtomList& atoms,
    const std::vector<basis_data::AtomBasis>& table,
    BasisSet basis_id)
{
    const int natm = static_cast<int>(atoms.atoms.size());

    // First pass: count AOs and validate element coverage.
    int n_basis = 0;
    for (int a = 0; a < natm; ++a) {
        const int Z = atoms.atoms[a].Z;
        if (Z < 1 || Z >= static_cast<int>(table.size()) ||
            table[static_cast<std::size_t>(Z)].shells.empty()) {
            throw std::runtime_error(
                std::string("cpu::scf::build_shell_list: '") +
                basis_set_name(basis_id) + "' has no data for Z=" +
                std::to_string(Z));
        }
        for (const auto& sh : table[static_cast<std::size_t>(Z)].shells) {
            if (sh.l == 0)      n_basis += 1;
            else if (sh.l == 1) n_basis += 3;
            else {
                throw std::runtime_error(
                    "cpu::scf::build_shell_list: l=" + std::to_string(sh.l) +
                    " not yet supported (Stage B2b adds d-shells).");
            }
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
                auto& store = out.shell_storage[static_cast<std::size_t>(shell_idx)];
                store.resize(static_cast<std::size_t>(2 * np));
                for (int p = 0; p < np; ++p) {
                    store[static_cast<std::size_t>(p)] = sh.alpha[static_cast<std::size_t>(p)];
                    store[static_cast<std::size_t>(np + p)] =
                        sh.coeff[static_cast<std::size_t>(p)] * norm_at(sh.alpha[static_cast<std::size_t>(p)]);
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

            if (sh.l == 0) {
                emit_shell(0, 0, native::phys_norm);
            } else {
                emit_shell(1, 0, native::phys_norm_p);
                emit_shell(1, 1, native::phys_norm_p);
                emit_shell(1, 2, native::phys_norm_p);
            }
        }
    }
    return out;
}

} // namespace

ShellList build_shell_list(const AtomList& atoms, BasisSet basis) {
    switch (basis) {
        case BasisSet::STO3G:
            return build_shell_list_sto3g(atoms);
        case BasisSet::Pople631G:
            return build_shell_list_from_table(
                atoms, basis_data::pople_631g(), basis);
        case BasisSet::Pople631Gs:
        case BasisSet::Def2SVP:
        case BasisSet::Def2TZVP:
            throw std::runtime_error(
                std::string("cpu::scf::build_shell_list: basis '") +
                basis_set_name(basis) +
                "' is not yet implemented (Stage B2b will add it).");
    }
    throw std::runtime_error("cpu::scf::build_shell_list: unknown BasisSet");
}

} // namespace scf
} // namespace cpu
} // namespace ATV
