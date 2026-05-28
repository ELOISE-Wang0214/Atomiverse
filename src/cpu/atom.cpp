/// @file src/cpu/atom.cpp
/// @brief cpu::Atom -- impl. Constructor + grow + map_set.

#include "atomiverse/cpu/atom.h"

#include <algorithm>
#include <cstdio>
#include <stdexcept>

namespace ATV {
namespace cpu {

Atom::Atom(int initial_nmax) {
    if (initial_nmax > 0)
        grow(initial_nmax);
}

void Atom::grow(int target_nmax) {
    if (target_nmax <= nmax)
        return;

    // LAMMPS DELTA growth: double the capacity until we cover the target.
    // Floor at 64 to avoid the degenerate first-call-with-1 path. The
    // resulting `new_nmax` is always >= target_nmax and >= 2 * nmax.
    int new_nmax = (nmax > 0) ? nmax : 64;
    while (new_nmax < target_nmax)
        new_nmax *= 2;

    const std::size_t n = static_cast<std::size_t>(new_nmax);

    // Helper to resize a 1D xtensor preserving the first `nmax` entries
    // and zero-initialising the trailing rows. xt::xtensor::resize does
    // NOT zero-initialise on growth, so we explicitly zero the tail.
    auto resize1d = [&](auto& arr) {
        using value_t = typename std::decay_t<decltype(arr)>::value_type;
        arr.resize({n});
        // Zero-init [nmax, new_nmax). On first grow nmax==0 so this is
        // the entire array.
        for (std::size_t i = static_cast<std::size_t>(nmax); i < n; ++i)
            arr[i] = value_t{0};
    };
    auto resize2d = [&](auto& arr, std::size_t cols) {
        using value_t = typename std::decay_t<decltype(arr)>::value_type;
        arr.resize({n, cols});
        for (std::size_t i = static_cast<std::size_t>(nmax) * cols;
             i < n * cols; ++i)
            arr.data()[i] = value_t{0};
    };

    resize1d(x);
    resize1d(y);
    resize1d(z);
    resize1d(vx);
    resize1d(vy);
    resize1d(vz);
    resize1d(fx);
    resize1d(fy);
    resize1d(fz);
    resize1d(mass);
    resize1d(charge);
    resize1d(type);

    // `tag` is zero-initialised here; downstream materialization writes
    // canonical tag values into [0, nlocal). LAMMPS-style: tag is owned
    // by Atom::data_atoms (the data-file reader), not Atom::grow itself.
    resize1d(tag);

    resize2d(image, 3);

    resize1d(num_bond);
    resize2d(bond_atom, static_cast<std::size_t>(kMaxBondPerAtom));
    resize2d(bond_type, static_cast<std::size_t>(kMaxBondPerAtom));

    resize1d(num_angle);
    resize2d(angle_atom, static_cast<std::size_t>(2 * kMaxAnglePerAtom));
    resize2d(angle_type, static_cast<std::size_t>(kMaxAnglePerAtom));

    resize1d(num_dihedral);
    resize2d(dihedral_atom,
             static_cast<std::size_t>(3 * kMaxDihedralPerAtom));
    resize2d(dihedral_type, static_cast<std::size_t>(kMaxDihedralPerAtom));

    resize2d(nspecial, 3);
    resize2d(special, static_cast<std::size_t>(kMaxSpecial));

    nmax = new_nmax;
}

void Atom::map_set(int begin, int end) {
    map.reserve(map.size() + static_cast<std::size_t>(end - begin));
    for (int i = begin; i < end; ++i) {
        map.emplace(tag[i], i);
    }
}

namespace {

[[noreturn]] void atom_topology_overflow(const char* what, int slot,
                                         int cap, std::int64_t tag) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "cpu::Atom topology overflow: tag %lld atom %s slot %d >= cap %d",
                  static_cast<long long>(tag), what, slot, cap);
    throw std::runtime_error(msg);
}

} // namespace

void Atom::add_bond(std::int64_t ti, std::int64_t tj, int type) {
    // LAMMPS Newton bond ON: smaller-tag atom owns. Bond appears once.
    const std::int64_t owner_tag = std::min(ti, tj);
    const std::int64_t partner_tag = std::max(ti, tj);
    auto it = map.find(owner_tag);
    if (it == map.end())
        return; // not on this rank
    const int i = it->second;
    const int slot = num_bond[i];
    if (slot >= kMaxBondPerAtom)
        atom_topology_overflow("bond", slot, kMaxBondPerAtom, owner_tag);
    bond_atom(i, slot) = partner_tag;
    bond_type(i, slot) = type;
    num_bond[i] = slot + 1;
}

void Atom::add_angle(std::int64_t ti, std::int64_t tj, std::int64_t tk,
                     int type) {
    // Center atom (tj) owns.
    auto it = map.find(tj);
    if (it == map.end())
        return;
    const int i = it->second;
    const int slot = num_angle[i];
    if (slot >= kMaxAnglePerAtom)
        atom_topology_overflow("angle", slot, kMaxAnglePerAtom, tj);
    angle_atom(i, 2 * slot + 0) = ti;
    angle_atom(i, 2 * slot + 1) = tk;
    angle_type(i, slot) = type;
    num_angle[i] = slot + 1;
}

void Atom::add_dihedral(std::int64_t ti, std::int64_t tj, std::int64_t tk,
                        std::int64_t tl, int type) {
    // Atom-j (tj, the 2nd atom) owns.
    auto it = map.find(tj);
    if (it == map.end())
        return;
    const int i = it->second;
    const int slot = num_dihedral[i];
    if (slot >= kMaxDihedralPerAtom)
        atom_topology_overflow("dihedral", slot, kMaxDihedralPerAtom, tj);
    dihedral_atom(i, 3 * slot + 0) = ti;
    dihedral_atom(i, 3 * slot + 1) = tk;
    dihedral_atom(i, 3 * slot + 2) = tl;
    dihedral_type(i, slot) = type;
    num_dihedral[i] = slot + 1;
}

void Atom::build_special_list() {
    // BFS over the bond graph to populate 1-2 / 1-3 / 1-4 partners.
    // Single-rank assumption: every atom's full bond list is locally
    // available via `map`. We scan over OWNER atoms' bond lists (which
    // only contain the smaller-tag side per Newton-bond-ON convention)
    // and reflect bonds to BOTH endpoints when building the local
    // 1-2 neighbor list — special-list is symmetric.
    //
    // First pass: build a per-atom 1-2 list (`special` rows 0..nspecial[0]).
    // Reset.
    for (int i = 0; i < nlocal; ++i) {
        nspecial(i, 0) = 0;
        nspecial(i, 1) = 0;
        nspecial(i, 2) = 0;
    }

    auto append_special = [&](int i, std::int64_t partner_tag,
                              int& count_slot) {
        const int slot = count_slot;
        if (slot >= kMaxSpecial)
            return; // silently truncate (rare in practice)
        // Skip duplicates.
        for (int s = 0; s < slot; ++s)
            if (special(i, s) == partner_tag)
                return;
        special(i, slot) = partner_tag;
        ++count_slot;
    };

    // 1-2 list: scan owner bonds, reflect to both endpoints.
    for (int i = 0; i < nlocal; ++i) {
        const std::int64_t ti = tag[i];
        const int nb = num_bond[i];
        for (int b = 0; b < nb; ++b) {
            const std::int64_t tj = bond_atom(i, b);
            // Append tj to i's special list.
            {
                int c = nspecial(i, 0);
                append_special(i, tj, c);
                nspecial(i, 0) = c;
                nspecial(i, 1) = c;
                nspecial(i, 2) = c;
            }
            // Append ti to j's special list (j must be local for single-rank).
            auto it = map.find(tj);
            if (it != map.end()) {
                int j = it->second;
                int c = nspecial(j, 0);
                append_special(j, ti, c);
                nspecial(j, 0) = c;
                nspecial(j, 1) = c;
                nspecial(j, 2) = c;
            }
        }
    }

    // 1-3 list: for each 1-2 neighbor's 1-2 neighbors, append if not in 1-2.
    for (int i = 0; i < nlocal; ++i) {
        const std::int64_t ti = tag[i];
        const int n12 = nspecial(i, 0);
        int n13 = n12;
        for (int s = 0; s < n12; ++s) {
            const std::int64_t tj = special(i, s);
            auto it = map.find(tj);
            if (it == map.end())
                continue;
            const int j = it->second;
            const int nj12 = nspecial(j, 0);
            for (int t = 0; t < nj12; ++t) {
                const std::int64_t tk = special(j, t);
                if (tk == ti)
                    continue;
                // Check not already in i's 1-2 list (would double-count).
                bool found = false;
                for (int u = 0; u < n13; ++u)
                    if (special(i, u) == tk) {
                        found = true;
                        break;
                    }
                if (found)
                    continue;
                if (n13 >= kMaxSpecial)
                    break;
                special(i, n13++) = tk;
            }
        }
        nspecial(i, 1) = n13;
        nspecial(i, 2) = n13;
    }

    // 1-4 list: for each 1-3 neighbor's 1-2 neighbors, append if not in
    // 1-2 / 1-3 yet.
    for (int i = 0; i < nlocal; ++i) {
        const std::int64_t ti = tag[i];
        const int n12 = nspecial(i, 0);
        const int n13 = nspecial(i, 1);
        int n14 = n13;
        for (int s = n12; s < n13; ++s) {
            const std::int64_t tj = special(i, s);
            auto it = map.find(tj);
            if (it == map.end())
                continue;
            const int j = it->second;
            const int nj12 = nspecial(j, 0);
            for (int t = 0; t < nj12; ++t) {
                const std::int64_t tk = special(j, t);
                if (tk == ti)
                    continue;
                bool found = false;
                for (int u = 0; u < n14; ++u)
                    if (special(i, u) == tk) {
                        found = true;
                        break;
                    }
                if (found)
                    continue;
                if (n14 >= kMaxSpecial)
                    break;
                special(i, n14++) = tk;
            }
        }
        nspecial(i, 2) = n14;
    }
}

} // namespace cpu
} // namespace ATV
