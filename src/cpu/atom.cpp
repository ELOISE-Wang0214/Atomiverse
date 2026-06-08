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

    // tag default = identity (row index): preserves the legacy
    // `cpu::Frame::Atoms` ctor invariant `tag(i) == i` so callers that
    // skip materialize and write directly into [0, nlocal) still get
    // sensible tags. Downstream `materialize` / `unpack_*` overwrite tag
    // for rows they touch.
    resize1d(tag);
    for (std::size_t i = static_cast<std::size_t>(nmax); i < n; ++i)
        tag[i] = static_cast<std::int64_t>(i);

    resize2d(image, 3);

    // Count arrays stay allocated-and-zeroed always (cheap, ~24 B/atom) so
    // every reader bounded by them is a safe no-op on a bonded-free system.
    // (resize2d zeroes the tail [nmax, new), so these run while nmax is still
    // the OLD capacity — nmax is updated last.)
    resize1d(num_bond);
    resize1d(num_angle);
    resize1d(num_dihedral);
    resize2d(nspecial, 3);

    // The wide partner/type/special arrays (~0.5 KB/atom) are grown only once
    // bonded topology exists (lazy atom_style molecular). A bonded-free system
    // leaves them empty and pays nothing.
    if (has_bonded_topology_) {
        resize2d(bond_atom, static_cast<std::size_t>(kMaxBondPerAtom));
        resize2d(bond_type, static_cast<std::size_t>(kMaxBondPerAtom));
        resize2d(angle_atom, static_cast<std::size_t>(2 * kMaxAnglePerAtom));
        resize2d(angle_type, static_cast<std::size_t>(kMaxAnglePerAtom));
        resize2d(dihedral_atom,
                 static_cast<std::size_t>(3 * kMaxDihedralPerAtom));
        resize2d(dihedral_type, static_cast<std::size_t>(kMaxDihedralPerAtom));
        resize2d(special, static_cast<std::size_t>(kMaxSpecial));
    }

    nmax = new_nmax;
}

void Atom::set_capacity(int new_nmax) {
    if (new_nmax < 0)
        new_nmax = 0;
    if (new_nmax == nmax)
        return;

    const std::size_t n = static_cast<std::size_t>(new_nmax);
    const std::size_t keep =
        std::min(static_cast<std::size_t>(nmax), n); // rows to preserve

    // Reallocate into fresh zeroed storage and copy the preserved prefix
    // explicitly — independent of whether xt::xtensor::resize preserves data,
    // and correct for shrink as well as grow.
    auto realloc1d = [&](auto& arr) {
        using arr_t = std::decay_t<decltype(arr)>;
        using value_t = typename arr_t::value_type;
        arr_t fresh = xt::zeros<value_t>({n});
        for (std::size_t i = 0; i < keep; ++i)
            fresh(i) = arr(i);
        arr = std::move(fresh);
    };
    auto realloc2d = [&](auto& arr, std::size_t cols) {
        using arr_t = std::decay_t<decltype(arr)>;
        using value_t = typename arr_t::value_type;
        arr_t fresh = xt::zeros<value_t>({n, cols});
        const std::size_t row_keep = keep * cols;
        for (std::size_t i = 0; i < row_keep; ++i)
            fresh.data()[i] = arr.data()[i];
        arr = std::move(fresh);
    };

    realloc1d(x);
    realloc1d(y);
    realloc1d(z);
    realloc1d(vx);
    realloc1d(vy);
    realloc1d(vz);
    realloc1d(fx);
    realloc1d(fy);
    realloc1d(fz);
    realloc1d(mass);
    realloc1d(charge);
    realloc1d(type);
    realloc1d(tag);
    // tag default = identity for any newly-added rows (matches grow()).
    for (std::size_t i = keep; i < n; ++i)
        tag(i) = static_cast<std::int64_t>(i);
    realloc2d(image, 3);
    // Count arrays always present; wide bonded arrays only when topology exists.
    realloc1d(num_bond);
    realloc1d(num_angle);
    realloc1d(num_dihedral);
    realloc2d(nspecial, 3);
    if (has_bonded_topology_) {
        realloc2d(bond_atom, static_cast<std::size_t>(kMaxBondPerAtom));
        realloc2d(bond_type, static_cast<std::size_t>(kMaxBondPerAtom));
        realloc2d(angle_atom, static_cast<std::size_t>(2 * kMaxAnglePerAtom));
        realloc2d(angle_type, static_cast<std::size_t>(kMaxAnglePerAtom));
        realloc2d(dihedral_atom,
                  static_cast<std::size_t>(3 * kMaxDihedralPerAtom));
        realloc2d(dihedral_type, static_cast<std::size_t>(kMaxDihedralPerAtom));
        realloc2d(special, static_cast<std::size_t>(kMaxSpecial));
    }

    nmax = new_nmax;
}

void Atom::ensure_bonded_storage() {
    if (has_bonded_topology_)
        return;
    has_bonded_topology_ = true;
    // Allocate the wide arrays to the current capacity, zero-filled. No prior
    // data to preserve (they were never written while bonded was off).
    if (nmax <= 0)
        return; // sized by the next grow()/set_capacity() once atoms exist
    const std::size_t n = static_cast<std::size_t>(nmax);
    bond_atom = xt::zeros<std::int64_t>({n, static_cast<std::size_t>(kMaxBondPerAtom)});
    bond_type = xt::zeros<int>({n, static_cast<std::size_t>(kMaxBondPerAtom)});
    angle_atom = xt::zeros<std::int64_t>(
        {n, static_cast<std::size_t>(2 * kMaxAnglePerAtom)});
    angle_type = xt::zeros<int>({n, static_cast<std::size_t>(kMaxAnglePerAtom)});
    dihedral_atom = xt::zeros<std::int64_t>(
        {n, static_cast<std::size_t>(3 * kMaxDihedralPerAtom)});
    dihedral_type =
        xt::zeros<int>({n, static_cast<std::size_t>(kMaxDihedralPerAtom)});
    special = xt::zeros<std::int64_t>({n, static_cast<std::size_t>(kMaxSpecial)});
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
    ensure_bonded_storage(); // lazy atom_style molecular on first topology
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
    ensure_bonded_storage();
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
    ensure_bonded_storage();
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

void Atom::reshuffle_topology_to_tags() {
    // Bonded-free system: counts are all zero and the wide arrays are empty —
    // nothing to reorder (and the snapshot below would index empty arrays).
    if (!has_bonded_topology_)
        return;
    // Snapshot the current per-atom topology, indexed by row.
    const int n = nlocal;
    std::vector<int> snap_num_bond(static_cast<std::size_t>(n));
    std::vector<int> snap_num_angle(static_cast<std::size_t>(n));
    std::vector<int> snap_num_dihedral(static_cast<std::size_t>(n));
    std::vector<std::int64_t> snap_bond_atom(
        static_cast<std::size_t>(n) * kMaxBondPerAtom);
    std::vector<int> snap_bond_type(static_cast<std::size_t>(n)
                                    * kMaxBondPerAtom);
    std::vector<std::int64_t> snap_angle_atom(
        static_cast<std::size_t>(n) * 2 * kMaxAnglePerAtom);
    std::vector<int> snap_angle_type(static_cast<std::size_t>(n)
                                     * kMaxAnglePerAtom);
    std::vector<std::int64_t> snap_dihedral_atom(
        static_cast<std::size_t>(n) * 3 * kMaxDihedralPerAtom);
    std::vector<int> snap_dihedral_type(static_cast<std::size_t>(n)
                                        * kMaxDihedralPerAtom);
    std::vector<int> snap_nspecial(static_cast<std::size_t>(n) * 3);
    std::vector<std::int64_t> snap_special(static_cast<std::size_t>(n)
                                           * kMaxSpecial);

    for (int i = 0; i < n; ++i) {
        snap_num_bond[static_cast<std::size_t>(i)] = num_bond[i];
        snap_num_angle[static_cast<std::size_t>(i)] = num_angle[i];
        snap_num_dihedral[static_cast<std::size_t>(i)] = num_dihedral[i];
        for (int m = 0; m < kMaxBondPerAtom; ++m) {
            snap_bond_atom[static_cast<std::size_t>(i) * kMaxBondPerAtom + m] =
                bond_atom(i, m);
            snap_bond_type[static_cast<std::size_t>(i) * kMaxBondPerAtom + m] =
                bond_type(i, m);
        }
        for (int m = 0; m < 2 * kMaxAnglePerAtom; ++m)
            snap_angle_atom[static_cast<std::size_t>(i) * 2 * kMaxAnglePerAtom
                            + m] = angle_atom(i, m);
        for (int m = 0; m < kMaxAnglePerAtom; ++m)
            snap_angle_type[static_cast<std::size_t>(i) * kMaxAnglePerAtom + m] =
                angle_type(i, m);
        for (int m = 0; m < 3 * kMaxDihedralPerAtom; ++m)
            snap_dihedral_atom[static_cast<std::size_t>(i) * 3
                                   * kMaxDihedralPerAtom
                               + m] = dihedral_atom(i, m);
        for (int m = 0; m < kMaxDihedralPerAtom; ++m)
            snap_dihedral_type[static_cast<std::size_t>(i) * kMaxDihedralPerAtom
                               + m] = dihedral_type(i, m);
        for (int m = 0; m < 3; ++m)
            snap_nspecial[static_cast<std::size_t>(i) * 3 + m] =
                nspecial(i, m);
        for (int m = 0; m < kMaxSpecial; ++m)
            snap_special[static_cast<std::size_t>(i) * kMaxSpecial + m] =
                special(i, m);
    }

    // For each new row i, find the source row whose original atom-tag
    // matches `tag[i]`. Since the snapshot was indexed by original row
    // (and the populator assumed tag == row at populate time), the
    // source row is just `tag[i]`.
    for (int i = 0; i < n; ++i) {
        const std::int64_t target_tag = tag[i];
        if (target_tag < 0 || target_tag >= n) {
            throw std::runtime_error(
                "cpu::Atom::reshuffle_topology_to_tags: tag out of "
                "range — caller must guarantee tag[0..nlocal) is a "
                "permutation of [0, nlocal)");
        }
        const std::size_t src = static_cast<std::size_t>(target_tag);
        num_bond[i] = snap_num_bond[src];
        num_angle[i] = snap_num_angle[src];
        num_dihedral[i] = snap_num_dihedral[src];
        for (int m = 0; m < kMaxBondPerAtom; ++m) {
            bond_atom(i, m) = snap_bond_atom[src * kMaxBondPerAtom + m];
            bond_type(i, m) = snap_bond_type[src * kMaxBondPerAtom + m];
        }
        for (int m = 0; m < 2 * kMaxAnglePerAtom; ++m)
            angle_atom(i, m) =
                snap_angle_atom[src * 2 * kMaxAnglePerAtom + m];
        for (int m = 0; m < kMaxAnglePerAtom; ++m)
            angle_type(i, m) = snap_angle_type[src * kMaxAnglePerAtom + m];
        for (int m = 0; m < 3 * kMaxDihedralPerAtom; ++m)
            dihedral_atom(i, m) =
                snap_dihedral_atom[src * 3 * kMaxDihedralPerAtom + m];
        for (int m = 0; m < kMaxDihedralPerAtom; ++m)
            dihedral_type(i, m) =
                snap_dihedral_type[src * kMaxDihedralPerAtom + m];
        for (int m = 0; m < 3; ++m)
            nspecial(i, m) = snap_nspecial[src * 3 + m];
        for (int m = 0; m < kMaxSpecial; ++m)
            special(i, m) = snap_special[src * kMaxSpecial + m];
    }
}

} // namespace cpu
} // namespace ATV
