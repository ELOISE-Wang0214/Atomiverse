/// @file src/cpu/frame/atom.cpp
/// @brief cpu::Atom -- impl. Constructor + grow + map_set.

#include "atomiverse/cpu/frame/atom.h"

#include <algorithm>
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

    // LAMMPS DELTA-style chunked growth. The first allocation is exact
    // (callers size the initial frame to their atom count); subsequent
    // grows step by max(nmax/8, 4096) rows. Doubling here (the previous
    // policy) left up to 100% dead headroom on EVERY per-atom array —
    // and the cell list, force scratch, and halo lists all size off this
    // capacity, so the waste multiplied: measured ~2x per-atom RSS vs
    // LAMMPS in multi-rank runs whose ghost influx (~20-35% of owners)
    // pushed the frame just past its initial exact capacity.
    //
    // The reallocation itself is set_capacity's explicit copy-into-fresh-
    // storage. grow() used to trust xt::xtensor::resize to preserve the
    // prefix; it does NOT — the first GCMC insertion past the initial
    // capacity zeroed every live atom, collapsed them onto the origin, and
    // froze the Markov chain on a NaN energy. Pinned by
    // AtomCopyRow/AtomGrow tests in tests/cpu/test_atom_variable_n.cpp.
    int new_nmax;
    if (nmax == 0) {
        new_nmax = std::max(target_nmax, 64);
    } else {
        const int step = std::max(nmax / 8, 4096);
        new_nmax = std::max(target_nmax, nmax + step);
    }
    set_capacity(new_nmax);
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
        tag(i) = static_cast<TagInt>(i);
    realloc2d(image, 3);
    // Count arrays always present; wide bonded arrays only when topology exists.
    realloc1d(num_bond);
    realloc1d(num_angle);
    realloc1d(num_dihedral);
    realloc2d(nspecial, 3);
    if (has_bonded_topology_) {
        realloc2d(bond_atom, static_cast<std::size_t>(bond_per_atom));
        realloc2d(bond_type, static_cast<std::size_t>(bond_per_atom));
        realloc2d(angle_atom, static_cast<std::size_t>(2 * angle_per_atom));
        realloc2d(angle_type, static_cast<std::size_t>(angle_per_atom));
        realloc2d(dihedral_atom,
                  static_cast<std::size_t>(3 * dihedral_per_atom));
        realloc2d(dihedral_type, static_cast<std::size_t>(dihedral_per_atom));
        realloc2d(special, static_cast<std::size_t>(maxspecial));
    }

    nmax = new_nmax;
}

void Atom::ensure_bonded_storage() {
    // The wide arrays themselves are sized by the stride growers — at this
    // point every stride may still be 0 (empty arrays). Just flip the flag so
    // grow()/set_capacity() keep whatever the growers allocate.
    has_bonded_topology_ = true;
}

namespace {

// Re-stride a wide per-atom 2D array from `unit * old_cols` to
// `unit * new_cols` columns, preserving each row's prefix (`unit` is the
// entries-per-slot multiplier: 1 for bond/type arrays, 2 for angle_atom,
// 3 for dihedral_atom). Trailing entries are zero-filled.
template <typename T>
void re_stride(xt::xtensor<T, 2>& arr, std::size_t rows, std::size_t unit,
               std::size_t old_cols, std::size_t new_cols) {
    xt::xtensor<T, 2> fresh = xt::zeros<T>({rows, unit * new_cols});
    const std::size_t old_w = unit * old_cols;
    const std::size_t new_w = unit * new_cols;
    const T* src = arr.data();
    T* dst = fresh.data();
    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t m = 0; m < old_w; ++m)
            dst[i * new_w + m] = src[i * old_w + m];
    arr = std::move(fresh);
}

} // namespace

void Atom::grow_bond_per_atom(int target) {
    if (target <= bond_per_atom)
        return;
    ensure_bonded_storage();
    const std::size_t n = static_cast<std::size_t>(std::max(nmax, 0));
    re_stride(bond_atom, n, 1, static_cast<std::size_t>(bond_per_atom),
              static_cast<std::size_t>(target));
    re_stride(bond_type, n, 1, static_cast<std::size_t>(bond_per_atom),
              static_cast<std::size_t>(target));
    bond_per_atom = target;
}

void Atom::grow_angle_per_atom(int target) {
    if (target <= angle_per_atom)
        return;
    ensure_bonded_storage();
    const std::size_t n = static_cast<std::size_t>(std::max(nmax, 0));
    re_stride(angle_atom, n, 2, static_cast<std::size_t>(angle_per_atom),
              static_cast<std::size_t>(target));
    re_stride(angle_type, n, 1, static_cast<std::size_t>(angle_per_atom),
              static_cast<std::size_t>(target));
    angle_per_atom = target;
}

void Atom::grow_dihedral_per_atom(int target) {
    if (target <= dihedral_per_atom)
        return;
    ensure_bonded_storage();
    const std::size_t n = static_cast<std::size_t>(std::max(nmax, 0));
    re_stride(dihedral_atom, n, 3, static_cast<std::size_t>(dihedral_per_atom),
              static_cast<std::size_t>(target));
    re_stride(dihedral_type, n, 1, static_cast<std::size_t>(dihedral_per_atom),
              static_cast<std::size_t>(target));
    dihedral_per_atom = target;
}

void Atom::grow_maxspecial(int target) {
    if (target <= maxspecial)
        return;
    ensure_bonded_storage();
    const std::size_t n = static_cast<std::size_t>(std::max(nmax, 0));
    re_stride(special, n, 1, static_cast<std::size_t>(maxspecial),
              static_cast<std::size_t>(target));
    maxspecial = target;
}

void Atom::reserve_topology(int bond, int angle, int dihedral, int special_n) {
    ensure_bonded_storage();
    grow_bond_per_atom(bond);
    grow_angle_per_atom(angle);
    grow_dihedral_per_atom(dihedral);
    grow_maxspecial(special_n);
}

void Atom::copy_row(int dst, int src) {
    if (dst == src)
        return;
    x[dst] = x[src];
    y[dst] = y[src];
    z[dst] = z[src];
    vx[dst] = vx[src];
    vy[dst] = vy[src];
    vz[dst] = vz[src];
    fx[dst] = fx[src];
    fy[dst] = fy[src];
    fz[dst] = fz[src];
    mass[dst] = mass[src];
    charge[dst] = charge[src];
    type[dst] = type[src];
    tag[dst] = tag[src];
    image(dst, 0) = image(src, 0);
    image(dst, 1) = image(src, 1);
    image(dst, 2) = image(src, 2);

    // Per-atom topology — must follow the row. Count arrays are always
    // present; the wide partner/type/special arrays exist only for
    // atom_style-molecular systems (lazy bonded storage), so guard their
    // row-move on has_bonded().
    num_bond(dst) = num_bond(src);
    num_angle(dst) = num_angle(src);
    num_dihedral(dst) = num_dihedral(src);
    nspecial(dst, 0) = nspecial(src, 0);
    nspecial(dst, 1) = nspecial(src, 1);
    nspecial(dst, 2) = nspecial(src, 2);
    if (has_bonded()) {
        for (int m = 0; m < bond_per_atom; ++m) {
            bond_atom(dst, m) = bond_atom(src, m);
            bond_type(dst, m) = bond_type(src, m);
        }
        for (int m = 0; m < 2 * angle_per_atom; ++m)
            angle_atom(dst, m) = angle_atom(src, m);
        for (int m = 0; m < angle_per_atom; ++m)
            angle_type(dst, m) = angle_type(src, m);
        for (int m = 0; m < 3 * dihedral_per_atom; ++m)
            dihedral_atom(dst, m) = dihedral_atom(src, m);
        for (int m = 0; m < dihedral_per_atom; ++m)
            dihedral_type(dst, m) = dihedral_type(src, m);
        for (int m = 0; m < maxspecial; ++m)
            special(dst, m) = special(src, m);
    }
}

TagInt Atom::acquire_tag() {
    if (next_tag_ < 0) {
        TagInt top = -1;
        for (int i = 0; i < n_total(); ++i)
            top = tag[i] > top ? tag[i] : top;
        next_tag_ = top + 1;
    }
    return next_tag_++;
}

void Atom::map_set(int begin, int end) {
    map.reserve(map.size() + static_cast<std::size_t>(end - begin));
    for (int i = begin; i < end; ++i) {
        map.emplace(tag[i], i);
    }
}

void Atom::add_bond(TagInt ti, TagInt tj, int type) {
    ensure_bonded_storage(); // lazy atom_style molecular on first topology
    // LAMMPS Newton bond ON: smaller-tag atom owns. Bond appears once.
    const TagInt owner_tag = std::min(ti, tj);
    const TagInt partner_tag = std::max(ti, tj);
    auto it = map.find(owner_tag);
    if (it == map.end())
        return; // not on this rank
    const int i = it->second;
    const int slot = num_bond[i];
    grow_bond_per_atom(slot + 1); // no-op unless this is the widest atom yet
    bond_atom(i, slot) = partner_tag;
    bond_type(i, slot) = type;
    num_bond[i] = slot + 1;
}

void Atom::add_angle(TagInt ti, TagInt tj, TagInt tk,
                     int type) {
    ensure_bonded_storage();
    // Center atom (tj) owns.
    auto it = map.find(tj);
    if (it == map.end())
        return;
    const int i = it->second;
    const int slot = num_angle[i];
    grow_angle_per_atom(slot + 1);
    angle_atom(i, 2 * slot + 0) = ti;
    angle_atom(i, 2 * slot + 1) = tk;
    angle_type(i, slot) = type;
    num_angle[i] = slot + 1;
}

void Atom::add_dihedral(TagInt ti, TagInt tj, TagInt tk,
                        TagInt tl, int type) {
    ensure_bonded_storage();
    // Atom-j (tj, the 2nd atom) owns.
    auto it = map.find(tj);
    if (it == map.end())
        return;
    const int i = it->second;
    const int slot = num_dihedral[i];
    grow_dihedral_per_atom(slot + 1);
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

    auto append_special = [&](int i, TagInt partner_tag,
                              int& count_slot) {
        const int slot = count_slot;
        grow_maxspecial(slot + 1); // runtime stride: never truncates
        // Skip duplicates.
        for (int s = 0; s < slot; ++s)
            if (special(i, s) == partner_tag)
                return;
        special(i, slot) = partner_tag;
        ++count_slot;
    };

    // 1-2 list: scan owner bonds, reflect to both endpoints.
    for (int i = 0; i < nlocal; ++i) {
        const TagInt ti = tag[i];
        const int nb = num_bond[i];
        for (int b = 0; b < nb; ++b) {
            const TagInt tj = bond_atom(i, b);
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
        const TagInt ti = tag[i];
        const int n12 = nspecial(i, 0);
        int n13 = n12;
        for (int s = 0; s < n12; ++s) {
            const TagInt tj = special(i, s);
            auto it = map.find(tj);
            if (it == map.end())
                continue;
            const int j = it->second;
            const int nj12 = nspecial(j, 0);
            for (int t = 0; t < nj12; ++t) {
                const TagInt tk = special(j, t);
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
                grow_maxspecial(n13 + 1);
                special(i, n13++) = tk;
            }
        }
        nspecial(i, 1) = n13;
        nspecial(i, 2) = n13;
    }

    // 1-4 list: for each 1-3 neighbor's 1-2 neighbors, append if not in
    // 1-2 / 1-3 yet.
    for (int i = 0; i < nlocal; ++i) {
        const TagInt ti = tag[i];
        const int n12 = nspecial(i, 0);
        const int n13 = nspecial(i, 1);
        int n14 = n13;
        for (int s = n12; s < n13; ++s) {
            const TagInt tj = special(i, s);
            auto it = map.find(tj);
            if (it == map.end())
                continue;
            const int j = it->second;
            const int nj12 = nspecial(j, 0);
            for (int t = 0; t < nj12; ++t) {
                const TagInt tk = special(j, t);
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
                grow_maxspecial(n14 + 1);
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
    const int MB = bond_per_atom;
    const int MA = angle_per_atom;
    const int MD = dihedral_per_atom;
    const int MS = maxspecial;
    std::vector<int> snap_num_bond(static_cast<std::size_t>(n));
    std::vector<int> snap_num_angle(static_cast<std::size_t>(n));
    std::vector<int> snap_num_dihedral(static_cast<std::size_t>(n));
    std::vector<TagInt> snap_bond_atom(
        static_cast<std::size_t>(n) * MB);
    std::vector<int> snap_bond_type(static_cast<std::size_t>(n)
                                    * MB);
    std::vector<TagInt> snap_angle_atom(
        static_cast<std::size_t>(n) * 2 * MA);
    std::vector<int> snap_angle_type(static_cast<std::size_t>(n)
                                     * MA);
    std::vector<TagInt> snap_dihedral_atom(
        static_cast<std::size_t>(n) * 3 * MD);
    std::vector<int> snap_dihedral_type(static_cast<std::size_t>(n)
                                        * MD);
    std::vector<int> snap_nspecial(static_cast<std::size_t>(n) * 3);
    std::vector<TagInt> snap_special(static_cast<std::size_t>(n)
                                           * MS);

    for (int i = 0; i < n; ++i) {
        snap_num_bond[static_cast<std::size_t>(i)] = num_bond[i];
        snap_num_angle[static_cast<std::size_t>(i)] = num_angle[i];
        snap_num_dihedral[static_cast<std::size_t>(i)] = num_dihedral[i];
        for (int m = 0; m < MB; ++m) {
            snap_bond_atom[static_cast<std::size_t>(i) * MB + m] =
                bond_atom(i, m);
            snap_bond_type[static_cast<std::size_t>(i) * MB + m] =
                bond_type(i, m);
        }
        for (int m = 0; m < 2 * MA; ++m)
            snap_angle_atom[static_cast<std::size_t>(i) * 2 * MA
                            + m] = angle_atom(i, m);
        for (int m = 0; m < MA; ++m)
            snap_angle_type[static_cast<std::size_t>(i) * MA + m] =
                angle_type(i, m);
        for (int m = 0; m < 3 * MD; ++m)
            snap_dihedral_atom[static_cast<std::size_t>(i) * 3
                                   * MD
                               + m] = dihedral_atom(i, m);
        for (int m = 0; m < MD; ++m)
            snap_dihedral_type[static_cast<std::size_t>(i) * MD
                               + m] = dihedral_type(i, m);
        for (int m = 0; m < 3; ++m)
            snap_nspecial[static_cast<std::size_t>(i) * 3 + m] =
                nspecial(i, m);
        for (int m = 0; m < MS; ++m)
            snap_special[static_cast<std::size_t>(i) * MS + m] =
                special(i, m);
    }

    // For each new row i, find the source row whose original atom-tag
    // matches `tag[i]`. Since the snapshot was indexed by original row
    // (and the populator assumed tag == row at populate time), the
    // source row is just `tag[i]`.
    for (int i = 0; i < n; ++i) {
        const TagInt target_tag = tag[i];
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
        for (int m = 0; m < MB; ++m) {
            bond_atom(i, m) = snap_bond_atom[src * MB + m];
            bond_type(i, m) = snap_bond_type[src * MB + m];
        }
        for (int m = 0; m < 2 * MA; ++m)
            angle_atom(i, m) =
                snap_angle_atom[src * 2 * MA + m];
        for (int m = 0; m < MA; ++m)
            angle_type(i, m) = snap_angle_type[src * MA + m];
        for (int m = 0; m < 3 * MD; ++m)
            dihedral_atom(i, m) =
                snap_dihedral_atom[src * 3 * MD + m];
        for (int m = 0; m < MD; ++m)
            dihedral_type(i, m) =
                snap_dihedral_type[src * MD + m];
        for (int m = 0; m < 3; ++m)
            nspecial(i, m) = snap_nspecial[src * 3 + m];
        for (int m = 0; m < MS; ++m)
            special(i, m) = snap_special[src * MS + m];
    }
}

std::size_t Atom::memory_bytes() const {
    auto xb = [](const auto& a) {
        return a.size()
               * sizeof(typename std::decay_t<decltype(a)>::value_type);
    };
    std::size_t n = xb(x) + xb(y) + xb(z) + xb(vx) + xb(vy) + xb(vz) + xb(fx)
                    + xb(fy) + xb(fz) + xb(mass) + xb(charge) + xb(type)
                    + xb(tag) + xb(image) + xb(num_bond) + xb(num_angle)
                    + xb(num_dihedral) + xb(nspecial);
    // Lazily-allocated wide bonded arrays (zero on a pure-pairwise system).
    n += xb(bond_atom) + xb(bond_type) + xb(angle_atom) + xb(angle_type)
         + xb(dihedral_atom) + xb(dihedral_type) + xb(special);
    // tag → local map: bucket array + node payload (key+val+ptr, ~32 B/node).
    n += map.size() * 32 + map.bucket_count() * sizeof(void*);
    return n;
}

} // namespace cpu
} // namespace ATV
