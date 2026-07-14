/// @file src/cpu/comm/atom_pack.cpp
/// @brief Implementation of comm pack/unpack helpers.

#include "atomiverse/cpu/comm/atom_pack.h"

#include <algorithm>
#include <cmath>

#include "atomiverse/cpu/frame/atom.h"
#include "atomiverse/cpu/frame/box.h"

namespace ATV {
namespace cpu {
namespace comm {

int pack_comm(const Atom& atom, const int* sendlist, int n, double sx,
              double sy, double sz, double* buf) {
    const auto* x = atom.x.data();
    const auto* y = atom.y.data();
    const auto* z = atom.z.data();

    int m = 0;
    for (int k = 0; k < n; ++k) {
        const int i = sendlist[k];
        buf[m++] = x[i] + sx;
        buf[m++] = y[i] + sy;
        buf[m++] = z[i] + sz;
    }
    return m;
}

void unpack_comm(Atom& atom, int n, int first, const double* buf) {
    auto* x = atom.x.data();
    auto* y = atom.y.data();
    auto* z = atom.z.data();

    int m = 0;
    for (int k = 0; k < n; ++k) {
        const int i = first + k;
        x[i] = buf[m++];
        y[i] = buf[m++];
        z[i] = buf[m++];
    }
}

int pack_border(const Atom& atom, const int* sendlist, int n, double sx,
                double sy, double sz, double* buf) {
    const auto* x = atom.x.data();
    const auto* y = atom.y.data();
    const auto* z = atom.z.data();
    const auto* tag = atom.tag.data();
    const int* type = atom.type.data();

    int m = 0;
    for (int k = 0; k < n; ++k) {
        const int i = sendlist[k];
        buf[m++] = x[i] + sx;
        buf[m++] = y[i] + sy;
        buf[m++] = z[i] + sz;
        buf[m++] = static_cast<double>(tag[i]);
        buf[m++] = static_cast<double>(type[i]);
        buf[m++] = static_cast<double>(atom.image(i, 0));
        buf[m++] = static_cast<double>(atom.image(i, 1));
        buf[m++] = static_cast<double>(atom.image(i, 2));
    }
    return m;
}

void unpack_border(Atom& atom, int n, int first, const double* buf) {
    auto* x = atom.x.data();
    auto* y = atom.y.data();
    auto* z = atom.z.data();
    auto* tag = atom.tag.data();
    int* type = atom.type.data();

    int m = 0;
    for (int k = 0; k < n; ++k) {
        const int i = first + k;
        x[i] = buf[m++];
        y[i] = buf[m++];
        z[i] = buf[m++];
        tag[i] = static_cast<std::int64_t>(buf[m++]);
        type[i] = static_cast<int>(buf[m++]);
        atom.image(i, 0) = static_cast<std::int8_t>(buf[m++]);
        atom.image(i, 1) = static_cast<std::int8_t>(buf[m++]);
        atom.image(i, 2) = static_cast<std::int8_t>(buf[m++]);
    }
}

int pack_exchange(const Atom& atom, int i, double* buf) {
    // Registry-driven structural encoder (frame-block-01). The wire format is
    // unchanged — schema order IS wire order (core/frame_schema.h): the Core
    // scalars flat, then each counted topology group as `count` followed by
    // `count` items of its payload columns interleaved (partner tags + type),
    // then the nspecial-prefixed special list; exchanged extension columns
    // ride at the very end.
    ColumnRef cols[kMaxExchangeColumns];
    const int n_cols = atom.collect_columns(cols, kMaxExchangeColumns);

    int m = 0;
    // Core section (schema columns only — extras follow the topology).
    for (int k = 0; k < kAtomFieldCount; ++k) {
        const ColumnRef& col = cols[k];
        if (col.group != ColGroup::Core)
            continue;
        for (int c = 0; c < col.cols; ++c)
            buf[m++] = col.get(i, c);
    }
    // Counted groups, in wire order.
    for (const ColGroup g :
         {ColGroup::Bond, ColGroup::Angle, ColGroup::Dihedral,
          ColGroup::Special}) {
        int count = 0;
        for (int k = 0; k < kAtomFieldCount; ++k) {
            const ColumnRef& col = cols[k];
            if (col.group != g || !(col.flags & kColCount))
                continue;
            for (int c = 0; c < col.cols; ++c)
                buf[m++] = col.get(i, c);
            // The last count component is the item count (nspecial's
            // cumulative 1-2/1-3/1-4 triple ends in the total).
            count = static_cast<int>(col.get(i, col.cols - 1));
        }
        for (int item = 0; item < count; ++item)
            for (int k = 0; k < kAtomFieldCount; ++k) {
                const ColumnRef& col = cols[k];
                if (col.group != g || (col.flags & kColCount))
                    continue;
                for (int e = 0; e < col.epi; ++e)
                    buf[m++] = col.get(i, item * col.epi + e);
            }
    }
    // Exchanged extension columns.
    for (int k = kAtomFieldCount; k < n_cols; ++k)
        for (int c = 0; c < cols[k].cols; ++c)
            buf[m++] = cols[k].get(i, c);
    return m;
}

int unpack_exchange(Atom& atom, int slot, const double* buf) {
    // Mirror of pack_exchange, with the defensive stride growth: a migrating
    // atom may arrive with more topology items than the local stride holds,
    // so each group's stride is grown to the incoming count BEFORE its
    // payload is written — and the column views are re-collected afterwards
    // (growth reallocates the wide arrays).
    ColumnRef cols[kMaxExchangeColumns];
    int n_cols = atom.collect_columns(cols, kMaxExchangeColumns);

    int m = 0;
    for (int k = 0; k < kAtomFieldCount; ++k) {
        const ColumnRef& col = cols[k];
        if (col.group != ColGroup::Core)
            continue;
        for (int c = 0; c < col.cols; ++c)
            col.set(slot, c, buf[m++]);
    }
    for (const ColGroup g :
         {ColGroup::Bond, ColGroup::Angle, ColGroup::Dihedral,
          ColGroup::Special}) {
        int count = 0;
        for (int k = 0; k < kAtomFieldCount; ++k) {
            const ColumnRef& col = cols[k];
            if (col.group != g || !(col.flags & kColCount))
                continue;
            for (int c = 0; c < col.cols; ++c)
                col.set(slot, c, buf[m++]);
            count = static_cast<int>(col.get(slot, col.cols - 1));
        }
        if (count > atom.group_capacity(g)) {
            atom.grow_group(g, count);
            n_cols = atom.collect_columns(cols, kMaxExchangeColumns); // growth moved the storage
        }
        for (int item = 0; item < count; ++item)
            for (int k = 0; k < kAtomFieldCount; ++k) {
                const ColumnRef& col = cols[k];
                if (col.group != g || (col.flags & kColCount))
                    continue;
                for (int e = 0; e < col.epi; ++e)
                    col.set(slot, item * col.epi + e, buf[m++]);
            }
    }
    for (int k = kAtomFieldCount; k < n_cols; ++k)
        for (int c = 0; c < cols[k].cols; ++c)
            cols[k].set(slot, c, buf[m++]);
    return m;
}

int max_exchange_doubles(const Atom& atom) {
    ColumnRef cols[kMaxExchangeColumns];
    const int n_cols = atom.collect_columns(cols, kMaxExchangeColumns);
    int m = 0;
    for (int k = 0; k < n_cols; ++k) {
        const ColumnRef& col = cols[k];
        if (col.group == ColGroup::Core || (col.flags & kColCount))
            m += col.cols; // flat values (extras land here too)
        else
            m += col.epi * atom.group_capacity(col.group); // worst-case items
    }
    return m;
}

double compute_bonded_extension_local(const Atom& atom, const Box& box) {
    const auto& map = atom.map;
    if (map.empty())
        return 0.0;

    auto try_pair = [&](double& max_sq, std::int64_t ti, std::int64_t tj) {
        auto it_i = map.find(ti);
        auto it_j = map.find(tj);
        if (it_i == map.end() || it_j == map.end())
            return;
        const int i = it_i->second;
        const int j = it_j->second;
        double dx = atom.x[i] - atom.x[j];
        double dy = atom.y[i] - atom.y[j];
        double dz = atom.z[i] - atom.z[j];
        box.min_image(dx, dy, dz);
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > max_sq)
            max_sq = d2;
    };

    double max_sq = 0.0;

    for (int i = 0; i < atom.nlocal; ++i) {
        const std::int64_t ti = atom.tag[i];

        const int nb = atom.num_bond[i];
        for (int b = 0; b < nb; ++b)
            try_pair(max_sq, ti, atom.bond_atom(i, b));

        const int na = atom.num_angle[i];
        for (int a = 0; a < na; ++a) {
            const std::int64_t tk = atom.angle_atom(i, 2 * a + 0);
            const std::int64_t tm = atom.angle_atom(i, 2 * a + 1);
            try_pair(max_sq, ti, tk);
            try_pair(max_sq, ti, tm);
            try_pair(max_sq, tk, tm);
        }

        const int nd = atom.num_dihedral[i];
        for (int d = 0; d < nd; ++d) {
            const std::int64_t ti2 = atom.dihedral_atom(i, 3 * d + 0);
            const std::int64_t tk = atom.dihedral_atom(i, 3 * d + 1);
            const std::int64_t tl = atom.dihedral_atom(i, 3 * d + 2);
            try_pair(max_sq, ti2, ti); // atom-j is `i`; partner 0 is atom-i of dihedral
            try_pair(max_sq, ti, tk);
            try_pair(max_sq, ti, tl);
            try_pair(max_sq, ti2, tk);
            try_pair(max_sq, ti2, tl);
            try_pair(max_sq, tk, tl);
        }
    }

    return (max_sq > 0.0) ? std::sqrt(max_sq) : 0.0;
}

} // namespace comm
} // namespace cpu
} // namespace ATV
