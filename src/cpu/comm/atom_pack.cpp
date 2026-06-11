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
    const double* x = atom.x.data();
    const double* y = atom.y.data();
    const double* z = atom.z.data();

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
    double* x = atom.x.data();
    double* y = atom.y.data();
    double* z = atom.z.data();

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
    const double* x = atom.x.data();
    const double* y = atom.y.data();
    const double* z = atom.z.data();
    const std::int64_t* tag = atom.tag.data();
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
    double* x = atom.x.data();
    double* y = atom.y.data();
    double* z = atom.z.data();
    std::int64_t* tag = atom.tag.data();
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
    int m = 0;
    // Continuous fields
    buf[m++] = atom.x[i];
    buf[m++] = atom.y[i];
    buf[m++] = atom.z[i];
    buf[m++] = atom.vx[i];
    buf[m++] = atom.vy[i];
    buf[m++] = atom.vz[i];
    buf[m++] = atom.fx[i];
    buf[m++] = atom.fy[i];
    buf[m++] = atom.fz[i];
    buf[m++] = atom.mass[i];
    buf[m++] = atom.charge[i];

    // Integer-cast-as-double fields
    buf[m++] = static_cast<double>(atom.type[i]);
    buf[m++] = static_cast<double>(atom.tag[i]);
    buf[m++] = static_cast<double>(atom.image(i, 0));
    buf[m++] = static_cast<double>(atom.image(i, 1));
    buf[m++] = static_cast<double>(atom.image(i, 2));

    // Per-atom topology — bonds
    const int nb = atom.num_bond[i];
    buf[m++] = static_cast<double>(nb);
    for (int b = 0; b < nb; ++b) {
        buf[m++] = static_cast<double>(atom.bond_atom(i, b));
        buf[m++] = static_cast<double>(atom.bond_type(i, b));
    }

    // Angles
    const int na = atom.num_angle[i];
    buf[m++] = static_cast<double>(na);
    for (int a = 0; a < na; ++a) {
        buf[m++] = static_cast<double>(atom.angle_atom(i, 2 * a + 0));
        buf[m++] = static_cast<double>(atom.angle_atom(i, 2 * a + 1));
        buf[m++] = static_cast<double>(atom.angle_type(i, a));
    }

    // Dihedrals
    const int nd = atom.num_dihedral[i];
    buf[m++] = static_cast<double>(nd);
    for (int d = 0; d < nd; ++d) {
        buf[m++] = static_cast<double>(atom.dihedral_atom(i, 3 * d + 0));
        buf[m++] = static_cast<double>(atom.dihedral_atom(i, 3 * d + 1));
        buf[m++] = static_cast<double>(atom.dihedral_atom(i, 3 * d + 2));
        buf[m++] = static_cast<double>(atom.dihedral_type(i, d));
    }

    // Special-neighbor list. nspecial[i][2] is the total length.
    const int ns0 = atom.nspecial(i, 0);
    const int ns1 = atom.nspecial(i, 1);
    const int ns2 = atom.nspecial(i, 2);
    buf[m++] = static_cast<double>(ns0);
    buf[m++] = static_cast<double>(ns1);
    buf[m++] = static_cast<double>(ns2);
    for (int s = 0; s < ns2; ++s) {
        buf[m++] = static_cast<double>(atom.special(i, s));
    }
    return m;
}

int unpack_exchange(Atom& atom, int slot, const double* buf) {
    int m = 0;
    atom.x[slot] = buf[m++];
    atom.y[slot] = buf[m++];
    atom.z[slot] = buf[m++];
    atom.vx[slot] = buf[m++];
    atom.vy[slot] = buf[m++];
    atom.vz[slot] = buf[m++];
    atom.fx[slot] = buf[m++];
    atom.fy[slot] = buf[m++];
    atom.fz[slot] = buf[m++];
    atom.mass[slot] = buf[m++];
    atom.charge[slot] = buf[m++];

    atom.type[slot] = static_cast<int>(buf[m++]);
    atom.tag[slot] = static_cast<std::int64_t>(buf[m++]);
    atom.image(slot, 0) = static_cast<std::int8_t>(buf[m++]);
    atom.image(slot, 1) = static_cast<std::int8_t>(buf[m++]);
    atom.image(slot, 2) = static_cast<std::int8_t>(buf[m++]);

    // Defensive stride growth: strides converge to identical values on every
    // rank at load time (each rank scans the full global topology), so these
    // are no-ops in practice — they only fire if a frame was populated
    // rank-locally with heterogeneous valence.
    const int nb = static_cast<int>(buf[m++]);
    atom.grow_bond_per_atom(nb);
    atom.num_bond[slot] = nb;
    for (int b = 0; b < nb; ++b) {
        atom.bond_atom(slot, b) = static_cast<std::int64_t>(buf[m++]);
        atom.bond_type(slot, b) = static_cast<int>(buf[m++]);
    }

    const int na = static_cast<int>(buf[m++]);
    atom.grow_angle_per_atom(na);
    atom.num_angle[slot] = na;
    for (int a = 0; a < na; ++a) {
        atom.angle_atom(slot, 2 * a + 0) =
            static_cast<std::int64_t>(buf[m++]);
        atom.angle_atom(slot, 2 * a + 1) =
            static_cast<std::int64_t>(buf[m++]);
        atom.angle_type(slot, a) = static_cast<int>(buf[m++]);
    }

    const int nd = static_cast<int>(buf[m++]);
    atom.grow_dihedral_per_atom(nd);
    atom.num_dihedral[slot] = nd;
    for (int d = 0; d < nd; ++d) {
        atom.dihedral_atom(slot, 3 * d + 0) =
            static_cast<std::int64_t>(buf[m++]);
        atom.dihedral_atom(slot, 3 * d + 1) =
            static_cast<std::int64_t>(buf[m++]);
        atom.dihedral_atom(slot, 3 * d + 2) =
            static_cast<std::int64_t>(buf[m++]);
        atom.dihedral_type(slot, d) = static_cast<int>(buf[m++]);
    }

    const int ns0 = static_cast<int>(buf[m++]);
    const int ns1 = static_cast<int>(buf[m++]);
    const int ns2 = static_cast<int>(buf[m++]);
    atom.grow_maxspecial(ns2);
    atom.nspecial(slot, 0) = ns0;
    atom.nspecial(slot, 1) = ns1;
    atom.nspecial(slot, 2) = ns2;
    for (int s = 0; s < ns2; ++s) {
        atom.special(slot, s) = static_cast<std::int64_t>(buf[m++]);
    }
    return m;
}

int max_exchange_doubles(const Atom& atom) {
    return 11                            // x/y/z, v, f, mass, charge
           + 5                           // type + tag + image[3]
           + 1 + 2 * atom.bond_per_atom  // num_bond + (partner, type) pairs
           + 1 + 3 * atom.angle_per_atom // num_angle + (i, k, type) triples
           + 1 + 4 * atom.dihedral_per_atom // num_dihedral + (i,k,l,type)
           + 3 + atom.maxspecial;        // nspecial[3] + special list
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
