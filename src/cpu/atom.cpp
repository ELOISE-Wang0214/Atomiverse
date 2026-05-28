/// @file src/cpu/atom.cpp
/// @brief cpu::Atom -- impl. Constructor + grow + map_set.

#include "atomiverse/cpu/atom.h"

#include <algorithm>

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

} // namespace cpu
} // namespace ATV
