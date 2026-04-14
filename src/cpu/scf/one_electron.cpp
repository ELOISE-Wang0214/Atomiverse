/// @file cpu/scf/one_electron.cpp
/// @brief Native-engine implementations of the one-electron integral builders.

#include "atomiverse/cpu/scf/one_electron.h"

#include <cmath>

#include "atomiverse/cpu/scf/native_integrals.h"

namespace ATV {
namespace cpu {
namespace scf {

xt::xtensor<double, 2> build_overlap(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> S = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double v = native::overlap(shells.shells[i], shells.shells[j]);
            S(i, j) = v;
            S(j, i) = v;
        }
    }
    return S;
}

xt::xtensor<double, 2> build_kinetic(const ShellList& shells) {
    const int n = shells.n_basis;
    xt::xtensor<double, 2> T = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double v = native::kinetic(shells.shells[i], shells.shells[j]);
            T(i, j) = v;
            T(j, i) = v;
        }
    }
    return T;
}

xt::xtensor<double, 2> build_nuclear_attraction(
    const ShellList& shells, const AtomList& atoms)
{
    const int n = shells.n_basis;
    const int natm = static_cast<int>(atoms.atoms.size());
    xt::xtensor<double, 2> V = xt::zeros<double>({n, n});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double v = native::nuclear(
                shells.shells[i], shells.shells[j],
                atoms.atoms.data(), natm);
            V(i, j) = v;
            V(j, i) = v;
        }
    }
    return V;
}

double nuclear_repulsion_energy(const AtomList& atoms) {
    const int natm = static_cast<int>(atoms.atoms.size());
    double e = 0.0;
    for (int a = 0; a < natm; ++a) {
        for (int b = a + 1; b < natm; ++b) {
            const double dx = atoms.atoms[a].pos[0] - atoms.atoms[b].pos[0];
            const double dy = atoms.atoms[a].pos[1] - atoms.atoms[b].pos[1];
            const double dz = atoms.atoms[a].pos[2] - atoms.atoms[b].pos[2];
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            e += static_cast<double>(atoms.atoms[a].Z * atoms.atoms[b].Z) / r;
        }
    }
    return e;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
