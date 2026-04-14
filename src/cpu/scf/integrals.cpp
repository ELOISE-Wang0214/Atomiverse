/// @file cpu/scf/integrals.cpp
/// @brief `build_integrals_sto3g` — thin compositor over the decomposed
///        builders (basis_set.h + one_electron.h + two_electron.h).
///
/// This file exists only for backward compatibility. New code should call
/// the single-responsibility builders directly:
///
/// @code
///   auto atoms  = make_atom_list(frame, charge);
///   auto shells = build_shell_list(atoms, BasisSet::STO3G);
///   auto S = build_overlap(shells);
///   auto T = build_kinetic(shells);
///   auto V = build_nuclear_attraction(shells, atoms);
///   auto Hcore = xt::eval(T + V);
///   double v_nuc = nuclear_repulsion_energy(atoms);
///   // Optionally: auto ERI = build_eri_dense(shells);
/// @endcode

#include "atomiverse/cpu/scf/integrals.h"

#include "atomiverse/core/log.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/one_electron.h"
#include "atomiverse/cpu/scf/two_electron.h"

namespace ATV {
namespace cpu {
namespace scf {

IntegralPack build_integrals(const Frame& frame, BasisSet basis,
                             int charge, bool compute_eri) {
    const auto t_total = LogClock::now();

    auto atom_list  = make_atom_list(frame, charge);
    auto shell_list = build_shell_list(atom_list, basis);

    ATV_INFO(
        "[integrals] start natm={} basis={} n_basis={} electrons={} charge={}",
        static_cast<int>(atom_list.atoms.size()),
        basis_set_name(basis),
        shell_list.n_basis,
        atom_list.n_electrons,
        charge);

    IntegralPack pack;
    pack.n_basis = shell_list.n_basis;
    pack.n_electrons = atom_list.n_electrons;

    const auto t_1e = LogClock::now();
    pack.S = build_overlap(shell_list);
    pack.T = build_kinetic(shell_list);
    pack.V = build_nuclear_attraction(shell_list, atom_list);
    pack.Hcore = pack.T + pack.V;
    ATV_INFO("[integrals] one-electron done ms={:.1f}", elapsed_ms(t_1e));

    if (compute_eri) {
        const auto t_2e = LogClock::now();
        pack.ERI = build_eri_dense(shell_list);
        ATV_INFO("[integrals] two-electron done ms={:.1f}", elapsed_ms(t_2e));
    }

    pack.V_nuc = nuclear_repulsion_energy(atom_list);

    // IntegralPack preserves owned shell storage because downstream
    // consumers (UKS grid, gradient, ASPC) hold references to it.
    pack.atoms = atom_list.atoms;
    pack.shell_data = shell_list.shell_storage;
    pack.shells.resize(shell_list.n_basis);
    for (int i = 0; i < shell_list.n_basis; ++i) {
        pack.shells[i] = shell_list.shells[i];
        const int np = pack.shells[i].n_prim;
        pack.shells[i].alpha = pack.shell_data[i].data();
        pack.shells[i].coeff = pack.shell_data[i].data() + np;
    }

    ATV_INFO(
        "[integrals] done total_ms={:.1f} v_nuc={:.6f}",
        elapsed_ms(t_total), pack.V_nuc);

    return pack;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
