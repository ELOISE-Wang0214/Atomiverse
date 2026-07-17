/// @file cpu/scf/integrals.cpp
/// @brief IntegralPack constructor — production compositor over the
///        single-responsibility builders (AtomList / ShellList +
///        one_electron.h + two_electron.h).
///
/// This is the canonical entry point for assembling an `IntegralPack` from a
/// `Frame` and a `BasisSet`. Low-level callers that need a single integral
/// can call the underlying builders directly:
///
/// @code
///   AtomList atoms{frame, charge};
///   ShellList shells{atoms, cpu::scf::BasisSet::STO3G};
///   auto S = build_overlap(shells);
///   auto T = build_kinetic(shells);
///   auto V = build_nuclear_attraction(shells, atoms);
///   auto Hcore = xt::eval(T + V);
///   double v_nuc = nuclear_repulsion_energy(atoms);
///   // Optionally: auto ERI = build_eri_dense(shells);
/// @endcode

#include "atomiverse/cpu/scf/internal/integrals.h"

#include "atomiverse/core/log.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/internal/basis_workspace.h"
#include "atomiverse/cpu/scf/internal/one_electron.h"
#include "atomiverse/cpu/scf/internal/two_electron.h"

namespace ATV {
namespace cpu {
namespace scf {

IntegralPack::IntegralPack(const Frame& frame, BasisSet basis, int charge,
                           bool compute_eri) {
    const auto t_total = LogClock::now();

    AtomList atom_list{frame, charge};
    ShellList shell_list{atom_list, basis};

    ATV_INFO(
        "[integrals] start natm={} basis={} n_basis={} electrons={} charge={}",
        static_cast<int>(atom_list.atoms.size()), basis_set_name(basis),
        shell_list.n_basis, atom_list.n_electrons, charge);

    n_basis = shell_list.n_basis;
    n_electrons = atom_list.n_electrons;

    const auto t_1e = LogClock::now();
    S = build_overlap(shell_list);
    T = build_kinetic(shell_list);
    V = build_nuclear_attraction(shell_list, atom_list);
    Hcore = T + V;
    ATV_INFO("[integrals] one-electron done ms={:.1f}", elapsed_ms(t_1e));

    if (compute_eri) {
        const auto t_2e = LogClock::now();
        ERI = build_eri_dense(shell_list);
        ATV_INFO("[integrals] two-electron done ms={:.1f}", elapsed_ms(t_2e));
    }

    V_nuc = nuclear_repulsion_energy(atom_list);

    // IntegralPack preserves owned shell storage because downstream
    // consumers (UKS grid, gradient, ASPC) hold references to it.
    atoms = atom_list.atoms;
    shell_data = shell_list.shell_storage;
    shells.resize(shell_list.n_basis);
    for (int i = 0; i < shell_list.n_basis; ++i) {
        shells[i] = shell_list.shells[i];
        const int np = shells[i].n_prim;
        shells[i].alpha = shell_data[i].data();
        shells[i].coeff = shell_data[i].data() + np;
    }

    ATV_INFO("[integrals] done total_ms={:.1f} v_nuc={:.6f}",
             elapsed_ms(t_total), V_nuc);
}

} // namespace scf
} // namespace cpu
} // namespace ATV
