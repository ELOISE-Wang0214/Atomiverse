/// @file cpu/scf/integrals.cpp
/// @brief IntegralPack constructor — production compositor over the
///        private integral members (AtomList / ShellList +
///        IntegralPack::overlap/kinetic/nuclear_attraction/eri_dense).
///
/// This is the canonical entry point for assembling an `IntegralPack` from a
/// `Frame` and a `BasisSet`. The single-responsibility integral kernels are
/// private static members of `IntegralPack`, composed here:
///
/// @code
///   IntegralPack pack{frame, cpu::scf::BasisSet::STO3G};
///   // pack.S, pack.T, pack.V, pack.Hcore, pack.V_nuc are filled;
///   // pack.ERI additionally with IntegralPack{frame, basis, 0, true}.
/// @endcode

#include "atomiverse/cpu/scf/internal/integrals.h"

#include "atomiverse/core/log.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/internal/basis_workspace.h"
#include "atomiverse/cpu/scf/internal/one_electron.h"

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
    S = overlap(shell_list);
    T = kinetic(shell_list);
    V = nuclear_attraction(shell_list, atom_list);
    Hcore = T + V;
    ATV_INFO("[integrals] one-electron done ms={:.1f}", elapsed_ms(t_1e));

    if (compute_eri) {
        const auto t_2e = LogClock::now();
        ERI = eri_dense(shell_list);
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
