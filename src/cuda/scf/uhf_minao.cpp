/// @file cuda/scf/uhf_minao.cpp
/// @brief Host overlap projection of the PySCF MINAO guess.
///
/// pyscf.scf.hf.init_guess_by_minao, then dma = dmb = 1/2 P
/// (pyscf.scf.uhf.UHF.init_guess_by_minao; _break_dm_spin_symm is a no-op
/// when the spin is nonzero). Not the SCF density.

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/scf/basis_set.h"
#include "atomiverse/cpu/scf/internal/basis_workspace.h"
#include "atomiverse/cpu/scf/internal/integral_view.h"
#include "atomiverse/cpu/scf/internal/minao.h"
#include "atomiverse/cpu/scf/minao_data.h"
#include "atomiverse/cpu/scf/native_integrals.h"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor/containers/xtensor.hpp>

namespace ATV {
namespace cuda {
namespace scf {

/// @return Row-major 1/2 P_MINAO, or empty when an element has no MINAO table.
std::vector<double> uhf_minao_half(const int* Z, const double* x_A,
                                   const double* y_A, const double* z_A,
                                   int natm, cpu::scf::BasisSet basis,
                                   int charge, int n_ao) {
    if (natm <= 0 || Z == nullptr)
        return {};
    for (int i = 0; i < natm; ++i) {
        if (!cpu::scf::minao::has_minao(Z[i]))
            return {};
    }
    cpu::scf::AtomList atoms;
    atoms.atoms.resize(static_cast<std::size_t>(natm));
    atoms.basis_z.resize(static_cast<std::size_t>(natm));
    int zsum = 0;
    for (int i = 0; i < natm; ++i) {
        zsum += Z[i];
        atoms.basis_z[static_cast<std::size_t>(i)] = Z[i];
        atoms.atoms[static_cast<std::size_t>(i)].Z = Z[i];
        auto* pos = atoms.atoms[static_cast<std::size_t>(i)].pos;
        const double xyz[3] = {x_A[i], y_A[i], z_A[i]};
        for (int d = 0; d < 3; ++d)
            pos[d] = units::convert(xyz[d], units::Length::Angstrom,
                                    units::Length::Bohr);
    }
    cpu::scf::ShellList shells(atoms, basis);
    const int n = shells.n_basis;
    if (n != n_ao) {
        throw std::runtime_error("cuda::scf::UHF: MINAO AO count "
                                 + std::to_string(n)
                                 + " does not match the device basis ("
                                 + std::to_string(n_ao) + ").");
    }
    xt::xtensor<double, 2> S = xt::zeros<double>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            const double sij = cpu::scf::native::overlap(
                shells.shells[static_cast<std::size_t>(i)],
                shells.shells[static_cast<std::size_t>(j)]);
            S(i, j) = sij;
            S(j, i) = sij;
        }
    }
    struct Hold {
        xt::xtensor<double, 2> T, V, Hcore, W, c2s;
        xt::xtensor<double, 4> ERI;
        double V_nuc = 0;
        int n_basis = 0;
        int n_basis_cart = 0;
        int n_electrons = 0;
        std::string ecp;
        bool has_ecp = false;
        std::vector<cpu::scf::internal::EcpCenter> ecp_centers;
        std::vector<cpu::EmbeddingSite> embed;
        bool embed_self = false;
        std::array<double, 3> field{};
        double solvent_eps = 0, solvent_radii = 0, solvent_x = 0;
        cpu::scf::SolventModel solvent = cpu::scf::SolventModel::CPCM;
    } hold;
    hold.n_basis = n;
    hold.n_basis_cart = n;
    hold.n_electrons = zsum - charge;
    cpu::scf::internal::IntegralView view{S,
                                          hold.T,
                                          hold.V,
                                          hold.Hcore,
                                          hold.W,
                                          hold.ERI,
                                          hold.V_nuc,
                                          hold.n_basis,
                                          hold.n_basis_cart,
                                          hold.c2s,
                                          hold.n_electrons,
                                          hold.ecp,
                                          hold.has_ecp,
                                          hold.ecp_centers,
                                          hold.embed,
                                          hold.embed_self,
                                          hold.field,
                                          hold.solvent_eps,
                                          hold.solvent_radii,
                                          hold.solvent_x,
                                          hold.solvent,
                                          atoms.atoms,
                                          shells.shells,
                                          shells.shell_storage};
    const cpu::scf::RealMatrix Pmin = cpu::scf::build_minao_density(view);
    std::vector<double> half(static_cast<std::size_t>(n * n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            half[static_cast<std::size_t>(i * n + j)] = 0.5 * Pmin(i, j);
    return half;
}

} // namespace scf
} // namespace cuda
} // namespace ATV
