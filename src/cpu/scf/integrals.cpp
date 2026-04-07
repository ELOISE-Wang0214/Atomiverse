/// @file cpu/scf/integrals.cpp
/// @brief STO-3G molecular integrals via native C++ implementation.
///
/// Supports elements Z=1..8 (H through O) with s and sp shells.
/// Uses sto3g_data.h as the single source of truth for basis data.

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/cpu/scf/sto3g_data.h"
#include "atomiverse/core/constants.h"

#include <cmath>
#include <vector>

namespace ATV {
namespace cpu {
namespace scf {

IntegralPack build_integrals_sto3g(const Frame& frame, int charge) {
    const int natm = static_cast<int>(frame.atoms.size());

    int n_electrons = -charge;
    for (int i = 0; i < natm; ++i)
        n_electrons += frame.atoms.type[i];

    // Count total expanded basis functions
    int n_basis = 0;
    for (int i = 0; i < natm; ++i) {
        int Z = static_cast<int>(frame.atoms.type[i]);
        n_basis += sto3g_n_basis(Z);
    }

    IntegralPack pack;
    pack.n_basis = n_basis;
    pack.n_electrons = n_electrons;
    pack.S = xt::zeros<double>({n_basis, n_basis});
    pack.T = xt::zeros<double>({n_basis, n_basis});
    pack.V = xt::zeros<double>({n_basis, n_basis});
    pack.ERI = xt::zeros<double>({n_basis, n_basis, n_basis, n_basis});

    // Temporary storage for shell construction within this function.
    // norm_store owns the exponent+coefficient data that shells point into.
    std::vector<native::Shell> shells(n_basis);
    std::vector<native::Atom> atoms(natm);
    std::vector<std::array<double, 6>> norm_store(n_basis); // alpha[3] + norm_coeff[3]

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        int Z = static_cast<int>(frame.atoms.type[a]);
        double cx = static_cast<double>(frame.atoms.x[a]) * kAngstromToBohr;
        double cy = static_cast<double>(frame.atoms.y[a]) * kAngstromToBohr;
        double cz = static_cast<double>(frame.atoms.z[a]) * kAngstromToBohr;

        atoms[a].Z = Z;
        atoms[a].pos[0] = cx;
        atoms[a].pos[1] = cy;
        atoms[a].pos[2] = cz;

        const auto& basis = sto3g_data[Z];
        for (int s = 0; s < basis.n_raw_shells; ++s) {
            const auto& raw = basis.raw[s];

            if (raw.type == 0) {
                // Pure s shell: 1 basis function
                for (int p = 0; p < 3; ++p) {
                    norm_store[shell_idx][p] = raw.alpha[p];
                    norm_store[shell_idx][p + 3] = raw.coeff_s[p] * native::phys_norm(raw.alpha[p]);
                }
                shells[shell_idx].n_prim = 3;
                shells[shell_idx].alpha = &norm_store[shell_idx][0];
                shells[shell_idx].coeff = &norm_store[shell_idx][3];
                shells[shell_idx].center[0] = cx;
                shells[shell_idx].center[1] = cy;
                shells[shell_idx].center[2] = cz;
                shells[shell_idx].l = 0;
                shells[shell_idx].cart = 0;
                ++shell_idx;
            } else {
                // SP shell: expand to 4 basis functions (2s + 2px + 2py + 2pz)

                // s-component of the sp shell
                for (int p = 0; p < 3; ++p) {
                    norm_store[shell_idx][p] = raw.alpha[p];
                    norm_store[shell_idx][p + 3] = raw.coeff_s[p] * native::phys_norm(raw.alpha[p]);
                }
                shells[shell_idx].n_prim = 3;
                shells[shell_idx].alpha = &norm_store[shell_idx][0];
                shells[shell_idx].coeff = &norm_store[shell_idx][3];
                shells[shell_idx].center[0] = cx;
                shells[shell_idx].center[1] = cy;
                shells[shell_idx].center[2] = cz;
                shells[shell_idx].l = 0;
                shells[shell_idx].cart = 0;
                ++shell_idx;

                // p-components: px (cart=0), py (cart=1), pz (cart=2)
                for (int c = 0; c < 3; ++c) {
                    for (int p = 0; p < 3; ++p) {
                        norm_store[shell_idx][p] = raw.alpha[p];
                        norm_store[shell_idx][p + 3] = raw.coeff_p[p] * native::phys_norm_p(raw.alpha[p]);
                    }
                    shells[shell_idx].n_prim = 3;
                    shells[shell_idx].alpha = &norm_store[shell_idx][0];
                    shells[shell_idx].coeff = &norm_store[shell_idx][3];
                    shells[shell_idx].center[0] = cx;
                    shells[shell_idx].center[1] = cy;
                    shells[shell_idx].center[2] = cz;
                    shells[shell_idx].l = 1;
                    shells[shell_idx].cart = c;
                    ++shell_idx;
                }
            }
        }
    }

    // 1e integrals (general s+p functions)
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j <= i; ++j) {
            pack.S(i, j) = static_cast<double>(native::overlap(shells[i], shells[j]));
            pack.S(j, i) = pack.S(i, j);
            pack.T(i, j) = static_cast<double>(native::kinetic(shells[i], shells[j]));
            pack.T(j, i) = pack.T(i, j);
            pack.V(i, j) = static_cast<double>(native::nuclear(shells[i], shells[j], atoms.data(), natm));
            pack.V(j, i) = pack.V(i, j);
        }

    // 2e integrals
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j < n_basis; ++j)
            for (int k = 0; k < n_basis; ++k)
                for (int l = 0; l < n_basis; ++l)
                    pack.ERI(i, j, k, l) = static_cast<double>(
                        native::eri(shells[i], shells[j], shells[k], shells[l]));

    pack.Hcore = pack.T + pack.V;

    // Store atoms and shells for grid-based methods (UKS).
    // pack.shell_data owns the coefficient storage; pack.shells points into it.
    pack.atoms = atoms;
    pack.shell_data.resize(n_basis);
    pack.shells.resize(n_basis);
    for (int i = 0; i < n_basis; ++i) {
        for (int p = 0; p < 3; ++p) {
            pack.shell_data[i][p] = norm_store[i][p];         // alpha
            pack.shell_data[i][p + 3] = norm_store[i][p + 3]; // normalized coeff
        }
        pack.shells[i].n_prim = 3;
        pack.shells[i].alpha = &pack.shell_data[i][0];
        pack.shells[i].coeff = &pack.shell_data[i][3];
        pack.shells[i].center[0] = shells[i].center[0];
        pack.shells[i].center[1] = shells[i].center[1];
        pack.shells[i].center[2] = shells[i].center[2];
        pack.shells[i].l = shells[i].l;
        pack.shells[i].cart = shells[i].cart;
    }

    // Nuclear repulsion energy
    pack.V_nuc = 0;
    for (int a = 0; a < natm; ++a)
        for (int b = a + 1; b < natm; ++b) {
            double dx = (static_cast<double>(frame.atoms.x[a]) - static_cast<double>(frame.atoms.x[b])) * kAngstromToBohr;
            double dy = (static_cast<double>(frame.atoms.y[a]) - static_cast<double>(frame.atoms.y[b])) * kAngstromToBohr;
            double dz = (static_cast<double>(frame.atoms.z[a]) - static_cast<double>(frame.atoms.z[b])) * kAngstromToBohr;
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            pack.V_nuc += static_cast<double>(frame.atoms.type[a] * frame.atoms.type[b] / r);
        }

    return pack;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
