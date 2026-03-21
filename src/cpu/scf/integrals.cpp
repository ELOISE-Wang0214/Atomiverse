/// @file cpu/scf/integrals.cpp
/// @brief STO-3G molecular integrals via native C++ implementation.

#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/native_integrals.h"
#include "atomiverse/core/constants.h"

#include <cmath>
#include <vector>

namespace ATV {
namespace cpu {
namespace scf {

// STO-3G basis data: exponents and contraction coefficients per shell.
// Reference: Hehre, Stewart, Pople, J. Chem. Phys. 51, 2657 (1969)
struct STO3GShell {
    int l;             // angular momentum (0 = s)
    int n_prim;        // number of primitives
    double alpha[3];   // exponents
    double coeff[3];   // contraction coefficients (unnormalized)
};

struct STO3GBasis {
    int n_shells;
    STO3GShell shells[2]; // max 2 shells per atom (Li: 1s + 2s)
};

// Element-indexed STO-3G basis data (index 0 unused, 1=H, 2=He, 3=Li)
static const STO3GBasis sto3g_basis[] = {
    {0, {}},  // Z=0: unused
    // Z=1: H — 1 shell (1s)
    {1, {
        {0, 3, {3.42525091, 0.62391373, 0.16885540}, {0.15432897, 0.53532814, 0.44463454}},
    }},
    // Z=2: He — 1 shell (1s)
    {1, {
        {0, 3, {6.36242139, 1.15892300, 0.31364979}, {0.15432897, 0.53532814, 0.44463454}},
    }},
    // Z=3: Li — 2 shells (1s + 2s)
    {2, {
        {0, 3, {16.1195750, 2.93620070, 0.79465370}, { 0.15432897,  0.53532814,  0.44463454}},
        {0, 3, { 0.63628970, 0.14786010, 0.04808870}, {-0.09996723,  0.39951283,  0.70011547}},
    }},
};
static constexpr int kMaxZ = 3;

IntegralPack build_integrals_sto3g(const Frame& frame, int charge) {
    int natm = static_cast<int>(frame.n_atoms());
    Int n_electrons = -charge;  // sum(Z) - charge
    for (Int i = 0; i < natm; ++i)
        n_electrons += frame.type_id[i];

    // Count total basis functions (sum of shells per atom)
    Int n_basis = 0;
    for (Int i = 0; i < natm; ++i) {
        int Z = static_cast<int>(frame.type_id[i]);
        n_basis += sto3g_basis[Z].n_shells;
    }

    IntegralPack pack;
    pack.n_basis = n_basis;
    pack.n_electrons = n_electrons;
    pack.S = xt::zeros<Real>({n_basis, n_basis});
    pack.T = xt::zeros<Real>({n_basis, n_basis});
    pack.V = xt::zeros<Real>({n_basis, n_basis});
    pack.ERI = xt::zeros<Real>({n_basis, n_basis, n_basis, n_basis});

    // Build per-shell normalized coefficients with physical norm (2α/π)^{3/4}
    std::vector<std::array<double, 3>> norm_coeffs(n_basis);
    std::vector<native::Shell> shells(n_basis);
    std::vector<native::Atom> atoms(natm);

    int shell_idx = 0;
    for (int a = 0; a < natm; ++a) {
        int Z = static_cast<int>(frame.type_id[a]);
        double cx = static_cast<double>(frame.x[a]) * kAngstromToBohr;
        double cy = static_cast<double>(frame.y[a]) * kAngstromToBohr;
        double cz = static_cast<double>(frame.z[a]) * kAngstromToBohr;

        atoms[a].Z = Z;
        atoms[a].pos[0] = cx; atoms[a].pos[1] = cy; atoms[a].pos[2] = cz;

        const auto& basis = sto3g_basis[Z];
        for (int s = 0; s < basis.n_shells; ++s) {
            const auto& sh = basis.shells[s];
            for (int p = 0; p < 3; ++p)
                norm_coeffs[shell_idx][p] = sh.coeff[p] * native::phys_norm(sh.alpha[p]);
            shells[shell_idx].n_prim = 3;
            shells[shell_idx].alpha = sh.alpha;
            shells[shell_idx].coeff = norm_coeffs[shell_idx].data();
            shells[shell_idx].center[0] = cx;
            shells[shell_idx].center[1] = cy;
            shells[shell_idx].center[2] = cz;
            ++shell_idx;
        }
    }

    // 1e integrals
    for (int i = 0; i < n_basis; ++i)
        for (int j = 0; j <= i; ++j) {
            pack.S(i, j) = static_cast<Real>(native::overlap_ss(shells[i], shells[j]));
            pack.S(j, i) = pack.S(i, j);
            pack.T(i, j) = static_cast<Real>(native::kinetic_ss(shells[i], shells[j]));
            pack.T(j, i) = pack.T(i, j);
            pack.V(i, j) = static_cast<Real>(native::nuclear_ss(shells[i], shells[j], atoms.data(), natm));
            pack.V(j, i) = pack.V(i, j);
        }

    // 2e integrals
    for (int i = 0; i < n_basis; ++i)
      for (int j = 0; j < n_basis; ++j)
        for (int k = 0; k < n_basis; ++k)
          for (int l = 0; l < n_basis; ++l)
            pack.ERI(i, j, k, l) = static_cast<Real>(
                native::eri_ssss(shells[i], shells[j], shells[k], shells[l]));

    pack.Hcore = pack.T + pack.V;

    // Store atoms and shells for grid-based methods (UKS)
    pack.atoms = atoms;
    pack.shell_data.resize(n_basis);
    pack.shells.resize(n_basis);
    for (int i = 0; i < n_basis; ++i) {
        for (int p = 0; p < 3; ++p) {
            pack.shell_data[i][p] = shells[i].alpha[p];
            pack.shell_data[i][p + 3] = norm_coeffs[i][p];
        }
        pack.shells[i].n_prim = 3;
        pack.shells[i].alpha = &pack.shell_data[i][0];
        pack.shells[i].coeff = &pack.shell_data[i][3];
        pack.shells[i].center[0] = shells[i].center[0];
        pack.shells[i].center[1] = shells[i].center[1];
        pack.shells[i].center[2] = shells[i].center[2];
    }

    // Nuclear repulsion energy
    pack.V_nuc = 0;
    for (Int a = 0; a < natm; ++a)
        for (Int b = a + 1; b < natm; ++b) {
            double dx = (static_cast<double>(frame.x[a]) - static_cast<double>(frame.x[b])) * kAngstromToBohr;
            double dy = (static_cast<double>(frame.y[a]) - static_cast<double>(frame.y[b])) * kAngstromToBohr;
            double dz = (static_cast<double>(frame.z[a]) - static_cast<double>(frame.z[b])) * kAngstromToBohr;
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            pack.V_nuc += static_cast<Real>(frame.type_id[a] * frame.type_id[b] / r);
        }

    return pack;
}

} // namespace scf
} // namespace cpu
} // namespace ATV
