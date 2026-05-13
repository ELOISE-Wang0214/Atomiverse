/// @file cpu/pw/solve_si.cpp
#include "atomiverse/cpu/pw/solve_si.h"
#include "atomiverse/cpu/pw/internal/crystal.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include "atomiverse/cpu/pw/internal/kpoint.h"
#include "atomiverse/cpu/pw/internal/ewald.h"
#include "atomiverse/cpu/pw/internal/upf_reader.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include "atomiverse/cpu/pw/pw_poisson.h"
#include "atomiverse/cpu/pw/pw_lda.h"
#include "atomiverse/cpu/pw/pw_hamiltonian.h"
#include "atomiverse/cpu/pw/pw_veff.h"
#include "atomiverse/cpu/pw/pw_davidson.h"
#include "atomiverse/cpu/pw/pw_density.h"
#include "atomiverse/cpu/pw/pw_convergence.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/core/box.h"
#include "atomiverse/core/units.h"
#include <iostream>
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {

double solve_si(const std::string& upf_path, double ecut_ry,
                int k_mesh, double tol, int max_iter) {
    // 1. Build Si crystal
    auto frame = internal::build_si_diamond(5.43);
    Box box;
    box.cell = frame.box;

    // 2. Load UPF or use toy V_loc
    Real ecut_ht = static_cast<Real>(ecut_ry) * kRy2Hartree;
    int n_elec = 8; // 2 Si atoms × 4 valence
    int n_bands = n_elec / 2; // closed shell

    // 3. G-vector list (Gamma-point for simplicity)
    auto gv = internal::enumerate_g_vectors(box, ecut_ht, Vec3::zeros());

    // 4. k-point mesh
    auto kg = internal::generate_mp_mesh(k_mesh, k_mesh, k_mesh);

    // 5. Ewald energy
    std::vector<Vec3> pos;
    std::vector<Real> charges;
    for (int i = 0; i < frame.n_atoms(); ++i) {
        pos.push_back({frame.x[i], frame.y[i], frame.z[i]});
        charges.push_back(4); // Z_val for Si
    }
    internal::EwaldSummation ewald(pos, charges, box, ecut_ht);

    // 6. Setup workspace
    PWWorkspace ws;
    ws.n_g = gv.count();
    ws.gvec = &gv;
    auto ng = static_cast<std::size_t>(ws.n_g);
    ws.rho_G.resize({ng});
    ws.rho_G_old.resize({ng});
    ws.V_H_G.resize({ng});
    ws.V_xc_G.resize({ng});
    ws.V_eff_G.resize({ng});
    ws.V_loc_G_real.resize({ng});
    ws.rho_r.resize({ng});
    ws.V_xc_r.resize({ng});
    ws.T_G.resize({ng});
    ws.mix_beta = Real(0.4);

    // Initialize density: G=0 component = n_electron / volume
    for (int i = 0; i < ws.n_g; ++i) {
        ws.rho_G[i] = std::complex<Real>(0, 0);
        // Toy local potential: simple Coulomb tail
        Real g2 = gv.g[i].norm_squared();
        if (g2 < Real(1e-16))
            ws.V_loc_G_real[i] = Real(0);
        else
            ws.V_loc_G_real[i] = Real(-4.0) * Real(4.0 * 3.141592653589793) / g2;
    }
    // G=0 density: n_electrons / volume
    ws.rho_G[0] = std::complex<Real>(Real(n_elec) / box.volume(), 0);

    // 7. Create context with workspace
    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame, &ws};

    // 8. Run SCF loop manually
    Real e_total = 0;
    for (int iter = 0; iter < max_iter; ++iter) {
        // Build Hamiltonian phase
        PW_Poisson po; po.build(bctx); po.launch({ctx, frame, 0, iter});
        PW_LDA lda; lda.build(bctx); lda.launch({ctx, frame, 0, iter});
        PW_Hamiltonian h; h.build(bctx); h.launch({ctx, frame, 0, iter});
        PW_Veff v; v.build(bctx); v.launch({ctx, frame, 0, iter});

        // Diagonalize
        PW_Davidson d{n_bands, Real(tol), max_iter};
        d.build(bctx); d.launch({ctx, frame, 1, iter});

        // Build density
        PW_Density den{Real(0.4), 5};
        den.build(bctx); den.launch({ctx, frame, 2, iter});

        // Check convergence
        PW_Convergence conv{Real(tol), 2};
        conv.build(bctx); conv.launch({ctx, frame, 3, iter});

        if (ws.converged) {
            // Compute total energy from eigenvalues
            for (int b = 0; b < n_bands; ++b)
                e_total += static_cast<double>(ws.eigenvalues[b]) * 2.0; // 2 electrons per band
            e_total += static_cast<double>(ewald.total_energy()) / static_cast<double>(kRy2eV);
            break;
        }
    }

    return e_total;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
