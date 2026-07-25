/// @file cpu/planewave/internal/restricted_scf.cpp
/// @brief Shared restricted fixed-occupation PW SCF loop (bands/DOS/Berry).

#include "atomiverse/cpu/planewave/internal/restricted_scf.h"

#include "atomiverse/cpu/planewave/internal/convergence.h"
#include "atomiverse/cpu/planewave/internal/davidson.h"
#include "atomiverse/cpu/planewave/internal/density_builder.h"
#include "atomiverse/cpu/planewave/internal/poisson.h"
#include "atomiverse/cpu/planewave/internal/veff_builder.h"
#include "atomiverse/cpu/planewave/internal/xc_builder.h"
#include "atomiverse/cpu/planewave/qe_units.h"

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

PWResult run_restricted_scf(const cpu::Frame& frame, const PW& method,
                            PWProblem& prob) {
    PWResult result{};

    setup_pw_problem(frame, method.upf_path, method.ecut_ry, method.k_mesh,
                     prob);
    PWWorkspace& ws = prob.ws;
    auto& gv0 = prob.gv0;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const int n_bands_occ = prob.n_bands;
    const AccReal vol_bohr = prob.vol_bohr;
    ws.mix_beta = AccReal(method.mix_beta);
    ws.pulay_depth = method.pulay_depth;

    PoissonSolver poisson;
    XCBuilder xc;
    VeffBuilder veff;
    DensityBuilder density(AccReal(method.mix_beta), method.pulay_depth);
    DavidsonSolver davidson(n_bands_occ, AccReal(method.tol) * AccReal(0.1),
                            40);

    xc.prepare(ws);
    veff.prepare(ws);
    density.prepare(ws);
    density.reset(ws);
    davidson.prepare(ws);
    ConvergenceCheck convergence{AccReal(method.tol), 2};

    for (int iter = 0; iter < method.max_iter; ++iter) {
        poisson.step(ws);
        xc.step(ws, method.xc);
        veff.step(ws);

        density.begin(ws);
        AccReal e_band = 0;
        auto* saved_gvec = ws.gvec;
        for (int ik = 0; ik < n_kpts; ++ik) {
            auto& gvk = prob.gv_k[static_cast<std::size_t>(ik)];
            ws.gvec = &gvk;
            ws.n_g = gvk.count();
            ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.T_G[ig] = gvk.ekin[ig];
            if (!prob.kb_proj_k.empty())
                ws.kb_proj = prob.kb_proj_k[static_cast<std::size_t>(ik)];
            if (!prob.atom_wfc_k.empty())
                ws.atom_wfc = prob.atom_wfc_k[static_cast<std::size_t>(ik)];

            davidson.step(ws);

            const AccReal weight =
                prob.kg.weights[static_cast<std::size_t>(ik)];
            for (int b = 0; b < n_bands_occ; ++b)
                e_band += weight * ws.eigenvalues[b] * AccReal(2);
            density.accumulate_kpoint(ws, weight, n_bands_occ);

            ws.gvec = saved_gvec;
            ws.n_g = gv0.count();
        }

        density.finish(ws);
        result.iterations = iter + 1;
        const auto terms =
            iteration_energy(ws, e_band, prob.ew_total, vol_bohr);
        result.energy_ha = static_cast<double>(terms.total);
        result.energy = result.energy_ha / qe_units::rydberg_to_hartree;
        log_iteration_energy(terms, ws, prob.ew_total, iter);

        if (convergence.step(ws)) {
            result.converged = true;
            break;
        }
    }

    // Final V_eff consistent with mixed density for NSCF.
    poisson.step(ws);
    xc.step(ws, method.xc);
    veff.step(ws);
    return result;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
