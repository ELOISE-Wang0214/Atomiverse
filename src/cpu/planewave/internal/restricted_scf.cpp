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

#include <algorithm>
#include <cmath>
#include <stdexcept>

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

PWResult run_collinear_scf(const cpu::Frame& frame, const PW& method,
                           PWProblem& prob, int& n_alpha, int& n_beta) {
    if (!method.spin_polarized)
        throw std::invalid_argument(
            "run_collinear_scf: requires method.spin_polarized");
    if (method.smearing != Smearing::none)
        throw std::invalid_argument(
            "run_collinear_scf: metallic smearing not supported");
    if (method.multiplicity < 1)
        throw std::invalid_argument(
            "run_collinear_scf: multiplicity must be >= 1");

    PWResult result{};
    setup_pw_problem(frame, method.upf_path, method.ecut_ry, method.k_mesh,
                     prob);
    PWWorkspace& ws = prob.ws;
    auto& gv0 = prob.gv0;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const int n_elec = prob.n_elec;
    const int mag = method.multiplicity - 1;
    if ((n_elec + mag) % 2 != 0 || mag > n_elec)
        throw std::invalid_argument(
            "run_collinear_scf: n_elec and multiplicity have incompatible "
            "parity");
    n_alpha = (n_elec + mag) / 2;
    n_beta = (n_elec - mag) / 2;
    const int n_occ[2] = {n_alpha, n_beta};
    const int n_bands = std::max(std::max(n_alpha, n_beta), 1);
    const AccReal vol_bohr = prob.vol_bohr;
    const int ng0 = gv0.count();

    ws.mix_beta = AccReal(method.mix_beta);
    ws.pulay_depth = method.pulay_depth;
    ws.nspin = 2;

    const std::size_t ngs = static_cast<std::size_t>(ng0);
    ws.rho_spin_G.resize(2);
    ws.rho_spin_G_old.resize(2);
    ws.V_xc_spin_G.resize(2);
    for (int s = 0; s < 2; ++s) {
        ws.rho_spin_G[s].resize({ngs});
        ws.rho_spin_G_old[s].resize({ngs});
        ws.V_xc_spin_G[s].resize({ngs});
    }
    const AccReal fa =
        (n_elec > 0) ? AccReal(n_alpha) / AccReal(n_elec) : AccReal(0.5);
    const AccReal fb =
        (n_elec > 0) ? AccReal(n_beta) / AccReal(n_elec) : AccReal(0.5);
    for (int ig = 0; ig < ng0; ++ig) {
        ws.rho_spin_G[0][ig] = fa * ws.rho_G[ig];
        ws.rho_spin_G[1][ig] = fb * ws.rho_G[ig];
        ws.rho_spin_G_old[0][ig] = ws.rho_spin_G[0][ig];
        ws.rho_spin_G_old[1][ig] = ws.rho_spin_G[1][ig];
    }

    PoissonSolver poisson;
    XCBuilder xc;
    VeffBuilder veff;
    DensityBuilder density(AccReal(method.mix_beta), method.pulay_depth);
    DavidsonSolver davidson(n_bands, AccReal(method.tol) * AccReal(0.1), 40);

    xc.prepare(ws);
    veff.prepare(ws);
    density.prepare(ws);
    density.reset(ws);
    davidson.prepare(ws);

    for (int iter = 0; iter < method.max_iter; ++iter) {
        for (int ig = 0; ig < ng0; ++ig)
            ws.rho_G[ig] = ws.rho_spin_G[0][ig] + ws.rho_spin_G[1][ig];
        poisson.step(ws);
        xc.step_spin(ws, method.xc);

        xt::xarray<std::complex<AccReal>> rho_out[2];
        AccReal e_band = 0;
        auto* saved_gvec = ws.gvec;
        for (int s = 0; s < 2; ++s) {
            for (int ig = 0; ig < ng0; ++ig)
                ws.V_xc_G[ig] = ws.V_xc_spin_G[s][ig];
            veff.step(ws);
            density.begin(ws);
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
                for (int b = 0; b < n_occ[s]; ++b)
                    e_band += weight * ws.eigenvalues[b];
                density.accumulate_kpoint(ws, weight, n_occ[s], AccReal(1));
                ws.gvec = saved_gvec;
                ws.n_g = ng0;
            }
            density.extract(ws, rho_out[s]);
        }
        for (int s = 0; s < 2; ++s) {
            for (int ig = 0; ig < ng0; ++ig) {
                ws.rho_spin_G_old[s][ig] = ws.rho_spin_G[s][ig];
                ws.rho_spin_G[s][ig] =
                    ws.rho_spin_G_old[s][ig]
                    + AccReal(method.mix_beta)
                          * (rho_out[s][ig] - ws.rho_spin_G_old[s][ig]);
            }
        }
        AccReal dr = 0;
        for (int ig = 0; ig < ng0; ++ig) {
            const auto old = ws.rho_G[ig];
            ws.rho_G[ig] = ws.rho_spin_G[0][ig] + ws.rho_spin_G[1][ig];
            const auto d = ws.rho_G[ig] - old;
            dr += d.real() * d.real() + d.imag() * d.imag();
        }
        result.iterations = iter + 1;
        // Diagnostic total energy from restricted formula (spin energy
        // decomposition is not required for Berry NSCF).
        const auto terms =
            iteration_energy(ws, e_band, prob.ew_total, vol_bohr);
        result.energy_ha = static_cast<double>(terms.total);
        result.energy = result.energy_ha / qe_units::rydberg_to_hartree;
        result.magnetization = static_cast<double>(n_alpha - n_beta);
        if (std::sqrt(static_cast<double>(dr)) < method.tol) {
            result.converged = true;
            break;
        }
    }

    // Freeze total density + spin XC for NSCF; last V_eff is spin β —
    // callers re-bind per channel before diagonalization.
    for (int ig = 0; ig < ng0; ++ig)
        ws.rho_G[ig] = ws.rho_spin_G[0][ig] + ws.rho_spin_G[1][ig];
    poisson.step(ws);
    xc.step_spin(ws, method.xc);
    for (int ig = 0; ig < ng0; ++ig)
        ws.V_xc_G[ig] = ws.V_xc_spin_G[1][ig];
    veff.step(ws);
    return result;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
