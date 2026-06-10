/// @file cpu/pw/pw_scf.cpp
/// @brief pw::pw_scf() — generic plane-wave DFT total energy.
///
/// Problem assembly (Ewald, G-vectors, V_loc, projectors, initial density)
/// lives in internal/pw_setup.cpp and is shared with the CUDA backend
/// (`cuda::pw::pw_scf`); this file owns only the CPU SCF loop.

#include "atomiverse/cpu/pw/pw_scf.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/pw/internal/g_vector.h"
#include "atomiverse/cpu/pw/internal/pw_setup.h"
#include "atomiverse/cpu/pw/pw_convergence.h"
#include "atomiverse/cpu/pw/pw_davidson.h"
#include "atomiverse/cpu/pw/pw_hamiltonian.h"
#include "atomiverse/cpu/pw/pw_lda.h"
#include "atomiverse/cpu/pw/pw_poisson.h"
#include "atomiverse/cpu/pw/pw_veff.h"
#include "atomiverse/cpu/pw/pw_workspace.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>

namespace ATV {
namespace cpu {
namespace pw {

PWResult pw_scf(const cpu::Frame& frame, const std::string& upf_path,
                double ecut_ry, int k_mesh, double tol, int max_iter) {
    PWResult result;

    internal::PWProblem prob;
    internal::setup_pw_problem(frame, upf_path, ecut_ry, k_mesh, prob);
    PWWorkspace& ws = prob.ws;
    auto& gv0 = prob.gv0;
    auto& gv_k = prob.gv_k;
    int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    int n_bands = prob.n_bands;
    Real vol_bohr = prob.vol_bohr;

    // === Context ===
    cpu::Context ctx;
    cpu::Frame mutable_frame = frame;
    cpu::BuildContext bctx{ctx, mutable_frame, &ws};

    // === SCF loop ===
    // Components are built once and launched per iteration (lifecycle
    // rule). PW_Convergence in particular carries its warmup counter
    // across iterations — constructing it inside the loop (as the code
    // originally did) resets the counter and the check never fires.
    PW_Poisson po;
    po.build(bctx);
    PW_LDA lda;
    lda.build(bctx);
    PW_Veff v;
    v.build(bctx);
    PW_Convergence conv{Real(tol), 2};
    conv.build(bctx);

    for (int iter = 0; iter < max_iter; ++iter) {
        po.launch({ctx, mutable_frame, 0, iter});
        lda.launch({ctx, mutable_frame, 0, iter});
        v.launch({ctx, mutable_frame, 0, iter});

        // Accumulate density over k-points in real space
        auto* rho_r_kpt = static_cast<Real*>(
            fftw_malloc(sizeof(Real) * static_cast<size_t>(ws.nr)));
        for (int i = 0; i < ws.nr; ++i)
            rho_r_kpt[i] = Real(0);

        Real e_band = 0;
        for (int ik = 0; ik < n_kpts; ++ik) {
            int n_g_k = gv_k[ik].count();
            auto* saved_gvec = ws.gvec;
            ws.gvec = &gv_k[ik];
            ws.n_g = n_g_k;
            // Nonlocal/wfc tables live on q = k+G: bind this k's tables.
            if (!prob.kb_proj_k.empty())
                ws.kb_proj = prob.kb_proj_k[ik];
            if (!prob.atom_wfc_k.empty())
                ws.atom_wfc = prob.atom_wfc_k[ik];
            // T_G is resized + filled by PW_Hamiltonian::build below; the
            // old pre-fill here wrote past the Γ-sphere-sized buffer.

            PW_Hamiltonian h;
            h.build(bctx);
            h.launch({ctx, mutable_frame, 1, iter});
            PW_Davidson d{n_bands, Real(tol) * Real(0.1), 40};
            d.build(bctx);
            d.launch({ctx, mutable_frame, 1, iter});

            Real weight = prob.kg.weights[ik];
            for (int b = 0; b < n_bands; ++b)
                e_band += weight * ws.eigenvalues[b] * Real(2);

            // Accumulate density from this k-point in real space
            if (ws.nr > 0 && ws.fft) {
                int nb =
                    std::min(n_bands, static_cast<int>(ws.eigenvalues.size()));
                int ng_k = ws.n_g;
                int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
                auto& gvk = *ws.gvec;
                auto* psi_cc = static_cast<fftw_complex*>(fftw_malloc(
                    sizeof(fftw_complex) * static_cast<size_t>(ws.nr)));
                for (int b = 0; b < nb; ++b) {
                    for (int i = 0; i < ws.nr; ++i)
                        psi_cc[i][0] = psi_cc[i][1] = 0;
                    for (int ig = 0; ig < ng_k; ++ig) {
                        int ix = ((gvk.h[ig] % nx) + nx) % nx;
                        int iy = ((gvk.k[ig] % ny) + ny) % ny;
                        int iz = ((gvk.l[ig] % nz) + nz) % nz;
                        int idx = ix + nx * (iy + ny * iz);
                        if (idx >= 0 && idx < ws.nr) {
                            auto val = ws.wfc(b, ig);
                            psi_cc[idx][0] = static_cast<double>(val.real());
                            psi_cc[idx][1] = static_cast<double>(val.imag());
                        }
                    }
                    ws.fft->c2c_bwd(psi_cc);
                    for (int i = 0; i < ws.nr; ++i) {
                        Real psi_sq =
                            static_cast<Real>(psi_cc[i][0] * psi_cc[i][0]
                                              + psi_cc[i][1] * psi_cc[i][1]);
                        // ψ_phys(r) = ψ_fft(r)/√Ω, so |ψ|²_phys = |ψ_fft|²/Ω
                        rho_r_kpt[i] += weight * Real(2) * psi_sq / ws.omega;
                    }
                }
                fftw_free(psi_cc);
            }

            ws.gvec = saved_gvec;
            ws.n_g = gv0.count();
        }

        // FFT total density to G-space; slow mixing preserves atomic density
        if (ws.nr > 0 && ws.fft) {
            int nc = ws.nc;
            auto* rho_g_full = static_cast<fftw_complex*>(
                fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nc)));
            ws.fft->r2c(rho_r_kpt, rho_g_full);

            for (int i = 0; i < ws.n_g; ++i)
                ws.rho_G_old[i] = ws.rho_G[i];

            int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz,
                nxc = nx / 2 + 1;
            Real inv_nr = Real(1) / Real(ws.nr);
            for (int ig = 0; ig < ws.n_g; ++ig) {
                int h = gv0.h[ig], k = gv0.k[ig], l = gv0.l[ig];
                int kx = ((h % nx) + nx) % nx, ky = ((k % ny) + ny) % ny,
                    kz = ((l % nz) + nz) % nz;
                bool conj = false;
                if (kx > nx / 2) {
                    kx = nx - kx;
                    ky = (ny - ky) % ny;
                    kz = (nz - kz) % nz;
                    conj = true;
                }
                int idx = kx + nxc * (ky + ny * kz);
                if (idx >= 0 && idx < nc) {
                    ws.rho_G[ig] = std::complex<Real>(
                        static_cast<Real>(rho_g_full[idx][0]) * inv_nr,
                        static_cast<Real>(conj ? -rho_g_full[idx][1]
                                               : rho_g_full[idx][1])
                            * inv_nr);
                } else
                    ws.rho_G[ig] = std::complex<Real>(0, 0);
            }
            fftw_free(rho_g_full);

            internal::pulay_mix(ws);
        }
        fftw_free(rho_r_kpt);

        conv.launch({ctx, mutable_frame, 3, iter});

        result.iterations = iter + 1;

        Real e_total_ha = internal::report_iteration_energy(
            ws, e_band, prob.ew_total, vol_bohr, iter, n_bands);
        result.energy = static_cast<double>(e_total_ha) * 2.0;

        if (ws.converged) {
            result.converged = true;
            std::cout << "\n>>> Converged in " << (iter + 1) << " iters"
                      << std::endl;
            break;
        }
    }

    if (!result.converged)
        std::cout << "\n>>> NOT converged after " << max_iter << " iters"
                  << std::endl;

    return result;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
