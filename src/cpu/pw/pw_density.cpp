/// @file cpu/pw/pw_density.cpp
/// @brief Build charge density from wavefunctions using FFT + Pulay mixing.
#include "atomiverse/cpu/pw/pw_density.h"
#include <cmath>
#include <complex>
#include <cstring>
#include <fftw3.h>

namespace ATV {
namespace cpu {
namespace pw {

void PW_Density::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
    ws_->pulay_depth = pulay_depth_;
    ws_->mix_beta = mix_beta_;
}

void PW_Density::launch(const cpu::LaunchContext&) {
    int ng = ws_->n_g;
    int nb = static_cast<int>(ws_->eigenvalues.size());

    for (int i = 0; i < ng; ++i)
        ws_->rho_G_old[i] = ws_->rho_G[i];

    if (nb > 0 && ws_->wfc.size() > 0 && ws_->fft && ws_->nr > 0) {
        int nr = ws_->nr;
        int nc = ws_->nc;
        int nx = ws_->fft_nx, ny = ws_->fft_ny, nz = ws_->fft_nz;

        auto& wfc = ws_->wfc;
        int ng_wfc = static_cast<int>(wfc.shape(1));
        auto& gv = *ws_->gvec;

        // Allocate full complex grid for c2c transforms
        auto* psi_cc = static_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nr)));
        auto* rho_r = static_cast<Real*>(
            fftw_malloc(sizeof(Real) * static_cast<size_t>(nr)));
        auto* rho_g_half = static_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nc)));

        for (int i = 0; i < nr; ++i) rho_r[i] = Real(0);

        // Build map: (h,k,l) on wfc grid → flat index on density c2c grid
        // The c2c grid uses a FULL complex array (nx*ny*nz), not half-complex.
        // Index: ix + nx*(iy + ny*iz) with ix=0..nx-1, iy=0..ny-1, iz=0..nz-1
        for (int b = 0; b < nb; ++b) {
            for (int i = 0; i < nr; ++i)
                psi_cc[i][0] = psi_cc[i][1] = 0;

            int count = 0;
            for (int ig = 0; ig < ng_wfc; ++ig) {
                int h = gv.h[ig], k = gv.k[ig], l = gv.l[ig];
                int ix = ((h % nx) + nx) % nx;
                int iy = ((k % ny) + ny) % ny;
                int iz = ((l % nz) + nz) % nz;
                int idx = ix + nx * (iy + ny * iz);
                if (idx >= 0 && idx < nr) {
                    auto val = wfc(b, ig);
                    psi_cc[idx][0] = static_cast<double>(val.real());
                    psi_cc[idx][1] = static_cast<double>(val.imag());
                    ++count;
                }
            }

            // c2c backward: G → r. FFTW normalizes by 1/N (N=nr).
            ws_->fft->c2c_bwd(psi_cc);

            // Accumulate: ρ(r) += 2 * |ψ(r)|²
            // ψ_fft(r) = (1/N) Σ_G ψ(G) exp(iGr)
            // ψ_phys(r) = (1/√Ω) Σ_G ψ(G) exp(iGr) = √Ω/N · ψ_fft(r)
            // |ψ_phys(r)|² = Ω/N² · |ψ_fft(r)|²
            // ρ_phys(r) = 2 Σ_b |ψ_phys(r)|² = 2 Ω/N² · Σ_b |ψ_fft(r)|²
            // We'll fix normalization at the end.
            for (int i = 0; i < nr; ++i) {
                Real v = static_cast<Real>(psi_cc[i][0] * psi_cc[i][0]
                                         + psi_cc[i][1] * psi_cc[i][1]);
                rho_r[i] += Real(2) * v;
            }
        }

        // r2c: ρ(r) → ρ(G) on half-complex grid.
        // r2c returns unnormalized sum: C'[G] = Σ_r ρ(r) exp(-iGr)
        ws_->fft->r2c(rho_r, rho_g_half);

        // Extract ρ(G) at G-vector sphere positions
        // r2c half-complex format: kx=0..nx/2
        int nxc = nx/2 + 1;
        Real inv_nr = Real(1) / static_cast<Real>(nr);
        for (int ig = 0; ig < ng; ++ig) {
            int h = gv.h[ig], k = gv.k[ig], l = gv.l[ig];
            int kx = ((h % nx) + nx) % nx;
            int ky = ((k % ny) + ny) % ny;
            int kz = ((l % nz) + nz) % nz;
            bool conj = false;
            if (kx > nx/2) {
                kx = nx - kx;
                ky = (ny - ky) % ny;
                kz = (nz - kz) % nz;
                conj = true;
            }
            int idx = kx + nxc * (ky + ny * kz);
            if (idx >= 0 && idx < nc) {
                Real re = static_cast<Real>(rho_g_half[idx][0]) * inv_nr;
                Real im = static_cast<Real>(rho_g_half[idx][1]) * inv_nr;
                if (conj) im = -im;
                ws_->rho_G[ig] = std::complex<Real>(re, im);
            } else {
                ws_->rho_G[ig] = std::complex<Real>(0, 0);
            }
        }

        // Renormalize: enforce ρ(G=0) = n_elec / Ω
        Real current_rho0 = ws_->rho_G[0].real();
        Real target_rho0 = ws_->rho_G_old[0].real();
        if (std::abs(current_rho0) > Real(1e-16)) {
            Real scale = target_rho0 / current_rho0;
            for (int ig = 0; ig < ng; ++ig)
                ws_->rho_G[ig] *= scale;
        }

        fftw_free(psi_cc);
        fftw_free(rho_r);
        fftw_free(rho_g_half);

        // === USPP Augmentation charge ===
        // Δρ(G) = Σ_a S_a(G) * F(G) * Σ_ij Q_ij Σ_b f_b c_{b,i}^{a*} c_{b,j}^a
        if (!ws_->q_ij.empty() && !ws_->kb_proj.empty()
            && ws_->aug_form_factor_G.size() > 0) {
            int n_proj = ws_->n_proj_total;
            Real E_aug = 0;

            for (size_t ia = 0; ia < ws_->kb_proj.size(); ++ia) {
                auto& atom_proj = ws_->kb_proj[ia];
                Real ax = ws_->atom_pos_bohr[ia][0];
                Real ay = ws_->atom_pos_bohr[ia][1];
                Real az = ws_->atom_pos_bohr[ia][2];
                int np = static_cast<int>(atom_proj.size());

                // Compute overlap sums per band + projector
                // S_ij = Σ_b f_b c_i^{a*} c_j^a
                std::vector<std::complex<Real>> S_ij(np * np, std::complex<Real>(0,0));
                for (int b = 0; b < nb; ++b) {
                    std::vector<std::complex<Real>> c(np, std::complex<Real>(0,0));
                    // Reconstruct spherically-averaged beta for augmentation:
                    // β_avg(G) = √(4π)/√Ω * f_l(G) * exp(-iG·R)
                    Real prefac = std::sqrt(Real(4.0*M_PI)/ws_->omega);
                    for (int ip = 0; ip < np; ++ip) {
                        for (int ig = 0; ig < ng_wfc; ++ig) {
                            Real Gx=gv.g[ig][0], Gy=gv.g[ig][1], Gz=gv.g[ig][2];
                            Real phase = Gx*ax + Gy*ay + Gz*az;
                            std::complex<Real> beta_avg(
                                prefac * atom_proj[ip].f_G[ig] * std::cos(phase),
                                -prefac * atom_proj[ip].f_G[ig] * std::sin(phase));
                            c[ip] += std::conj(beta_avg) * wfc(b, ig);
                        }
                    }
                    for (int ip = 0; ip < np; ++ip)
                        for (int jp = 0; jp < np; ++jp)
                            S_ij[ip * np + jp] += Real(2) * std::conj(c[ip]) * c[jp];
                }

                // Compute S = Σ_ij Q_ij S_ij for this atom
                std::complex<Real> S_atom(0, 0);
                for (int ip = 0; ip < np; ++ip) {
                    for (int jp = 0; jp < np; ++jp) {
                        int q_idx = ip * n_proj + jp;
                        if (q_idx < static_cast<int>(ws_->q_ij.size()))
                            S_atom += ws_->q_ij[q_idx] * S_ij[ip * np + jp];
                    }
                }
                E_aug += S_atom.real(); // augmentation energy contribution

                // Add Δρ_aug(G) = S_a(G) * F(G) * S_atom
                for (int ig = 0; ig < ng; ++ig) {
                    if (std::abs(ws_->aug_form_factor_G[ig]) < Real(1e-30)) continue;
                    Real Gx = gv.g[ig][0], Gy = gv.g[ig][1], Gz = gv.g[ig][2];
                    Real phase = Gx*ax + Gy*ay + Gz*az;
                    std::complex<Real> Sa(std::cos(phase), -std::sin(phase));
                    ws_->rho_G[ig] += ws_->aug_form_factor_G[ig] * Sa * S_atom;
                }
            }
            ws_->E_aug = E_aug;
        }
    }

    // Pulay mixing
    xt::xarray<std::complex<Real>> res = ws_->rho_G - ws_->rho_G_old;

    if (iter_ == 0) {
        for (int i = 0; i < ng; ++i)
            ws_->rho_G[i] = ws_->rho_G_old[i] + ws_->mix_beta * res[i];
    } else {
        ws_->pulay_res.push_back(res);
        ws_->pulay_rho.push_back(ws_->rho_G_old);
        if (static_cast<int>(ws_->pulay_res.size()) > ws_->pulay_depth) {
            ws_->pulay_res.erase(ws_->pulay_res.begin());
            ws_->pulay_rho.erase(ws_->pulay_rho.begin());
        }
        for (int i = 0; i < ng; ++i)
            ws_->rho_G[i] = ws_->rho_G_old[i] + ws_->mix_beta * res[i];
    }
    ++iter_;
}

} // namespace pw
} // namespace cpu
} // namespace ATV
