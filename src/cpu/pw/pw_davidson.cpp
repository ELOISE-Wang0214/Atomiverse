/// @file cpu/pw/pw_davidson.cpp
/// @brief Iterative block Davidson with FFT-based H·ψ and KB non-local potential.
///
/// Subspace utilities (Y_lm, Jacobi eigensolver, initial subspace) live in
/// internal/pw_subspace.cpp, shared with the CUDA band solver.
#include "atomiverse/cpu/pw/pw_davidson.h"
#include "atomiverse/cpu/pw/internal/pw_subspace.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <fftw3.h>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {

using internal::compute_ylm;
using internal::n_m_for_l;

namespace {

// ── Pad G-sphere vector to full c2c FFT grid ──

void pad_to_fft(const PWWorkspace& ws, const std::complex<Real>* psi_g,
                fftw_complex* psi_fft) {
    int nr = ws.nr;
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    auto& gv = *ws.gvec;
    for (int i = 0; i < nr; ++i)
        psi_fft[i][0] = psi_fft[i][1] = 0;
    for (int ig = 0; ig < ws.n_g; ++ig) {
        int ix = ((gv.h[ig] % nx) + nx) % nx;
        int iy = ((gv.k[ig] % ny) + ny) % ny;
        int iz = ((gv.l[ig] % nz) + nz) % nz;
        int idx = ix + nx * (iy + ny * iz);
        if (idx >= 0 && idx < nr) {
            psi_fft[idx][0] = static_cast<double>(psi_g[ig].real());
            psi_fft[idx][1] = static_cast<double>(psi_g[ig].imag());
        }
    }
}

// ── Extract from FFT grid back to G-sphere ──

void extract_from_fft(const PWWorkspace& ws, const fftw_complex* psi_fft,
                      std::complex<Real>* psi_g) {
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    auto& gv = *ws.gvec;
    for (int ig = 0; ig < ws.n_g; ++ig) {
        int ix = ((gv.h[ig] % nx) + nx) % nx;
        int iy = ((gv.k[ig] % ny) + ny) % ny;
        int iz = ((gv.l[ig] % nz) + nz) % nz;
        int idx = ix + nx * (iy + ny * iz);
        psi_g[ig] = std::complex<Real>(static_cast<Real>(psi_fft[idx][0]),
                                       static_cast<Real>(psi_fft[idx][1]));
    }
}

// ── Apply Hamiltonian: (Hψ)(G) = T(G)ψ(G) + FFT[V_eff(r)·IFFT[ψ(G)]] + V_NL·ψ ──

void apply_hamiltonian(const PWWorkspace& ws, const std::complex<Real>* psi_g,
                       std::complex<Real>* Hpsi_g, fftw_complex* work_buf) {
    int ng = ws.n_g, nr = ws.nr;

    // Step 1: IFFT ψ(G) → ψ(r)  (c2c backward)
    pad_to_fft(ws, psi_g, work_buf);
    ws.fft->c2c_bwd(work_buf); // now work_buf = ψ(r)

    // Step 2: ψ(r) *= V_eff(r)
    for (int i = 0; i < nr; ++i) {
        Real V = ws.V_eff_r[i];
        work_buf[i][0] *= V;
        work_buf[i][1] *= V;
    }

    // Step 3: FFT V_eff(r)·ψ(r) → (Vψ)(G)  (c2c forward)
    ws.fft->c2c_fwd(work_buf); // now work_buf = FFT[V_eff·ψ]

    // QE: tgwave_r2g divides by N. FFTW c2c_fwd does not → divide by N.
    // This gives the correct off-diagonal V_eff for density development.
    extract_from_fft(ws, work_buf, Hpsi_g);
    Real inv_N = Real(1) / Real(nr);
    for (int ig = 0; ig < ng; ++ig)
        Hpsi_g[ig] = Hpsi_g[ig] * inv_N + ws.T_G[ig] * psi_g[ig];

    // Step 5: Non-local KB potential: V_NL = Σ_a Σ_lm D_l β_lm^a ⟨β_lm^a|ψ⟩
    Real four_pi_over_sqrt_omega = Real(4.0 * M_PI) / std::sqrt(ws.omega);
    for (size_t ia = 0; ia < ws.kb_proj.size(); ++ia) {
        auto& atom_proj = ws.kb_proj[ia];
        Real ax = ws.atom_pos_bohr[ia][0], ay = ws.atom_pos_bohr[ia][1],
             az = ws.atom_pos_bohr[ia][2];
        auto& gv = *ws.gvec;

        for (size_t ip = 0; ip < atom_proj.size(); ++ip) {
            int l = atom_proj[ip].l;
            Real D = atom_proj[ip].d;
            if (std::abs(D) < Real(1e-16))
                continue;
            int nm = n_m_for_l(l);

            // Compute c_lm = ⟨β_lm^a|ψ⟩ and apply V_NL in one pass
            // β_lm(G) = (4π/√Ω)·(-i)^l·Y_lm(Ĝ)·f_l(G)·exp(-iG·R)
            std::vector<std::complex<Real>> c_lm(nm, std::complex<Real>(0, 0));
            for (int ig = 0; ig < ng; ++ig) {
                Real Gx = gv.g[ig][0], Gy = gv.g[ig][1], Gz = gv.g[ig][2];
                Real Gr = std::sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
                Real fl_val = atom_proj[ip].f_G[ig];
                if (std::abs(fl_val) < Real(1e-30))
                    continue;
                Real phase = Gx * ax + Gy * ay + Gz * az;
                Real cos_p = std::cos(phase), sin_p = std::sin(phase);

                std::vector<Real> ylm(nm);
                compute_ylm(Gx, Gy, Gz, Gr, l, ylm.data());

                std::complex<Real> im_l(1, 0);
                if (l == 1)
                    im_l = std::complex<Real>(0, -1);
                else if (l == 2)
                    im_l = std::complex<Real>(-1, 0);

                for (int m = 0; m < nm; ++m) {
                    std::complex<Real> beta_lm =
                        four_pi_over_sqrt_omega * im_l * ylm[m] * fl_val
                        * std::complex<Real>(cos_p, -sin_p);
                    c_lm[m] += std::conj(beta_lm) * psi_g[ig];
                }
            }

            // Apply V_NL: Hψ(G) += Σ_m D·β_lm(G)·c_lm
            for (int ig = 0; ig < ng; ++ig) {
                Real Gx = gv.g[ig][0], Gy = gv.g[ig][1], Gz = gv.g[ig][2];
                Real Gr = std::sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
                Real fl_val = atom_proj[ip].f_G[ig];
                if (std::abs(fl_val) < Real(1e-30))
                    continue;
                Real phase = Gx * ax + Gy * ay + Gz * az;
                Real cos_p = std::cos(phase), sin_p = std::sin(phase);

                std::vector<Real> ylm(nm);
                compute_ylm(Gx, Gy, Gz, Gr, l, ylm.data());

                std::complex<Real> im_l(1, 0);
                if (l == 1)
                    im_l = std::complex<Real>(0, -1);
                else if (l == 2)
                    im_l = std::complex<Real>(-1, 0);

                for (int m = 0; m < nm; ++m) {
                    std::complex<Real> beta_lm =
                        four_pi_over_sqrt_omega * im_l * ylm[m] * fl_val
                        * std::complex<Real>(cos_p, -sin_p);
                    Hpsi_g[ig] += D * beta_lm * c_lm[m];
                }
            }
        }
    }
}

} // namespace

// ── Public interface ──

void PW_Davidson::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_Davidson::launch(const cpu::LaunchContext&) {
    int ng = ws_->n_g;
    int nb = n_bands_;
    int nr = ws_->nr;

    ws_->wfc.resize(
        {static_cast<std::size_t>(nb), static_cast<std::size_t>(ng)});
    ws_->eigenvalues.resize({static_cast<std::size_t>(nb)});

    // Allocate work buffer for FFT-based H·ψ
    auto* work_buf = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nr)));

    // ── Initialize subspace (atomic wavefunctions or lowest-T(G) PWs) ──
    std::vector<std::complex<Real>> sub_vectors;
    int n_sub = internal::init_subspace(*ws_, nb, sub_vectors);
    std::vector<std::complex<Real>> H_sub_vectors(n_sub * ng);
    std::vector<bool> H_computed(n_sub, false);

    // ── Davidson iteration ──
    int max_david = std::min(max_iter_, 20);
    Real tol_sq = Real(1e-16); // very tight, always add residuals

    for (int idavid = 0; idavid < max_david; ++idavid) {
        // Compute H·v_j for new subspace vectors
        for (int j = 0; j < n_sub; ++j) {
            if (H_computed[j])
                continue;
            apply_hamiltonian(*ws_, &sub_vectors[j * ng],
                              &H_sub_vectors[j * ng], work_buf);
            H_computed[j] = true;
        }

        // Build subspace Hamiltonian H_sub[μ,ν] = ⟨v_μ|H|v_ν⟩
        std::vector<Real> Hs(n_sub * n_sub, Real(0));
        for (int mu = 0; mu < n_sub; ++mu) {
            for (int nu = mu; nu < n_sub; ++nu) {
                std::complex<Real> sum(0, 0);
                for (int ig = 0; ig < ng; ++ig)
                    sum += std::conj(sub_vectors[mu * ng + ig])
                           * H_sub_vectors[nu * ng + ig];
                Hs[mu * n_sub + nu] = Hs[nu * n_sub + mu] = sum.real();
            }
        }

        // Diagonalize subspace
        std::vector<Real> sub_eig(n_sub), sub_vecs(n_sub * n_sub);
        internal::jacobi_eigh(n_sub, Hs.data(), sub_eig.data(),
                              sub_vecs.data());

        // Build Ritz vectors: ψ_i = Σ_μ x_{μ,i} v_μ
        std::vector<std::complex<Real>> ritz(nb * ng);
        for (int b = 0; b < nb; ++b) {
            for (int ig = 0; ig < ng; ++ig)
                ritz[b * ng + ig] = std::complex<Real>(0, 0);
            for (int mu = 0; mu < n_sub; ++mu)
                for (int ig = 0; ig < ng; ++ig)
                    ritz[b * ng + ig] +=
                        sub_vecs[mu * n_sub + b] * sub_vectors[mu * ng + ig];
        }

        // Compute residuals and check convergence
        Real max_res = 0;
        for (int b = 0; b < nb; ++b) {
            // H·ψ_b
            std::vector<std::complex<Real>> Hpsi(ng);
            apply_hamiltonian(*ws_, &ritz[b * ng], Hpsi.data(), work_buf);

            // R = H·ψ - ε·ψ
            Real res_norm = 0;
            for (int ig = 0; ig < ng; ++ig) {
                auto R = Hpsi[ig] - sub_eig[b] * ritz[b * ng + ig];
                res_norm += std::norm(R);
            }
            res_norm = std::sqrt(res_norm);
            if (res_norm > max_res)
                max_res = res_norm;

            // Store eigenvalue
            ws_->eigenvalues[b] = sub_eig[b];

            // Copy Ritz vector to wfc
            for (int ig = 0; ig < ng; ++ig)
                ws_->wfc(b, ig) = ritz[b * ng + ig];
        }

        // Always add residuals to build crystal field (QE approach)

        // Add preconditioned residuals for unconverged bands
        if (idavid < max_david - 1 && n_sub < ng - 1) {
            int n_added = 0;
            for (int b = 0; b < nb && n_sub < ng; ++b) {
                // Compute residual R = H·ψ - ε·ψ
                std::vector<std::complex<Real>> Hpsi(ng);
                apply_hamiltonian(*ws_, &ritz[b * ng], Hpsi.data(), work_buf);

                std::vector<std::complex<Real>> R(ng);
                Real res_norm = 0;
                for (int ig = 0; ig < ng; ++ig) {
                    R[ig] = Hpsi[ig] - sub_eig[b] * ritz[b * ng + ig];
                    res_norm += std::norm(R[ig]);
                }

                if (res_norm < tol_sq)
                    continue; // already converged

                // Precondition: δ(G) = -R(G) / (T(G) - ε). The clamp keeps
                // the T(G) ≈ ε crossing from amplifying residual noise by
                // orders of magnitude (the preconditioner only needs to be
                // roughly right for Davidson to converge).
                std::vector<std::complex<Real>> delta(ng);
                Real eps = sub_eig[b];
                for (int ig = 0; ig < ng; ++ig) {
                    Real denom = ws_->T_G[ig] - eps;
                    Real damp = std::abs(denom);
                    if (damp < Real(1e-2))
                        denom = (denom >= 0 ? Real(1e-2) : Real(-1e-2));
                    delta[ig] = -R[ig] / denom;
                }

                // Orthogonalize against existing subspace
                for (int mu = 0; mu < n_sub; ++mu) {
                    std::complex<Real> ovlp(0, 0);
                    for (int ig = 0; ig < ng; ++ig)
                        ovlp +=
                            std::conj(sub_vectors[mu * ng + ig]) * delta[ig];
                    for (int ig = 0; ig < ng; ++ig)
                        delta[ig] -= ovlp * sub_vectors[mu * ng + ig];
                }

                // Normalize. A post-orthogonalization norm at rounding
                // level means the direction is already spanned —
                // normalizing it would inject pure noise into the
                // subspace and put a floor under the eigenvalues.
                Real delta_norm = 0;
                for (int ig = 0; ig < ng; ++ig)
                    delta_norm += std::norm(delta[ig]);
                if (delta_norm < Real(1e-12))
                    continue;
                Real inv_norm = Real(1) / std::sqrt(delta_norm);
                for (int ig = 0; ig < ng; ++ig)
                    delta[ig] *= inv_norm;

                // Add to subspace
                int j = n_sub++;
                sub_vectors.resize(static_cast<size_t>(n_sub) * ng);
                for (int ig = 0; ig < ng; ++ig)
                    sub_vectors[j * ng + ig] = delta[ig];
                H_computed.push_back(false);
                ++n_added;
            }
            // Resize arrays
            sub_vectors.resize(static_cast<size_t>(n_sub) * ng);
            H_sub_vectors.resize(static_cast<size_t>(n_sub) * ng);
            H_computed.resize(n_sub, false);
            if (n_added == 0)
                break; // all bands converged / spanned — done
        }
    }

    fftw_free(work_buf);
}

} // namespace pw
} // namespace cpu
} // namespace ATV
