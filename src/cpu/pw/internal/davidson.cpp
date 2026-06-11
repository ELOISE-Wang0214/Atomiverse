/// @file cpu/pw/internal/davidson.cpp
/// @brief Iterative block Davidson with FFT-based H·ψ and KB non-local potential.
#include "atomiverse/cpu/pw/internal/davidson.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

namespace {

// ── Real spherical harmonics Y_lm(Ĝ) ──

Real y00(Real, Real, Real) { return Real(1) / std::sqrt(Real(4.0 * M_PI)); }

void y1m(Real x, Real y, Real z, Real r, Real* out) {
    Real n = std::sqrt(Real(3) / (Real(4.0 * M_PI)));
    Real ir = (r > Real(1e-16)) ? Real(1) / r : Real(0);
    out[0] = n * y * ir; // m=-1 (p_y)
    out[1] = n * z * ir; // m=0  (p_z)
    out[2] = n * x * ir; // m=+1 (p_x)
}

void y2m(Real x, Real y, Real z, Real r, Real* out) {
    Real ir2 = (r > Real(1e-16)) ? Real(1) / (r * r) : Real(0);
    out[0] =
        std::sqrt(Real(15) / (Real(4.0 * M_PI))) * x * y * ir2; // m=-2 (d_xy)
    out[1] =
        std::sqrt(Real(15) / (Real(4.0 * M_PI))) * y * z * ir2; // m=-1 (d_yz)
    out[2] = std::sqrt(Real(5) / (Real(16.0 * M_PI)))
             * (Real(3) * z * z * ir2 - Real(1)); // m=0 (d_z²)
    out[3] =
        std::sqrt(Real(15) / (Real(4.0 * M_PI))) * x * z * ir2; // m=+1 (d_xz)
    out[4] = std::sqrt(Real(15) / (Real(16.0 * M_PI))) * (x * x - y * y)
             * ir2; // m=+2 (d_x²-y²)
}

int n_m_for_l(int l) { return 2 * l + 1; }

void compute_ylm(Real x, Real y, Real z, Real r, int l, Real* out) {
    switch (l) {
    case 0:
        out[0] = y00(x, y, z);
        break;
    case 1:
        y1m(x, y, z, r, out);
        break;
    case 2:
        y2m(x, y, z, r, out);
        break;
    default:
        for (int m = 0; m < n_m_for_l(l); ++m)
            out[m] = 0;
        break;
    }
}

/// Deterministic LCG in [-0.5, 0.5) for the atomic+random initial-guess
/// perturbation. Seeded per subspace-vector index so the result is
/// independent of which MPI rank handles which k-point (std::rand's
/// global stream would make the wavefunction guess rank-layout
/// dependent).
struct TinyRng {
    std::uint32_t state;
    explicit TinyRng(std::uint32_t seed) : state(seed * 2654435761u + 1u) {}
    Real next() {
        state = state * 1664525u + 1013904223u;
        return Real(state) / Real(4294967296.0) - Real(0.5);
    }
};

// ── Subspace Jacobi diagonalization (for small matrices) ──

void jacobi_eigh(int n, Real* A, Real* eig, Real* vecs, int max_sweep = 20) {
    // Initialize eigenvectors to identity
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j)
            vecs[i * n + j] = (i == j) ? Real(1) : Real(0);
        eig[i] = A[i * n + i];
    }

    for (int sweep = 0; sweep < max_sweep; ++sweep) {
        Real max_off = 0;
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                max_off = std::max(max_off, std::abs(A[i * n + j]));
        if (max_off < Real(1e-14))
            break;

        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                Real apq = A[p * n + q];
                if (std::abs(apq) < Real(1e-14))
                    continue;
                Real app = A[p * n + p], aqq = A[q * n + q];
                Real theta = Real(0.5) * std::atan2(Real(2) * apq, aqq - app);
                Real c = std::cos(theta), s = std::sin(theta);

                // Rotate A
                A[p * n + p] =
                    c * c * app + s * s * aqq - Real(2) * c * s * apq;
                A[q * n + q] =
                    s * s * app + c * c * aqq + Real(2) * c * s * apq;
                A[p * n + q] = A[q * n + p] = Real(0);
                for (int r = 0; r < n; ++r) {
                    if (r == p || r == q)
                        continue;
                    Real arp = A[r * n + p], arq = A[r * n + q];
                    A[r * n + p] = A[p * n + r] = c * arp - s * arq;
                    A[r * n + q] = A[q * n + r] = s * arp + c * arq;
                }
                // Rotate eigenvectors
                for (int r = 0; r < n; ++r) {
                    Real vrp = vecs[r * n + p], vrq = vecs[r * n + q];
                    vecs[r * n + p] = c * vrp - s * vrq;
                    vecs[r * n + q] = s * vrp + c * vrq;
                }
            }
        }
    }
    for (int i = 0; i < n; ++i)
        eig[i] = A[i * n + i];
}

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

void DavidsonSolver::prepare(const PWWorkspace& ws) {
    if (ws.nr > 0)
        work_buf_.resize(static_cast<std::size_t>(ws.nr));
}

void DavidsonSolver::step(PWWorkspace& ws) {
    int ng = ws.n_g;
    int nb = n_bands_;
    int nr = ws.nr;

    ws.wfc.resize(
        {static_cast<std::size_t>(nb), static_cast<std::size_t>(ng)});
    ws.eigenvalues.resize({static_cast<std::size_t>(nb)});

    // FFT work buffer for H·ψ — pre-allocated in prepare(); guard against
    // direct low-level use that skipped prepare().
    if (work_buf_.size() < static_cast<std::size_t>(nr))
        work_buf_.resize(static_cast<std::size_t>(nr));
    auto* work_buf = work_buf_.data();

    // ── Initialize subspace: lowest T(G) plane waves ──
    struct DiagElem {
        Real val;
        int idx;
    };
    std::vector<DiagElem> diag(ng);
    for (int i = 0; i < ng; ++i)
        diag[i] = {ws.T_G[i], i};
    std::sort(
        diag.begin(), diag.end(),
        [](const DiagElem& a, const DiagElem& b) { return a.val < b.val; });

    // Count atomic wavefunction states: n_atoms × Σ (2l+1)
    int n_atom_states = 0;
    for (size_t ia = 0; ia < ws.atom_wfc.size(); ++ia)
        for (size_t iw = 0; iw < ws.atom_wfc[ia].size(); ++iw)
            n_atom_states += 2 * ws.atom_wfc[ia][iw].l + 1;

    int n_sub;
    if (n_atom_states >= nb) {
        n_sub = n_atom_states; // QE: start from atomic wavefunctions
    } else {
        n_sub = std::min(ng, nb * 4);
    }
    std::vector<std::complex<Real>> sub_vectors(n_sub * ng);
    std::vector<std::complex<Real>> H_sub_vectors(n_sub * ng);
    std::vector<bool> H_computed(n_sub, false);

    if (n_atom_states >= nb) {
        // Initialize subspace with atomic wavefunctions (QE atomic_wfc_mod.f90)
        // ψ_{nlm}(G) = i^l · Y_lm(Ĝ) · chi_G(G) · exp(-iG·R)
        int j = 0;
        auto& gv = *ws.gvec;
        for (size_t ia = 0; ia < ws.atom_wfc.size() && j < n_sub; ++ia) {
            Real ax = ws.atom_pos_bohr[ia][0], ay = ws.atom_pos_bohr[ia][1],
                 az = ws.atom_pos_bohr[ia][2];
            for (size_t iw = 0; iw < ws.atom_wfc[ia].size() && j < n_sub;
                 ++iw) {
                int l = ws.atom_wfc[ia][iw].l;
                int nm = 2 * l + 1;
                for (int m = 0; m < nm && j < n_sub; ++m) {
                    Real norm2 = 0;
                    for (int ig = 0; ig < ng; ++ig) {
                        Real Gx = gv.g[ig][0], Gy = gv.g[ig][1],
                             Gz = gv.g[ig][2];
                        Real Gr = std::sqrt(Gx * Gx + Gy * Gy + Gz * Gz);
                        Real chi = ws.atom_wfc[ia][iw].chi_G[ig];
                        if (std::abs(chi) < Real(1e-30)) {
                            sub_vectors[j * ng + ig] = std::complex<Real>(0, 0);
                            continue;
                        }
                        // i^l factor (QE: lphase = (0,1)^l)
                        std::complex<Real> il(1, 0);
                        if (l == 1)
                            il = std::complex<Real>(0, 1);
                        else if (l == 2)
                            il = std::complex<Real>(-1, 0);
                        // Y_lm(Ĝ)
                        std::vector<Real> ylm(nm);
                        compute_ylm(Gx, Gy, Gz, Gr, l, ylm.data());
                        // Structure factor S(G) = exp(-iG·R)
                        Real phase = Gx * ax + Gy * ay + Gz * az;
                        std::complex<Real> S(std::cos(phase), -std::sin(phase));
                        sub_vectors[j * ng + ig] = il * ylm[m] * chi * S;
                        norm2 += std::norm(sub_vectors[j * ng + ig]);
                    }
                    // QE startingwfc='atomic+random': add random perturbation
                    // (deterministic per state index — see TinyRng).
                    TinyRng rng(static_cast<std::uint32_t>(j + 1));
                    for (int ig = 0; ig < ng; ++ig) {
                        Real rnd = rng.next();
                        sub_vectors[j * ng + ig] +=
                            std::complex<Real>(rnd * Real(1e-4), Real(0));
                    }
                    ++j;
                }
            }
        }
        n_sub = j; // actual count
    } else {
        // Fallback: lowest T(G) plane waves
        for (int j = 0; j < n_sub; ++j) {
            int g_idx = diag[j].idx;
            for (int ig = 0; ig < ng; ++ig)
                sub_vectors[j * ng + ig] = (ig == g_idx)
                                               ? std::complex<Real>(1, 0)
                                               : std::complex<Real>(0, 0);
        }
    }

    // ── Davidson iteration ──
    int max_david = std::min(max_iter_, 20);
    Real tol_sq = Real(1e-16); // very tight, always add residuals

    for (int idavid = 0; idavid < max_david; ++idavid) {
        // Compute H·v_j for new subspace vectors
        for (int j = 0; j < n_sub; ++j) {
            if (H_computed[j])
                continue;
            apply_hamiltonian(ws, &sub_vectors[j * ng],
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
        jacobi_eigh(n_sub, Hs.data(), sub_eig.data(), sub_vecs.data());

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
            apply_hamiltonian(ws, &ritz[b * ng], Hpsi.data(), work_buf);

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
            ws.eigenvalues[b] = sub_eig[b];

            // Copy Ritz vector to wfc
            for (int ig = 0; ig < ng; ++ig)
                ws.wfc(b, ig) = ritz[b * ng + ig];
        }

        // Always add residuals to build crystal field (QE approach)

        // Add preconditioned residuals for unconverged bands
        if (idavid < max_david - 1 && n_sub < ng - 1) {
            for (int b = 0; b < nb && n_sub < ng; ++b) {
                // Compute residual R = H·ψ - ε·ψ
                std::vector<std::complex<Real>> Hpsi(ng);
                apply_hamiltonian(ws, &ritz[b * ng], Hpsi.data(), work_buf);

                std::vector<std::complex<Real>> R(ng);
                Real res_norm = 0;
                for (int ig = 0; ig < ng; ++ig) {
                    R[ig] = Hpsi[ig] - sub_eig[b] * ritz[b * ng + ig];
                    res_norm += std::norm(R[ig]);
                }

                if (res_norm < tol_sq)
                    continue; // already converged

                // Precondition: δ(G) = -R(G) / (T(G) - ε)
                std::vector<std::complex<Real>> delta(ng);
                Real eps = sub_eig[b];
                for (int ig = 0; ig < ng; ++ig) {
                    Real denom = ws.T_G[ig] - eps;
                    // Damp near-degenerate contributions
                    Real damp = std::abs(denom);
                    if (damp < Real(1e-6))
                        denom = (denom >= 0 ? Real(1e-6) : Real(-1e-6));
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

                // Normalize
                Real delta_norm = 0;
                for (int ig = 0; ig < ng; ++ig)
                    delta_norm += std::norm(delta[ig]);
                if (delta_norm < Real(1e-30))
                    continue;
                Real inv_norm = Real(1) / std::sqrt(delta_norm);
                for (int ig = 0; ig < ng; ++ig)
                    delta[ig] *= inv_norm;

                // Add to subspace — grow the storage BEFORE writing the
                // new column. (The historical PW_Davidson wrote slot
                // n_sub first and resized afterwards: a heap-buffer
                // overflow of one full column per expansion, latent
                // until the allocator layout exposed it.)
                int j = n_sub++;
                sub_vectors.resize(static_cast<std::size_t>(n_sub) * ng);
                H_sub_vectors.resize(static_cast<std::size_t>(n_sub) * ng);
                for (int ig = 0; ig < ng; ++ig)
                    sub_vectors[j * ng + ig] = delta[ig];
                H_computed.push_back(false);
            }
        }
    }
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
