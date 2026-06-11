/// @file cpu/pw/internal/density_builder.cpp
/// @brief Build charge density from wavefunctions using FFT + Pulay mixing.
#include "atomiverse/cpu/pw/internal/density_builder.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

void DensityBuilder::prepare(const PWWorkspace& ws) {
    if (ws.nr <= 0)
        return;
    rho_r_accum_.resize(static_cast<std::size_t>(ws.nr));
    psi_cc_.resize(static_cast<std::size_t>(ws.nr));
    rho_g_full_.resize(static_cast<std::size_t>(ws.nc));
}

void DensityBuilder::reset(PWWorkspace& ws) {
    ws.pulay_res.clear();
    ws.pulay_rho.clear();
}

void DensityBuilder::begin(const PWWorkspace& ws) {
    auto* acc = rho_r_accum_.data();
    for (int i = 0; i < ws.nr; ++i)
        acc[i] = Real(0);
}

void DensityBuilder::accumulate_kpoint(const PWWorkspace& ws, Real weight,
                                       int n_bands) {
    if (ws.nr <= 0 || !ws.fft)
        return;

    int nb = std::min(n_bands, static_cast<int>(ws.eigenvalues.size()));
    int ng_k = ws.n_g;
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    auto& gvk = *ws.gvec;
    auto* psi_cc = psi_cc_.data();
    auto* acc = rho_r_accum_.data();

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
            Real psi_sq = static_cast<Real>(psi_cc[i][0] * psi_cc[i][0]
                                            + psi_cc[i][1] * psi_cc[i][1]);
            // ψ_phys(r) = ψ_fft(r)/√Ω, so |ψ|²_phys = |ψ_fft|²/Ω
            acc[i] += weight * Real(2) * psi_sq / ws.omega;
        }
    }
}

void DensityBuilder::finish(PWWorkspace& ws) {
    if (ws.nr <= 0 || !ws.fft)
        return;

    int nc = ws.nc;
    auto* rho_g_full = rho_g_full_.data();
    ws.fft->r2c(rho_r_accum_.data(), rho_g_full);

    for (int i = 0; i < ws.n_g; ++i)
        ws.rho_G_old[i] = ws.rho_G[i];

    // Extract ρ(G) at the density G-sphere from the half-complex grid.
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz, nxc = nx / 2 + 1;
    auto& gv = *ws.gvec;
    Real inv_nr = Real(1) / Real(ws.nr);
    for (int ig = 0; ig < ws.n_g; ++ig) {
        int h = gv.h[ig], k = gv.k[ig], l = gv.l[ig];
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

    // ── Pulay/DIIS mixing (QE's Anderson/Broyden simplified) ──
    // Store residual F = ρ_out - ρ_in; keep a sliding window of the
    // newest `pulay_depth_` entries.
    xt::xarray<std::complex<Real>> F = ws.rho_G - ws.rho_G_old;
    ws.pulay_res.push_back(F);
    ws.pulay_rho.push_back(ws.rho_G_old);
    while (static_cast<int>(ws.pulay_res.size()) > pulay_depth_) {
        ws.pulay_res.erase(ws.pulay_res.begin());
        ws.pulay_rho.erase(ws.pulay_rho.begin());
    }
    int depth = static_cast<int>(ws.pulay_res.size());

    if (depth >= 2) {
        // Build B_{ij} = ⟨F_i|F_j⟩
        std::vector<Real> B(depth * depth, Real(0));
        for (int i = 0; i < depth; ++i) {
            for (int j = i; j < depth; ++j) {
                Real dot = 0;
                for (int ig = 0; ig < ws.n_g; ++ig) {
                    auto& fi = ws.pulay_res[i];
                    auto& fj = ws.pulay_res[j];
                    dot += fi[ig].real() * fj[ig].real()
                           + fi[ig].imag() * fj[ig].imag();
                }
                B[i * depth + j] = B[j * depth + i] = dot;
            }
        }
        // Solve B·c = 1 (with constraint Σc=1 via Lagrangian)
        // Augmented system: [B  1] [c] = [0]
        //                    [1  0] [λ]   [1]
        int n = depth + 1;
        std::vector<Real> A(n * n, Real(0));
        for (int i = 0; i < depth; ++i)
            for (int j = 0; j < depth; ++j)
                A[i * n + j] = B[i * depth + j];
        for (int i = 0; i < depth; ++i)
            A[i * n + depth] = A[depth * n + i] = Real(1);

        std::vector<Real> rhs(n, Real(0));
        rhs[depth] = Real(1);

        // Solve via Gaussian elimination with partial pivoting
        std::vector<Real> c(n, Real(0));
        auto aug = A;
        auto b = rhs;
        for (int k = 0; k < n; ++k) {
            int imax = k;
            for (int i = k + 1; i < n; ++i)
                if (std::abs(aug[i * n + k]) > std::abs(aug[imax * n + k]))
                    imax = i;
            if (imax != k) {
                for (int j = 0; j < n; ++j)
                    std::swap(aug[k * n + j], aug[imax * n + j]);
                std::swap(b[k], b[imax]);
            }
            if (std::abs(aug[k * n + k]) < Real(1e-14))
                continue;
            for (int i = k + 1; i < n; ++i) {
                Real factor = aug[i * n + k] / aug[k * n + k];
                for (int j = k; j < n; ++j)
                    aug[i * n + j] -= factor * aug[k * n + j];
                b[i] -= factor * b[k];
            }
        }
        for (int k = n - 1; k >= 0; --k) {
            Real sum = b[k];
            for (int j = k + 1; j < n; ++j)
                sum -= aug[k * n + j] * c[j];
            if (std::abs(aug[k * n + k]) > Real(1e-14))
                c[k] = sum / aug[k * n + k];
        }
        // Build mixed density: ρ_new = Σ c_i (ρ_in^(i) + α·F^(i))
        Real alpha = mix_beta_;
        for (int ig = 0; ig < ws.n_g; ++ig)
            ws.rho_G[ig] = std::complex<Real>(0, 0);
        for (int i = 0; i < depth; ++i) {
            for (int ig = 0; ig < ws.n_g; ++ig)
                ws.rho_G[ig] +=
                    c[i] * (ws.pulay_rho[i][ig] + alpha * ws.pulay_res[i][ig]);
        }
    } else {
        // Simple linear mixing for the first iteration
        for (int i = 0; i < ws.n_g; ++i)
            ws.rho_G[i] = ws.rho_G_old[i] + mix_beta_ * F[i];
    }
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
