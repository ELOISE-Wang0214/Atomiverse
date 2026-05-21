/// @file cpu/pw/pw_lda.cpp
/// @brief LDA XC: compute V_xc(r) from ρ(r), FFT to G-space.
#include "atomiverse/cpu/pw/pw_lda.h"
#include "atomiverse/cpu/pw/internal/fft.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include <cmath>
#include <cstring>
#include <fftw3.h>

namespace ATV {
namespace cpu {
namespace pw {

void PW_LDA::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_LDA::launch(const cpu::LaunchContext&) {
    int ng = ws_->n_g;

    if (!ws_->fft || ws_->nr <= 0) {
        for (int i = 0; i < ng; ++i)
            ws_->V_xc_G[i] = std::complex<Real>(0, 0);
        return;
    }

    int nr = ws_->nr;
    int nc = ws_->nc;
    int nx = ws_->fft_nx, ny = ws_->fft_ny, nz = ws_->fft_nz;
    int nxc = nx / 2 + 1;
    auto& gv = *ws_->gvec;

    // Step 1: FFT ρ(G) → ρ(r)
    // Pad ρ(G) to half-complex FFT grid, then c2r
    auto* rho_g_half = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nc)));
    auto* rho_r =
        static_cast<Real*>(fftw_malloc(sizeof(Real) * static_cast<size_t>(nr)));

    for (int i = 0; i < nc; ++i)
        rho_g_half[i][0] = rho_g_half[i][1] = 0;

    for (int ig = 0; ig < ng; ++ig) {
        int h = gv.h[ig], k = gv.k[ig], l = gv.l[ig];
        int kx = ((h % nx) + nx) % nx;
        int ky = ((k % ny) + ny) % ny;
        int kz = ((l % nz) + nz) % nz;
        bool conj = false;
        if (kx > nx / 2) {
            kx = nx - kx;
            ky = (ny - ky) % ny;
            kz = (nz - kz) % nz;
            conj = true;
        }
        int idx = kx + nxc * (ky + ny * kz);
        if (idx >= 0 && idx < nc) {
            rho_g_half[idx][0] = static_cast<double>(ws_->rho_G[ig].real());
            rho_g_half[idx][1] = static_cast<double>(
                conj ? -ws_->rho_G[ig].imag() : ws_->rho_G[ig].imag());
        }
    }

    // FFTW does NOT normalize: c2r output = Σ_G input[G] exp(iGr)
    // ρ(r) = Σ_G ρ(G) exp(iGr), so no scaling needed.
    ws_->fft->c2r(rho_g_half, rho_r);

    // Step 2: Compute V_xc(r) and E_xc from ρ(r)
    if (ws_->V_xc_r.size() != static_cast<std::size_t>(nr))
        ws_->V_xc_r.resize({static_cast<std::size_t>(nr)});
    if (ws_->exc_r.size() != static_cast<std::size_t>(nr))
        ws_->exc_r.resize({static_cast<std::size_t>(nr)});

    Real dV = ws_->omega / static_cast<Real>(nr); // volume element
    Real E_xc = 0;

    for (int i = 0; i < nr; ++i) {
        Real rho = rho_r[i];
        if (rho < Real(1e-12)) {
            ws_->V_xc_r[i] = 0;
            ws_->exc_r[i] = 0;
        } else {
            auto out = cpu::scf::evaluate_xc(cpu::scf::XCType::PBE,
                                             rho / Real(2), rho / Real(2));
            ws_->V_xc_r[i] = static_cast<Real>(out.vrho_alpha);
            ws_->exc_r[i] = static_cast<Real>(out.exc);
            // E_xc = ∫ ρ(r) ε_xc(r) d³r = Σ ρ_i * exc_i * dV
            E_xc += rho * static_cast<Real>(out.exc) * dV;
        }
    }
    ws_->E_xc = E_xc;

    // Step 3: FFT V_xc(r) → V_xc(G)
    // r2c: out[G] = Σ_r in[r] exp(-iGr)  (unnormalized)
    // V_xc(G) = (1/N) Σ_r V_xc(r) exp(-iGr) = (1/N) * out[G]
    auto* vxc_g_half = static_cast<fftw_complex*>(
        fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nc)));
    ws_->fft->r2c(ws_->V_xc_r.data(), vxc_g_half);

    Real inv_nr = Real(1) / static_cast<Real>(nr);
    for (int ig = 0; ig < ng; ++ig) {
        int h = gv.h[ig], k = gv.k[ig], l = gv.l[ig];
        int kx = ((h % nx) + nx) % nx;
        int ky = ((k % ny) + ny) % ny;
        int kz = ((l % nz) + nz) % nz;
        bool conj = false;
        if (kx > nx / 2) {
            kx = nx - kx;
            ky = (ny - ky) % ny;
            kz = (nz - kz) % nz;
            conj = true;
        }
        int idx = kx + nxc * (ky + ny * kz);
        if (idx >= 0 && idx < nc) {
            ws_->V_xc_G[ig] = std::complex<Real>(
                static_cast<Real>(vxc_g_half[idx][0]) * inv_nr,
                static_cast<Real>(conj ? -vxc_g_half[idx][1]
                                       : vxc_g_half[idx][1])
                    * inv_nr);
        } else {
            ws_->V_xc_G[ig] = std::complex<Real>(0, 0);
        }
    }

    // Ensure V_xc(G=0) is real
    ws_->V_xc_G[0] = std::complex<Real>(ws_->V_xc_G[0].real(), Real(0));

    fftw_free(rho_g_half);
    fftw_free(rho_r);
    fftw_free(vxc_g_half);
}

} // namespace pw
} // namespace cpu
} // namespace ATV
