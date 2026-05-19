/// @file cpu/pw/pw_veff.cpp
/// @brief Assemble V_eff(G) = V_loc + V_H + V_xc, then FFT to real-space V_eff(r).
#include "atomiverse/cpu/pw/pw_veff.h"
#include <cstring>
#include <fftw3.h>

namespace ATV {
namespace cpu {
namespace pw {

void PW_Veff::build(const cpu::BuildContext& ctx) {
    ws_ = static_cast<PWWorkspace*>(ctx.pw_workspace);
}

void PW_Veff::launch(const cpu::LaunchContext&) {
    int ng = ws_->n_g;

    // Assemble in G-space
    for (int i = 0; i < ng; ++i)
        ws_->V_eff_G[i] = ws_->V_loc_G_real[i] + ws_->V_H_G[i] + ws_->V_xc_G[i];

    // FFT V_eff(G) → V_eff(r) for use in FFT-based H·ψ
    if (ws_->fft && ws_->nr > 0) {
        int nr = ws_->nr, nc = ws_->nc;
        int nx=ws_->fft_nx, ny=ws_->fft_ny, nz=ws_->fft_nz;
        int nxc = nx/2 + 1;
        auto& gv = *ws_->gvec;

        // Pad V_eff(G) to half-complex FFT grid
        auto* vg = static_cast<fftw_complex*>(
            fftw_malloc(sizeof(fftw_complex) * static_cast<size_t>(nc)));
        for (int i=0; i<nc; ++i) vg[i][0]=vg[i][1]=0;

        for (int ig=0; ig<ng; ++ig) {
            int h=gv.h[ig], k=gv.k[ig], l=gv.l[ig];
            int kx=((h%nx)+nx)%nx, ky=((k%ny)+ny)%ny, kz=((l%nz)+nz)%nz;
            bool conj=false;
            if (kx > nx/2) { kx=nx-kx; ky=(ny-ky)%ny; kz=(nz-kz)%nz; conj=true; }
            int idx = kx + nxc*(ky + ny*kz);
            if (idx>=0 && idx<nc) {
                vg[idx][0] = static_cast<double>(ws_->V_eff_G[ig].real());
                vg[idx][1] = static_cast<double>(conj ? -ws_->V_eff_G[ig].imag()
                                                      : ws_->V_eff_G[ig].imag());
            }
        }

        // Resize V_eff_r if needed
        if (ws_->V_eff_r.size() != static_cast<std::size_t>(nr))
            ws_->V_eff_r.resize({static_cast<std::size_t>(nr)});

        // c2r: V_eff(G) → V_eff(r). FFTW does NOT normalize.
        // V_eff(r) = Σ_G V_eff(G) exp(iGr)
        ws_->fft->c2r(vg, ws_->V_eff_r.data());

        fftw_free(vg);
    }
}

} // namespace pw
} // namespace cpu
} // namespace ATV
