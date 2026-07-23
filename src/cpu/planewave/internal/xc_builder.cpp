/// @file cpu/pw/internal/xc_builder.cpp
/// @brief XC potential: compute V_xc(r) from ρ(r), FFT to G-space.
#include "atomiverse/cpu/planewave/internal/xc_builder.h"
#include "atomiverse/cpu/scf/xc_functional.h"
#include <cmath>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

void XCBuilder::prepare(const PWWorkspace& ws) {
    if (ws.nr <= 0)
        return;
    rho_g_half_.resize(static_cast<std::size_t>(ws.nc));
    rho_r_.resize(static_cast<std::size_t>(ws.nr));
    vxc_g_half_.resize(static_cast<std::size_t>(ws.nc));
    if (ws.nspin == 2) {
        rho_up_r_.resize(static_cast<std::size_t>(ws.nr));
        rho_dn_r_.resize(static_cast<std::size_t>(ws.nr));
        vxc_up_r_.resize(static_cast<std::size_t>(ws.nr));
        vxc_dn_r_.resize(static_cast<std::size_t>(ws.nr));
    }
}

namespace {

/// Pad a G-sphere density into the half-complex FFT buffer and c2r → real grid.
void fft_rho_g_to_r(const PWWorkspace& ws, const xt::xarray<std::complex<AccReal>>& rho_G,
                    fftw_complex* rho_g_half, AccReal* rho_r) {
    int ng = ws.n_g;
    int nc = ws.nc;
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    int nxc = nx / 2 + 1;
    auto& gv = *ws.gvec;
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
            rho_g_half[idx][0] = static_cast<double>(rho_G[ig].real());
            rho_g_half[idx][1] = static_cast<double>(
                conj ? -rho_G[ig].imag() : rho_G[ig].imag());
        }
    }
    ws.fft->c2r(rho_g_half, rho_r);
}

/// r2c real Vxc → G-sphere coefficients.
void fft_vxc_r_to_g(const PWWorkspace& ws, AccReal* vxc_r, fftw_complex* vxc_g_half,
                    xt::xarray<std::complex<AccReal>>& V_xc_G) {
    int ng = ws.n_g;
    int nc = ws.nc;
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    int nxc = nx / 2 + 1;
    auto& gv = *ws.gvec;
    ws.fft->r2c(vxc_r, vxc_g_half);
    AccReal inv_nr = AccReal(1) / static_cast<AccReal>(ws.nr);
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
            V_xc_G[ig] = std::complex<AccReal>(
                static_cast<AccReal>(vxc_g_half[idx][0]) * inv_nr,
                static_cast<AccReal>(conj ? -vxc_g_half[idx][1]
                                          : vxc_g_half[idx][1])
                    * inv_nr);
        } else {
            V_xc_G[ig] = std::complex<AccReal>(0, 0);
        }
    }
    V_xc_G[0] = std::complex<AccReal>(V_xc_G[0].real(), AccReal(0));
}

} // namespace

void XCBuilder::step(PWWorkspace& ws) {
    int ng = ws.n_g;

    if (!ws.fft || ws.nr <= 0) {
        for (int i = 0; i < ng; ++i)
            ws.V_xc_G[i] = std::complex<AccReal>(0, 0);
        return;
    }

    int nr = ws.nr;
    int nc = ws.nc;
    int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    int nxc = nx / 2 + 1;
    auto& gv = *ws.gvec;

    // Step 1: FFT ρ(G) → ρ(r)
    // Pad ρ(G) to half-complex FFT grid, then c2r
    auto* rho_g_half = rho_g_half_.data();
    auto* rho_r = rho_r_.data();

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
            rho_g_half[idx][0] = static_cast<double>(ws.rho_G[ig].real());
            rho_g_half[idx][1] = static_cast<double>(
                conj ? -ws.rho_G[ig].imag() : ws.rho_G[ig].imag());
        }
    }

    // FFTW does NOT normalize: c2r output = Σ_G input[G] exp(iGr)
    // ρ(r) = Σ_G ρ(G) exp(iGr), so no scaling needed.
    ws.fft->c2r(rho_g_half, rho_r);

    // Step 2: Compute V_xc(r) and E_xc from ρ(r)
    if (ws.V_xc_r.size() != static_cast<std::size_t>(nr))
        ws.V_xc_r.resize({static_cast<std::size_t>(nr)});
    if (ws.exc_r.size() != static_cast<std::size_t>(nr))
        ws.exc_r.resize({static_cast<std::size_t>(nr)});

    AccReal dV = ws.omega / static_cast<AccReal>(nr); // volume element
    AccReal E_xc = 0;

    for (int i = 0; i < nr; ++i) {
        AccReal rho = rho_r[i];
        if (rho < AccReal(1e-12)) {
            ws.V_xc_r[i] = 0;
            ws.exc_r[i] = 0;
        } else {
            auto out = cpu::scf::evaluate_xc(
                cpu::scf::XCType::PBE, rho / AccReal(2), rho / AccReal(2));
            ws.V_xc_r[i] = static_cast<AccReal>(out.vrho_alpha);
            ws.exc_r[i] = static_cast<AccReal>(out.exc);
            // E_xc = ∫ ρ(r) ε_xc(r) d³r = Σ ρ_i * exc_i * dV
            E_xc += rho * static_cast<AccReal>(out.exc) * dV;
        }
    }
    ws.E_xc = E_xc;

    // Step 3: FFT V_xc(r) → V_xc(G)
    // r2c: out[G] = Σ_r in[r] exp(-iGr)  (unnormalized)
    // V_xc(G) = (1/N) Σ_r V_xc(r) exp(-iGr) = (1/N) * out[G]
    auto* vxc_g_half = vxc_g_half_.data();
    ws.fft->r2c(ws.V_xc_r.data(), vxc_g_half);

    AccReal inv_nr = AccReal(1) / static_cast<AccReal>(nr);
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
            ws.V_xc_G[ig] = std::complex<AccReal>(
                static_cast<AccReal>(vxc_g_half[idx][0]) * inv_nr,
                static_cast<AccReal>(conj ? -vxc_g_half[idx][1]
                                          : vxc_g_half[idx][1])
                    * inv_nr);
        } else {
            ws.V_xc_G[ig] = std::complex<AccReal>(0, 0);
        }
    }

    // Ensure V_xc(G=0) is real
    ws.V_xc_G[0] = std::complex<AccReal>(ws.V_xc_G[0].real(), AccReal(0));
}

void XCBuilder::step_spin(PWWorkspace& ws) {
    int ng = ws.n_g;
    if (!ws.fft || ws.nr <= 0 || ws.nspin != 2 || ws.rho_spin_G.size() != 2) {
        for (int i = 0; i < ng; ++i)
            ws.V_xc_G[i] = std::complex<AccReal>(0, 0);
        ws.E_xc = 0;
        ws.e_xc_int_spin = 0;
        return;
    }

    const int nr = ws.nr;
    auto* rho_g_half = rho_g_half_.data();
    auto* rho_up = rho_up_r_.data();
    auto* rho_dn = rho_dn_r_.data();
    auto* vxc_up = vxc_up_r_.data();
    auto* vxc_dn = vxc_dn_r_.data();

    fft_rho_g_to_r(ws, ws.rho_spin_G[0], rho_g_half, rho_up);
    fft_rho_g_to_r(ws, ws.rho_spin_G[1], rho_g_half, rho_dn);

    if (ws.V_xc_spin_G.size() != 2) {
        ws.V_xc_spin_G.resize(2);
        ws.V_xc_spin_G[0].resize({static_cast<std::size_t>(ng)});
        ws.V_xc_spin_G[1].resize({static_cast<std::size_t>(ng)});
    }
    for (int s = 0; s < 2; ++s) {
        if (ws.V_xc_spin_r[s].size() != static_cast<std::size_t>(nr))
            ws.V_xc_spin_r[s].resize({static_cast<std::size_t>(nr)});
    }
    if (ws.exc_r.size() != static_cast<std::size_t>(nr))
        ws.exc_r.resize({static_cast<std::size_t>(nr)});

    AccReal dV = ws.omega / static_cast<AccReal>(nr);
    AccReal E_xc = 0;
    AccReal e_xc_int = 0;

    for (int i = 0; i < nr; ++i) {
        AccReal ra = rho_up[i];
        AccReal rb = rho_dn[i];
        // Floor negatives from FFT ringing so LSDA stays defined.
        if (ra < AccReal(0))
            ra = AccReal(0);
        if (rb < AccReal(0))
            rb = AccReal(0);
        AccReal rtot = ra + rb;
        if (rtot < AccReal(1e-12)) {
            vxc_up[i] = 0;
            vxc_dn[i] = 0;
            ws.V_xc_spin_r[0][i] = 0;
            ws.V_xc_spin_r[1][i] = 0;
            ws.exc_r[i] = 0;
            continue;
        }
        // Same functional as restricted path (PBE with σ=0 ≡ PW92-like LSDA limit).
        auto out = cpu::scf::evaluate_xc(cpu::scf::XCType::PBE,
                                         static_cast<double>(ra),
                                         static_cast<double>(rb));
        vxc_up[i] = static_cast<AccReal>(out.vrho_alpha);
        vxc_dn[i] = static_cast<AccReal>(out.vrho_beta);
        ws.V_xc_spin_r[0][i] = vxc_up[i];
        ws.V_xc_spin_r[1][i] = vxc_dn[i];
        ws.exc_r[i] = static_cast<AccReal>(out.exc);
        E_xc += rtot * static_cast<AccReal>(out.exc) * dV;
        e_xc_int += (ra * vxc_up[i] + rb * vxc_dn[i]) * dV;
    }
    ws.E_xc = E_xc;
    ws.e_xc_int_spin = e_xc_int;

    auto* vxc_g_half = vxc_g_half_.data();
    fft_vxc_r_to_g(ws, vxc_up, vxc_g_half, ws.V_xc_spin_G[0]);
    fft_vxc_r_to_g(ws, vxc_dn, vxc_g_half, ws.V_xc_spin_G[1]);
    // Restricted-compatible total Vxc (not used by spin Davidson).
    for (int ig = 0; ig < ng; ++ig)
        ws.V_xc_G[ig] =
            AccReal(0.5) * (ws.V_xc_spin_G[0][ig] + ws.V_xc_spin_G[1][ig]);
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
