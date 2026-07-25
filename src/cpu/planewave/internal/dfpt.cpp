/// @file cpu/pw/internal/dfpt.cpp
/// @brief Γ-point DFPT force constants for the plane-wave method.
///
/// Implements the q = 0 dynamical-matrix assembly of
/// Baroni–de Gironcoli–Dal Corso–Giannozzi, Rev. Mod. Phys. 73, 515 (2001):
///
///   C[Iα][Jβ] = C_ion (analytic Ewald Hessian)
///             + Σ_{k,i} w f ⟨ψ_i|∂²V_ext^{(I)}/∂R_α∂R_β|ψ_i⟩ δ_IJ
///             + Σ_{k,i} w f · 2 Re⟨Δψ_i^{(Jβ)}|∂V_ext^{(Iα)}|ψ_i⟩,
///
/// with the self-consistent Sternheimer response
/// `P_c (H − ε_i) P_c |Δψ⟩ = −P_c (ΔV_bare + ΔV_H[Δn] + f_xc Δn)|ψ_i⟩`
/// solved by conduction-space projected preconditioned CG.
///
/// Grid/sphere conventions deliberately mirror the ground-state SCF
/// (`density_builder` / `veff_builder` / `apply_hamiltonian`): every
/// multiplicative potential is restricted to the Γ density G-sphere before
/// being applied on the FFT grid, so the DFPT Hessian is the exact second
/// derivative of the *discretized* total energy the analytic forces
/// differentiate once (FD-vs-DFPT gate in `reg_dfpt-phonon-si`).
///
/// The KB projector loops intentionally mirror `davidson.cpp` /
/// `pw_forces.cpp` (per-G Y_lm evaluation on q = k+G); folding all three
/// into one shared skeleton is a deliberate follow-up, not a drive-by.

#include "atomiverse/cpu/planewave/internal/dfpt.h"

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/log.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/planewave/internal/davidson.h"
#include "atomiverse/cpu/planewave/internal/density_builder.h"
#include "atomiverse/cpu/planewave/internal/poisson.h"
#include "atomiverse/cpu/planewave/internal/pw_setup.h"
#include "atomiverse/cpu/planewave/internal/pw_subspace.h"
#include "atomiverse/cpu/planewave/internal/scf_loop.h"
#include "atomiverse/cpu/planewave/internal/upf_reader.h"
#include "atomiverse/cpu/planewave/internal/veff_builder.h"
#include "atomiverse/cpu/planewave/internal/xc_builder.h"
#include "atomiverse/cpu/scf/xc_functional.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace ATV {
namespace cpu {
namespace pw {
namespace internal {

namespace {

using Cx = std::complex<AccReal>;

// ── G-sphere ↔ full c2c FFT-grid staging (density_builder conventions) ──

void pad_sphere_to_grid(const PWWorkspace& ws, const GVectorList& gv, int ng,
                        const Cx* f, fftw_complex* grid) {
    const int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz, nr = ws.nr;
    for (int i = 0; i < nr; ++i)
        grid[i][0] = grid[i][1] = 0;
    for (int ig = 0; ig < ng; ++ig) {
        const int ix = ((gv.h[ig] % nx) + nx) % nx;
        const int iy = ((gv.k[ig] % ny) + ny) % ny;
        const int iz = ((gv.l[ig] % nz) + nz) % nz;
        const int idx = ix + nx * (iy + ny * iz);
        if (idx >= 0 && idx < nr) {
            grid[idx][0] = static_cast<double>(f[ig].real());
            grid[idx][1] = static_cast<double>(f[ig].imag());
        }
    }
}

void extract_sphere_from_grid(const PWWorkspace& ws, const GVectorList& gv,
                              int ng, const fftw_complex* grid, Cx* f) {
    const int nx = ws.fft_nx, ny = ws.fft_ny, nz = ws.fft_nz;
    for (int ig = 0; ig < ng; ++ig) {
        const int ix = ((gv.h[ig] % nx) + nx) % nx;
        const int iy = ((gv.k[ig] % ny) + ny) % ny;
        const int iz = ((gv.l[ig] % nz) + nz) % nz;
        const int idx = ix + nx * (iy + ny * iz);
        f[ig] = Cx(static_cast<AccReal>(grid[idx][0]),
                   static_cast<AccReal>(grid[idx][1]));
    }
}

/// Real-space field V(r) = Σ_G V(G) e^{iGr} from Γ-sphere coefficients
/// (Hermitian input ⇒ real output; the residual imaginary part is rounding).
void sphere_to_real(const PWWorkspace& ws, const GVectorList& gv0, const Cx* vg,
                    AccReal* vr, fftw_complex* work) {
    pad_sphere_to_grid(ws, gv0, gv0.count(), vg, work);
    ws.fft->c2c_bwd(work);
    for (int i = 0; i < ws.nr; ++i)
        vr[i] = static_cast<AccReal>(work[i][0]);
}

/// Γ-sphere coefficients f(G) = (1/N) Σ_r f(r) e^{−iGr} of a real field.
void real_to_sphere(const PWWorkspace& ws, const GVectorList& gv0,
                    const AccReal* fr, Cx* fg, fftw_complex* work) {
    const int nr = ws.nr;
    for (int i = 0; i < nr; ++i) {
        work[i][0] = static_cast<double>(fr[i]);
        work[i][1] = 0.0;
    }
    ws.fft->c2c_fwd(work);
    std::vector<Cx> full(static_cast<std::size_t>(gv0.count()));
    extract_sphere_from_grid(ws, gv0, gv0.count(), work, full.data());
    const AccReal inv_nr = AccReal(1) / AccReal(nr);
    for (int ig = 0; ig < gv0.count(); ++ig)
        fg[ig] = full[static_cast<std::size_t>(ig)] * inv_nr;
}

/// (V·ψ)(G) on the bound k-sphere from a real-space multiplicative V(r)
/// (apply_hamiltonian's local step: bwd FFT, multiply, fwd FFT, 1/N).
void apply_mult_potential(const PWWorkspace& ws, const AccReal* vr,
                          const Cx* psi, Cx* out, fftw_complex* work) {
    const auto& gv = *ws.gvec;
    pad_sphere_to_grid(ws, gv, ws.n_g, psi, work);
    ws.fft->c2c_bwd(work);
    for (int i = 0; i < ws.nr; ++i) {
        const double v = static_cast<double>(vr[i]);
        work[i][0] *= v;
        work[i][1] *= v;
    }
    ws.fft->c2c_fwd(work);
    extract_sphere_from_grid(ws, gv, ws.n_g, work, out);
    const AccReal inv_nr = AccReal(1) / AccReal(ws.nr);
    for (int ig = 0; ig < ws.n_g; ++ig)
        out[ig] *= inv_nr;
}

std::complex<AccReal> im_l_factor(int l) {
    if (l == 1)
        return {0, -1};
    if (l == 2)
        return {-1, 0};
    return {1, 0};
}

/// Per-(band, atom, projector, m) KB contraction bundle on the bound
/// k-sphere: c = ⟨β|ψ⟩, s_α = Σ q_α β*ψ (∂c/∂R_α = i s_α) and the six
/// second-moment sums t_αβ = Σ q_α q_β β*ψ (∂²c/∂R_α∂R_β = −t_αβ).
struct KbMoments {
    Cx c;
    Cx s[3];
    Cx t[3][3];
};

/// Evaluate KB moments of `psi` for every (m) of one projector of one atom.
/// `want_t` skips the second moments when only c/s are needed.
void kb_moments(const PWWorkspace& ws, std::size_t ia, std::size_t ip,
                const Cx* psi, bool want_t, std::vector<KbMoments>& out) {
    const auto& gv = *ws.gvec;
    const auto& proj = ws.kb_proj[ia][ip];
    const int l = proj.l;
    const int nm = n_m_for_l(l);
    const AccReal pref = AccReal(4.0 * M_PI) / std::sqrt(ws.omega);
    const auto im_l = im_l_factor(l);
    const AccReal ax = ws.atom_pos_bohr[ia][0];
    const AccReal ay = ws.atom_pos_bohr[ia][1];
    const AccReal az = ws.atom_pos_bohr[ia][2];

    out.assign(static_cast<std::size_t>(nm), KbMoments{});
    std::vector<AccReal> ylm(static_cast<std::size_t>(nm));

    for (int ig = 0; ig < ws.n_g; ++ig) {
        const AccReal fl = proj.f_G[ig];
        if (std::abs(fl) < AccReal(1e-30))
            continue;
        const AccReal qx = gv.g[ig][0] + gv.kvec[0];
        const AccReal qy = gv.g[ig][1] + gv.kvec[1];
        const AccReal qz = gv.g[ig][2] + gv.kvec[2];
        const AccReal qr = std::sqrt(qx * qx + qy * qy + qz * qz);
        const AccReal phase = qx * ax + qy * ay + qz * az;
        const AccReal cos_p = std::cos(phase), sin_p = std::sin(phase);
        compute_ylm(qx, qy, qz, qr, l, ylm.data());
        const AccReal q[3] = {qx, qy, qz};

        for (int m = 0; m < nm; ++m) {
            const Cx beta = pref * im_l * ylm[static_cast<std::size_t>(m)] * fl
                            * Cx(cos_p, -sin_p);
            const Cx bp = std::conj(beta) * psi[ig];
            auto& mo = out[static_cast<std::size_t>(m)];
            mo.c += bp;
            for (int a = 0; a < 3; ++a)
                mo.s[a] += q[a] * bp;
            if (want_t)
                for (int a = 0; a < 3; ++a)
                    for (int b = a; b < 3; ++b)
                        mo.t[a][b] += q[a] * q[b] * bp;
        }
    }
    if (want_t)
        for (auto& mo : out)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < a; ++b)
                    mo.t[a][b] = mo.t[b][a];
}

/// Accumulate (∂V_nl^{(I,α)} ψ)(G) for atom `ia`, axis `alpha` into `out`:
/// Σ_{p,m} i D β(G) (s_α − q_α c).
void add_dnl_psi(const PWWorkspace& ws, std::size_t ia, int alpha,
                 const Cx* psi, Cx* out) {
    if (ws.kb_proj.empty() || ia >= ws.kb_proj.size())
        return;
    const auto& gv = *ws.gvec;
    const AccReal pref = AccReal(4.0 * M_PI) / std::sqrt(ws.omega);
    const AccReal ax = ws.atom_pos_bohr[ia][0];
    const AccReal ay = ws.atom_pos_bohr[ia][1];
    const AccReal az = ws.atom_pos_bohr[ia][2];

    std::vector<KbMoments> mo;
    for (std::size_t ip = 0; ip < ws.kb_proj[ia].size(); ++ip) {
        const auto& proj = ws.kb_proj[ia][ip];
        const AccReal D = proj.d;
        if (std::abs(D) < AccReal(1e-16))
            continue;
        const int l = proj.l;
        const int nm = n_m_for_l(l);
        const auto im_l = im_l_factor(l);
        kb_moments(ws, ia, ip, psi, /*want_t=*/false, mo);

        std::vector<AccReal> ylm(static_cast<std::size_t>(nm));
        for (int ig = 0; ig < ws.n_g; ++ig) {
            const AccReal fl = proj.f_G[ig];
            if (std::abs(fl) < AccReal(1e-30))
                continue;
            const AccReal qx = gv.g[ig][0] + gv.kvec[0];
            const AccReal qy = gv.g[ig][1] + gv.kvec[1];
            const AccReal qz = gv.g[ig][2] + gv.kvec[2];
            const AccReal qr = std::sqrt(qx * qx + qy * qy + qz * qz);
            const AccReal q[3] = {qx, qy, qz};
            const AccReal phase = qx * ax + qy * ay + qz * az;
            const AccReal cos_p = std::cos(phase), sin_p = std::sin(phase);
            compute_ylm(qx, qy, qz, qr, l, ylm.data());
            for (int m = 0; m < nm; ++m) {
                const Cx beta = pref * im_l * ylm[static_cast<std::size_t>(m)]
                                * fl * Cx(cos_p, -sin_p);
                const auto& mm = mo[static_cast<std::size_t>(m)];
                out[ig] +=
                    Cx(0, 1) * D * beta * (mm.s[alpha] - q[alpha] * mm.c);
            }
        }
    }
}

AccReal lda_fxc_at(AccReal rho, cpu::scf::XCType xc) {
    if (rho < AccReal(1e-10))
        return AccReal(0);
    const double r = static_cast<double>(rho);
    const double h = std::max(1e-6 * r, 1e-12);
    const double vp =
        cpu::scf::evaluate_xc(xc, (r + h) / 2.0, (r + h) / 2.0).vrho_alpha;
    const double vm =
        cpu::scf::evaluate_xc(xc, (r - h) / 2.0, (r - h) / 2.0).vrho_alpha;
    return static_cast<AccReal>((vp - vm) / (2.0 * h));
}

/// Project out the occupied manifold: x ← x − Σ_j ψ_j ⟨ψ_j|x⟩.
void project_conduction(const Cx* occ_wfc, int n_occ, int ng, Cx* x) {
    for (int j = 0; j < n_occ; ++j) {
        const Cx* psi = occ_wfc + static_cast<std::ptrdiff_t>(j) * ng;
        Cx ovlp(0, 0);
        for (int ig = 0; ig < ng; ++ig)
            ovlp += std::conj(psi[ig]) * x[ig];
        for (int ig = 0; ig < ng; ++ig)
            x[ig] -= ovlp * psi[ig];
    }
}

} // namespace

// ── Ewald force constants ────────────────────────────────────────────

std::vector<AccReal> ewald_dynamical_matrix(const std::vector<Vec3>& pos_ang,
                                            const std::vector<AccReal>& charges,
                                            const cpu::Box& cell,
                                            AccReal ecut_hartree) {
    const AccReal a2b = AccReal(1.8897259886);
    cpu::Box cell_bohr;
    cell_bohr.set_cell(cell.cell() * a2b);
    const Mat3& cb = cell_bohr.cell();

    const int n = static_cast<int>(charges.size());
    const int dim = 3 * n;
    std::vector<AccReal> D(static_cast<std::size_t>(dim) * dim, AccReal(0));
    if (n == 0)
        return D;

    std::vector<Vec3> pos(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        pos[static_cast<std::size_t>(i)] = {
            pos_ang[static_cast<std::size_t>(i)][0] * a2b,
            pos_ang[static_cast<std::size_t>(i)][1] * a2b,
            pos_ang[static_cast<std::size_t>(i)][2] * a2b};

    const AccReal vol = cell_bohr.volume();
    if (vol <= 0)
        return D;

    // α selection and cutoffs identical to EwaldSummation.
    const AccReal rcut_alpha = cell_bohr.lengths().norm() * AccReal(0.5);
    const AccReal alpha = std::sqrt(-std::log(AccReal(2e-17))) / rcut_alpha;
    const AccReal two_alpha_over_sqrt_pi =
        AccReal(2) * alpha / std::sqrt(AccReal(units::pi));

    auto at = [&](int i, int a, int j, int b) -> AccReal& {
        return D[static_cast<std::size_t>((3 * i + a) * dim) + 3 * j + b];
    };

    // ── Real space ──
    const AccReal r_cut = AccReal(5) / alpha;
    const auto lens = cell_bohr.lengths();
    const int ncell_r = static_cast<int>(std::ceil(
                            r_cut / std::min({lens[0], lens[1], lens[2]})))
                        + 1;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const AccReal zz = charges[static_cast<std::size_t>(i)]
                               * charges[static_cast<std::size_t>(j)];
            for (int hx = -ncell_r; hx <= ncell_r; ++hx)
                for (int hy = -ncell_r; hy <= ncell_r; ++hy)
                    for (int hz = -ncell_r; hz <= ncell_r; ++hz) {
                        // L = hx a + hy b + hz c (rows of cb = lattice vectors).
                        const AccReal Lx = AccReal(hx) * cb(0, 0)
                                           + AccReal(hy) * cb(1, 0)
                                           + AccReal(hz) * cb(2, 0);
                        const AccReal Ly = AccReal(hx) * cb(0, 1)
                                           + AccReal(hy) * cb(1, 1)
                                           + AccReal(hz) * cb(2, 1);
                        const AccReal Lz = AccReal(hx) * cb(0, 2)
                                           + AccReal(hy) * cb(1, 2)
                                           + AccReal(hz) * cb(2, 2);
                        const Vec3 dr = {
                            pos[static_cast<std::size_t>(i)][0]
                                - pos[static_cast<std::size_t>(j)][0] + Lx,
                            pos[static_cast<std::size_t>(i)][1]
                                - pos[static_cast<std::size_t>(j)][1] + Ly,
                            pos[static_cast<std::size_t>(i)][2]
                                - pos[static_cast<std::size_t>(j)][2] + Lz};
                        const AccReal r = dr.norm();
                        if (!(r > AccReal(0)) || r >= r_cut)
                            continue;
                        const AccReal inv_r = AccReal(1) / r;
                        const AccReal inv_r2 = inv_r * inv_r;
                        const AccReal erfc_ar = std::erfc(alpha * r);
                        const AccReal gauss = std::exp(-alpha * alpha * r * r);
                        // φ = erfc(αr)/r; φ′ and φ″ analytic.
                        const AccReal phi1 =
                            -erfc_ar * inv_r2
                            - two_alpha_over_sqrt_pi * gauss * inv_r;
                        const AccReal phi2 =
                            AccReal(2) * erfc_ar * inv_r2 * inv_r
                            + two_alpha_over_sqrt_pi * gauss
                                  * (AccReal(2) * inv_r2
                                     + AccReal(2) * alpha * alpha);
                        for (int a = 0; a < 3; ++a)
                            for (int b = 0; b < 3; ++b) {
                                const AccReal ra = dr[a] * inv_r;
                                const AccReal rb = dr[b] * inv_r;
                                const AccReal hab =
                                    phi2 * ra * rb
                                    + phi1 * inv_r
                                          * ((a == b ? AccReal(1) : AccReal(0))
                                             - ra * rb);
                                const AccReal val = zz * hab;
                                at(i, a, i, b) += val;
                                at(j, a, j, b) += val;
                                at(i, a, j, b) -= val;
                                at(j, a, i, b) -= val;
                            }
                    }
        }
    }

    // ── Reciprocal space ──
    const AccReal g_cut = std::sqrt(AccReal(8) * ecut_hartree);
    const auto recip = cell_bohr.recip_cell();
    const Vec3 b0{recip(0, 0), recip(0, 1), recip(0, 2)};
    const Vec3 b1{recip(1, 0), recip(1, 1), recip(1, 2)};
    const Vec3 b2{recip(2, 0), recip(2, 1), recip(2, 2)};
    const int gmax =
        static_cast<int>(g_cut / std::min({b0.norm(), b1.norm(), b2.norm()}))
        + 1;
    const AccReal alpha_sq_4 = AccReal(4) * alpha * alpha;
    const AccReal two_pi_over_vol = AccReal(2) * AccReal(units::pi) / vol;

    for (int hx = -gmax; hx <= gmax; ++hx)
        for (int hy = -gmax; hy <= gmax; ++hy)
            for (int hz = -gmax; hz <= gmax; ++hz) {
                if (hx == 0 && hy == 0 && hz == 0)
                    continue;
                const Vec3 g =
                    b0 * AccReal(hx) + b1 * AccReal(hy) + b2 * AccReal(hz);
                const AccReal g2 = g.norm_squared();
                if (g2 > g_cut * g_cut)
                    continue;
                const AccReal coeff =
                    two_pi_over_vol * std::exp(-g2 / alpha_sq_4) / g2;
                for (int i = 0; i < n; ++i)
                    for (int j = i + 1; j < n; ++j) {
                        const AccReal phase =
                            g[0]
                                * (pos[static_cast<std::size_t>(i)][0]
                                   - pos[static_cast<std::size_t>(j)][0])
                            + g[1]
                                  * (pos[static_cast<std::size_t>(i)][1]
                                     - pos[static_cast<std::size_t>(j)][1])
                            + g[2]
                                  * (pos[static_cast<std::size_t>(i)][2]
                                     - pos[static_cast<std::size_t>(j)][2]);
                        const AccReal zz =
                            charges[static_cast<std::size_t>(i)]
                            * charges[static_cast<std::size_t>(j)];
                        const AccReal fac =
                            coeff * AccReal(2) * zz * std::cos(phase);
                        for (int a = 0; a < 3; ++a)
                            for (int b = 0; b < 3; ++b) {
                                const AccReal val = fac * g[a] * g[b];
                                at(i, a, j, b) += val;
                                at(j, a, i, b) += val;
                                at(i, a, i, b) -= val;
                                at(j, a, j, b) -= val;
                            }
                    }
            }

    return D;
}

std::vector<std::complex<AccReal>>
ewald_dynamical_matrix_q(const std::vector<Vec3>& pos_ang,
                         const std::vector<AccReal>& charges,
                         const cpu::Box& cell, AccReal ecut_hartree,
                         const Vec3& q_frac) {
    using Cx = std::complex<AccReal>;
    const AccReal a2b = AccReal(1.8897259886);
    cpu::Box cell_bohr;
    cell_bohr.set_cell(cell.cell() * a2b);
    const Mat3& cb = cell_bohr.cell();

    const int n = static_cast<int>(charges.size());
    const int dim = 3 * n;
    std::vector<Cx> D(static_cast<std::size_t>(dim) * dim, Cx(0, 0));
    if (n == 0)
        return D;

    // Exact Γ matrix when q vanishes (identical α / cutoffs / accumulation).
    if (std::abs(q_frac[0]) < AccReal(1e-14) && std::abs(q_frac[1]) < AccReal(1e-14)
        && std::abs(q_frac[2]) < AccReal(1e-14)) {
        const auto Dg =
            ewald_dynamical_matrix(pos_ang, charges, cell, ecut_hartree);
        for (std::size_t i = 0; i < Dg.size(); ++i)
            D[i] = Dg[i];
        return D;
    }

    std::vector<Vec3> pos(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        pos[static_cast<std::size_t>(i)] = {
            pos_ang[static_cast<std::size_t>(i)][0] * a2b,
            pos_ang[static_cast<std::size_t>(i)][1] * a2b,
            pos_ang[static_cast<std::size_t>(i)][2] * a2b};

    const AccReal vol = cell_bohr.volume();
    if (vol <= 0)
        return D;

    // Cartesian q [1/Bohr] from fractional crystal coords (recip includes 2π).
    const auto recip = cell_bohr.recip_cell();
    const Vec3 b0{recip(0, 0), recip(0, 1), recip(0, 2)};
    const Vec3 b1{recip(1, 0), recip(1, 1), recip(1, 2)};
    const Vec3 b2{recip(2, 0), recip(2, 1), recip(2, 2)};
    const Vec3 q_cart = b0 * q_frac[0] + b1 * q_frac[1] + b2 * q_frac[2];
    const AccReal rcut_alpha = cell_bohr.lengths().norm() * AccReal(0.5);
    const AccReal alpha = std::sqrt(-std::log(AccReal(2e-17))) / rcut_alpha;
    const AccReal two_alpha_over_sqrt_pi =
        AccReal(2) * alpha / std::sqrt(AccReal(units::pi));

    auto at = [&](int i, int a, int j, int b) -> Cx& {
        return D[static_cast<std::size_t>((3 * i + a) * dim) + 3 * j + b];
    };

    // ── Real space: Φ(i, j; L) e^{i q · L} ──
    const AccReal r_cut = AccReal(5) / alpha;
    const auto lens = cell_bohr.lengths();
    const int ncell_r = static_cast<int>(std::ceil(
                            r_cut / std::min({lens[0], lens[1], lens[2]})))
                        + 1;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const AccReal zz = charges[static_cast<std::size_t>(i)]
                               * charges[static_cast<std::size_t>(j)];
            for (int hx = -ncell_r; hx <= ncell_r; ++hx)
                for (int hy = -ncell_r; hy <= ncell_r; ++hy)
                    for (int hz = -ncell_r; hz <= ncell_r; ++hz) {
                        if (i == j && hx == 0 && hy == 0 && hz == 0)
                            continue;
                        const Vec3 L = {
                            AccReal(hx) * cb(0, 0) + AccReal(hy) * cb(1, 0)
                                + AccReal(hz) * cb(2, 0),
                            AccReal(hx) * cb(0, 1) + AccReal(hy) * cb(1, 1)
                                + AccReal(hz) * cb(2, 1),
                            AccReal(hx) * cb(0, 2) + AccReal(hy) * cb(1, 2)
                                + AccReal(hz) * cb(2, 2)};
                        const Vec3 dr = {
                            pos[static_cast<std::size_t>(i)][0]
                                - pos[static_cast<std::size_t>(j)][0] - L[0],
                            pos[static_cast<std::size_t>(i)][1]
                                - pos[static_cast<std::size_t>(j)][1] - L[1],
                            pos[static_cast<std::size_t>(i)][2]
                                - pos[static_cast<std::size_t>(j)][2] - L[2]};
                        const AccReal r = dr.norm();
                        if (!(r > AccReal(0)) || r >= r_cut)
                            continue;
                        const AccReal inv_r = AccReal(1) / r;
                        const AccReal inv_r2 = inv_r * inv_r;
                        const AccReal erfc_ar = std::erfc(alpha * r);
                        const AccReal gauss = std::exp(-alpha * alpha * r * r);
                        const AccReal phi1 =
                            -erfc_ar * inv_r2
                            - two_alpha_over_sqrt_pi * gauss * inv_r;
                        const AccReal phi2 =
                            AccReal(2) * erfc_ar * inv_r2 * inv_r
                            + two_alpha_over_sqrt_pi * gauss
                                  * (AccReal(2) * inv_r2
                                     + AccReal(2) * alpha * alpha);
                        // phase = e^{i q · L}; L in Bohr, q_cart in 1/Bohr.
                        // recip_cell rows are 2π b*, so q_cart already has 2π.
                        const AccReal ph =
                            q_cart[0] * L[0] + q_cart[1] * L[1]
                            + q_cart[2] * L[2];
                        const Cx phase(std::cos(ph), std::sin(ph));
                        for (int a = 0; a < 3; ++a)
                            for (int b = 0; b < 3; ++b) {
                                const AccReal ra = dr[a] * inv_r;
                                const AccReal rb = dr[b] * inv_r;
                                const AccReal hab =
                                    phi2 * ra * rb
                                    + phi1 * inv_r
                                          * ((a == b ? AccReal(1) : AccReal(0))
                                             - ra * rb);
                                // Φ_{iα,jβ}(L) = −∂²φ/∂r_α∂r_β for the pair
                                // (force-constant convention matching Γ path
                                // off-diagonals); self terms via ASR below.
                                const Cx val = zz * hab * phase;
                                at(i, a, j, b) -= val;
                            }
                    }
        }
    }

    // Acoustic sum rule at the real-space level: diagonal self blocks so
    // that sum_j C_{iα,jβ}(q=0) = 0. Store row sums at q via separate
    // accumulation — for finite q, Gonze ASR uses the q=0 row sum on the
    // diagonal. Compute ASR from the q=0 real-space pair sum only.
    {
        // Recompute q=0 real-space diagonal correction: for each i,a,b
        // C_{ia,ib} += -sum_{j,L≠0 or j≠i} Φ_{ia,jb}(L) with phase=1.
        // Equivalent: after the loop above at general q, add
        //   C_{ia,ib} += -sum_j C^{RS}_{ia,jb}(q=0)
        // We rebuild the q=0 real-space matrix contribution on the diagonal.
        std::vector<AccReal> asr(static_cast<std::size_t>(n) * 9, AccReal(0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const AccReal zz = charges[static_cast<std::size_t>(i)]
                                   * charges[static_cast<std::size_t>(j)];
                for (int hx = -ncell_r; hx <= ncell_r; ++hx)
                    for (int hy = -ncell_r; hy <= ncell_r; ++hy)
                        for (int hz = -ncell_r; hz <= ncell_r; ++hz) {
                            if (i == j && hx == 0 && hy == 0 && hz == 0)
                                continue;
                            const Vec3 L = {
                                AccReal(hx) * cb(0, 0) + AccReal(hy) * cb(1, 0)
                                    + AccReal(hz) * cb(2, 0),
                                AccReal(hx) * cb(0, 1) + AccReal(hy) * cb(1, 1)
                                    + AccReal(hz) * cb(2, 1),
                                AccReal(hx) * cb(0, 2) + AccReal(hy) * cb(1, 2)
                                    + AccReal(hz) * cb(2, 2)};
                            const Vec3 dr = {
                                pos[static_cast<std::size_t>(i)][0]
                                    - pos[static_cast<std::size_t>(j)][0]
                                    - L[0],
                                pos[static_cast<std::size_t>(i)][1]
                                    - pos[static_cast<std::size_t>(j)][1]
                                    - L[1],
                                pos[static_cast<std::size_t>(i)][2]
                                    - pos[static_cast<std::size_t>(j)][2]
                                    - L[2]};
                            const AccReal r = dr.norm();
                            if (!(r > AccReal(0)) || r >= r_cut)
                                continue;
                            const AccReal inv_r = AccReal(1) / r;
                            const AccReal inv_r2 = inv_r * inv_r;
                            const AccReal erfc_ar = std::erfc(alpha * r);
                            const AccReal gauss =
                                std::exp(-alpha * alpha * r * r);
                            const AccReal phi1 =
                                -erfc_ar * inv_r2
                                - two_alpha_over_sqrt_pi * gauss * inv_r;
                            const AccReal phi2 =
                                AccReal(2) * erfc_ar * inv_r2 * inv_r
                                + two_alpha_over_sqrt_pi * gauss
                                      * (AccReal(2) * inv_r2
                                         + AccReal(2) * alpha * alpha);
                            for (int a = 0; a < 3; ++a)
                                for (int b = 0; b < 3; ++b) {
                                    const AccReal ra = dr[a] * inv_r;
                                    const AccReal rb = dr[b] * inv_r;
                                    const AccReal hab =
                                        phi2 * ra * rb
                                        + phi1 * inv_r
                                              * ((a == b ? AccReal(1)
                                                         : AccReal(0))
                                                 - ra * rb);
                                    // Φ_ij = -zz*hab; ASR: C_ii += -sum_j Φ_ij
                                    //      = sum_j zz*hab
                                    asr[static_cast<std::size_t>(i) * 9 + 3 * a
                                        + b] += zz * hab;
                                }
                        }
            }
        }
        for (int i = 0; i < n; ++i)
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    at(i, a, i, b) +=
                        asr[static_cast<std::size_t>(i) * 9 + 3 * a + b];
    }

    // ── Reciprocal space ──
    // (4π/Ω) Σ_{G+q ≠ 0} e^{-|G+q|²/4α²}/|G+q|² (G+q)_α (G+q)_β
    //   × Z_i Z_j exp(i (G+q)·(R_i − R_j))
    // plus the G-sum self terms that restore ASR at the reciprocal level.
    const AccReal g_cut = std::sqrt(AccReal(8) * ecut_hartree);
    // Expand search to cover |G+q|.
    const AccReal qn = q_cart.norm();
    const AccReal g_cut_eff = g_cut + qn + AccReal(1e-6);
    const int gmax =
        static_cast<int>(g_cut_eff
                         / std::min({b0.norm(), b1.norm(), b2.norm()}))
        + 1;
    const AccReal alpha_sq_4 = AccReal(4) * alpha * alpha;
    const AccReal four_pi_over_vol =
        AccReal(4) * AccReal(units::pi) / vol;

    // Structure-factor form: for each G, accumulate
    //   S_α(i) pieces via double loop over i,j.
    for (int hx = -gmax; hx <= gmax; ++hx)
        for (int hy = -gmax; hy <= gmax; ++hy)
            for (int hz = -gmax; hz <= gmax; ++hz) {
                const Vec3 G =
                    b0 * AccReal(hx) + b1 * AccReal(hy) + b2 * AccReal(hz);
                const Vec3 Gq = {G[0] + q_cart[0], G[1] + q_cart[1],
                                 G[2] + q_cart[2]};
                const AccReal gq2 = Gq.norm_squared();
                if (gq2 < AccReal(1e-16) || gq2 > g_cut_eff * g_cut_eff)
                    continue;
                const AccReal coeff =
                    four_pi_over_vol * std::exp(-gq2 / alpha_sq_4) / gq2;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j) {
                        const AccReal phase_arg =
                            Gq[0]
                                * (pos[static_cast<std::size_t>(i)][0]
                                   - pos[static_cast<std::size_t>(j)][0])
                            + Gq[1]
                                  * (pos[static_cast<std::size_t>(i)][1]
                                     - pos[static_cast<std::size_t>(j)][1])
                            + Gq[2]
                                  * (pos[static_cast<std::size_t>(i)][2]
                                     - pos[static_cast<std::size_t>(j)][2]);
                        const Cx phase(std::cos(phase_arg),
                                       std::sin(phase_arg));
                        const AccReal zz =
                            charges[static_cast<std::size_t>(i)]
                            * charges[static_cast<std::size_t>(j)];
                        const Cx fac = coeff * zz * phase;
                        for (int a = 0; a < 3; ++a)
                            for (int b = 0; b < 3; ++b)
                                at(i, a, j, b) += fac * Gq[a] * Gq[b];
                    }
            }

    // Hermitian projection (numerical noise).
    for (int p = 0; p < dim; ++p)
        for (int q = p; q < dim; ++q) {
            const Cx s = AccReal(0.5)
                         * (D[static_cast<std::size_t>(p) * dim + q]
                            + std::conj(
                                D[static_cast<std::size_t>(q) * dim + p]));
            D[static_cast<std::size_t>(p) * dim + q] = s;
            D[static_cast<std::size_t>(q) * dim + p] = std::conj(s);
        }

    return D;
}

// ── Sternheimer conduction-space CG ─────────────────────────────────

int sternheimer_cg(const PWWorkspace& ws, const Cx* occ_wfc, int n_occ,
                   AccReal eps, const Cx* rhs, Cx* x, AccReal tol, int max_iter,
                   fftw_complex* work_buf) {
    const int ng = ws.n_g;

    auto apply_A = [&](const Cx* in, Cx* out) {
        pw_apply_h(ws, in, out, work_buf);
        for (int ig = 0; ig < ng; ++ig)
            out[ig] -= eps * in[ig];
        project_conduction(occ_wfc, n_occ, ng, out);
    };
    auto precond = [&](const Cx* r, Cx* z) {
        for (int ig = 0; ig < ng; ++ig) {
            AccReal denom = ws.T_G[ig] - eps;
            if (denom < AccReal(0.35))
                denom = AccReal(0.35);
            z[ig] = r[ig] / denom;
        }
        project_conduction(occ_wfc, n_occ, ng, z);
    };
    auto dot_re = [&](const Cx* a, const Cx* b) {
        AccReal s = 0;
        for (int ig = 0; ig < ng; ++ig)
            s += (std::conj(a[ig]) * b[ig]).real();
        return s;
    };

    std::vector<Cx> b(static_cast<std::size_t>(ng));
    std::copy(rhs, rhs + ng, b.begin());
    project_conduction(occ_wfc, n_occ, ng, b.data());
    project_conduction(occ_wfc, n_occ, ng, x);

    const AccReal b_norm =
        std::sqrt(std::max(dot_re(b.data(), b.data()), AccReal(0)));
    if (b_norm < AccReal(1e-300)) {
        std::fill(x, x + ng, Cx(0, 0));
        return 0;
    }
    const AccReal stop = tol * b_norm;

    std::vector<Cx> r(static_cast<std::size_t>(ng)),
        z(static_cast<std::size_t>(ng)), p(static_cast<std::size_t>(ng)),
        Ap(static_cast<std::size_t>(ng));

    apply_A(x, Ap.data());
    for (int ig = 0; ig < ng; ++ig)
        r[static_cast<std::size_t>(ig)] =
            b[static_cast<std::size_t>(ig)] - Ap[static_cast<std::size_t>(ig)];
    precond(r.data(), z.data());
    p = z;
    AccReal rz = dot_re(r.data(), z.data());

    for (int it = 1; it <= max_iter; ++it) {
        apply_A(p.data(), Ap.data());
        const AccReal pAp = dot_re(p.data(), Ap.data());
        if (!(pAp > AccReal(0))) {
            ATV_WARN("[pw][dfpt] sternheimer CG lost positive-definiteness "
                     "(pAp={}); returning current iterate",
                     pAp);
            return it;
        }
        const AccReal a = rz / pAp;
        for (int ig = 0; ig < ng; ++ig) {
            x[ig] += a * p[static_cast<std::size_t>(ig)];
            r[static_cast<std::size_t>(ig)] -=
                a * Ap[static_cast<std::size_t>(ig)];
        }
        const AccReal r_norm =
            std::sqrt(std::max(dot_re(r.data(), r.data()), AccReal(0)));
        if (r_norm <= stop) {
            project_conduction(occ_wfc, n_occ, ng, x);
            return it;
        }
        precond(r.data(), z.data());
        const AccReal rz_new = dot_re(r.data(), z.data());
        const AccReal beta = rz_new / rz;
        rz = rz_new;
        for (int ig = 0; ig < ng; ++ig)
            p[static_cast<std::size_t>(ig)] =
                z[static_cast<std::size_t>(ig)]
                + beta * p[static_cast<std::size_t>(ig)];
    }
    project_conduction(occ_wfc, n_occ, ng, x);
    return max_iter + 1;
}

// ── Γ-point DFPT assembly ───────────────────────────────────────────

DfptGammaOutput
dfpt_gamma_force_constants(const cpu::Frame& frame, const PW& method,
                           AccReal response_tol, int response_max_iter,
                           AccReal cg_tol, int cg_max_iter, AccReal mix_beta) {
    if (method.upf_path.empty())
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: requires a UPF pseudopotential "
            "(toy V_loc path unsupported)");
    if (method.xc != cpu::scf::XCType::LSDA)
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: response kernel implemented for "
            "XCType::LSDA only");
    if (method.smearing != Smearing::none)
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: metallic smearing unsupported "
            "(fixed occupations only)");
    if (method.spin_polarized)
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: collinear spin unsupported");
    if (!(response_tol > 0) || !(cg_tol > 0) || response_max_iter < 1
        || cg_max_iter < 1 || !(mix_beta > 0) || mix_beta > 1)
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: invalid response/CG parameters");

    const int n = frame.nlocal();
    if (n <= 0)
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: frame has no atoms");
    const int dim = 3 * n;

    // ── Ground-state SCF (shared minimal loop; honors method.xc) ──
    int n_elec = 0;
    const int* types = frame.span(kAtom_type);
    for (int i = 0; i < n; ++i)
        n_elec += valence_from_z(types[i]);
    const int n_bands = std::max(n_elec / 2, 1);
    if (n_elec % 2 != 0)
        throw std::invalid_argument(
            "dfpt_gamma_force_constants: odd electron count needs "
            "spin/smearing (unsupported)");

    PoissonSolver poisson;
    XCBuilder xc;
    VeffBuilder veff;
    DensityBuilder density{AccReal(method.mix_beta), method.pulay_depth};
    DavidsonSolver davidson{n_bands, AccReal(method.tol) * AccReal(0.1), 40};

    DfptGammaOutput out;
    PWProblem prob;
    const PWResult scf = run_scf_fixed_occ(frame, method, prob, poisson, xc,
                                           veff, density, davidson, "dfpt");
    out.scf_converged = scf.converged;
    if (!scf.converged)
        ATV_WARN("[pw][dfpt] ground-state SCF not converged after {} iters; "
                 "response uses last density",
                 scf.iterations);

    PWWorkspace& ws = prob.ws;
    const GVectorList& gv0 = prob.gv0;
    ws.gvec = &prob.gv0;
    ws.n_g = gv0.count();
    const int ng0 = gv0.count();
    const int nr = ws.nr;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const AccReal omega = ws.omega;
    const AccReal dV = omega / AccReal(nr);
    const AccReal occ_f = AccReal(2);

    FFTWBuffer<fftw_complex> work_buf(static_cast<std::size_t>(nr));
    FFTWBuffer<fftw_complex> work_buf2(static_cast<std::size_t>(nr));

    // ── Freeze the ground-state H: V_eff from the final mixed density ──
    poisson.step(ws);
    xc.step(ws, method.xc);
    veff.step(ws);

    // Ground-state density on the real grid (for f_xc).
    std::vector<AccReal> rho_real(static_cast<std::size_t>(nr));
    {
        std::vector<Cx> rho_g(static_cast<std::size_t>(ng0));
        for (int ig = 0; ig < ng0; ++ig)
            rho_g[static_cast<std::size_t>(ig)] = ws.rho_G[ig];
        sphere_to_real(ws, gv0, rho_g.data(), rho_real.data(), work_buf.data());
    }
    std::vector<AccReal> fxc_r(static_cast<std::size_t>(nr));
    for (int i = 0; i < nr; ++i)
        fxc_r[static_cast<std::size_t>(i)] = lda_fxc_at(
            std::max(rho_real[static_cast<std::size_t>(i)], AccReal(0)),
            method.xc);

    // ── Eigenstates of the frozen H at every k (all ranks, all k) ──
    struct KState {
        std::vector<Cx> wfc;      // n_bands × ng_k (row-major)
        std::vector<AccReal> eps; // n_bands
        int ng = 0;
    };
    std::vector<KState> kstates(static_cast<std::size_t>(n_kpts));
    auto bind_k = [&](int ik) {
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
    };
    for (int ik = 0; ik < n_kpts; ++ik) {
        bind_k(ik);
        davidson.step(ws);
        auto& ks = kstates[static_cast<std::size_t>(ik)];
        ks.ng = ws.n_g;
        ks.eps.assign(ws.eigenvalues.begin(), ws.eigenvalues.begin() + n_bands);
        ks.wfc.resize(static_cast<std::size_t>(n_bands) * ws.n_g);
        for (int b = 0; b < n_bands; ++b)
            for (int ig = 0; ig < ws.n_g; ++ig)
                ks.wfc[static_cast<std::size_t>(b) * ws.n_g + ig] =
                    ws.wfc(b, ig);
    }
    ws.gvec = &prob.gv0;
    ws.n_g = ng0;

    // ── Per-atom radial form factor v(G) on the density sphere ──
    const UPFData upf = read_upf_v1(method.upf_path);
    const std::vector<AccReal> vg = upf_vloc_to_gspace(upf, gv0, omega);

    std::vector<Vec3> pos_ang_v(static_cast<std::size_t>(n));
    {
        const StorageReal* sx = frame.span(kAtom_x);
        const StorageReal* sy = frame.span(kAtom_y);
        const StorageReal* sz = frame.span(kAtom_z);
        for (int i = 0; i < n; ++i)
            pos_ang_v[static_cast<std::size_t>(i)] = {
                static_cast<AccReal>(sx[i]), static_cast<AccReal>(sy[i]),
                static_cast<AccReal>(sz[i])};
    }

    // Bare local perturbation potentials ΔV_loc^{(I,α)}(r) [Ha/Bohr].
    std::vector<std::vector<AccReal>> dvloc_r(
        static_cast<std::size_t>(dim),
        std::vector<AccReal>(static_cast<std::size_t>(nr), AccReal(0)));
    {
        std::vector<Cx> dvg(static_cast<std::size_t>(ng0));
        for (int ia = 0; ia < n; ++ia) {
            const AccReal ax =
                ws.atom_pos_bohr[static_cast<std::size_t>(ia)][0];
            const AccReal ay =
                ws.atom_pos_bohr[static_cast<std::size_t>(ia)][1];
            const AccReal az =
                ws.atom_pos_bohr[static_cast<std::size_t>(ia)][2];
            for (int a = 0; a < 3; ++a) {
                for (int ig = 0; ig < ng0; ++ig) {
                    const AccReal ga = gv0.g[ig][a];
                    const AccReal phase = gv0.g[ig][0] * ax + gv0.g[ig][1] * ay
                                          + gv0.g[ig][2] * az;
                    // −i G_α v(G) e^{−iG·R_I}
                    const Cx eir(std::cos(phase), -std::sin(phase));
                    dvg[static_cast<std::size_t>(ig)] =
                        Cx(0, -1) * ga * vg[static_cast<std::size_t>(ig)] * eir;
                }
                sphere_to_real(
                    ws, gv0, dvg.data(),
                    dvloc_r[static_cast<std::size_t>(3 * ia + a)].data(),
                    work_buf.data());
            }
        }
    }

    // ── Ion–ion Ewald Hessian ──
    std::vector<AccReal> valence(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        valence[static_cast<std::size_t>(i)] =
            static_cast<AccReal>(valence_from_z(types[i]));
    const AccReal ecut_ht = static_cast<AccReal>(method.ecut_ry) * AccReal(0.5);
    std::vector<AccReal> C =
        ewald_dynamical_matrix(pos_ang_v, valence, frame.simbox, ecut_ht);

    auto Cat = [&](int mu, int nu) -> AccReal& {
        return C[static_cast<std::size_t>(mu) * dim + nu];
    };

    // ── Expectation terms (same atom): local ∂²V_loc + KB ∂² ──
    for (int ia = 0; ia < n; ++ia) {
        const AccReal ax = ws.atom_pos_bohr[static_cast<std::size_t>(ia)][0];
        const AccReal ay = ws.atom_pos_bohr[static_cast<std::size_t>(ia)][1];
        const AccReal az = ws.atom_pos_bohr[static_cast<std::size_t>(ia)][2];
        for (int ig = 0; ig < ng0; ++ig) {
            const AccReal phase =
                gv0.g[ig][0] * ax + gv0.g[ig][1] * ay + gv0.g[ig][2] * az;
            const Cx eir(std::cos(phase), -std::sin(phase));
            const AccReal re = (std::conj(ws.rho_G[ig])
                                * (vg[static_cast<std::size_t>(ig)] * eir))
                                   .real();
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    Cat(3 * ia + a, 3 * ia + b) -=
                        omega * gv0.g[ig][a] * gv0.g[ig][b] * re;
        }
    }

    {
        std::vector<KbMoments> mo;
        for (int ik = 0; ik < n_kpts; ++ik) {
            bind_k(ik);
            const auto& ks = kstates[static_cast<std::size_t>(ik)];
            const AccReal w = prob.kg.weights[static_cast<std::size_t>(ik)];
            for (int b = 0; b < n_bands; ++b) {
                const Cx* psi =
                    ks.wfc.data() + static_cast<std::ptrdiff_t>(b) * ks.ng;
                for (std::size_t ia = 0; ia < ws.kb_proj.size(); ++ia) {
                    for (std::size_t ip = 0; ip < ws.kb_proj[ia].size(); ++ip) {
                        const AccReal Dc = ws.kb_proj[ia][ip].d;
                        if (std::abs(Dc) < AccReal(1e-16))
                            continue;
                        kb_moments(ws, ia, ip, psi, /*want_t=*/true, mo);
                        for (const auto& mm : mo) {
                            for (int a = 0; a < 3; ++a)
                                for (int bb = 0; bb < 3; ++bb) {
                                    const AccReal term =
                                        w * occ_f * Dc
                                        * (AccReal(-2)
                                               * (mm.t[a][bb] * std::conj(mm.c))
                                                     .real()
                                           + AccReal(2)
                                                 * (mm.s[a]
                                                    * std::conj(mm.s[bb]))
                                                       .real());
                                    Cat(3 * static_cast<int>(ia) + a,
                                        3 * static_cast<int>(ia) + bb) += term;
                                }
                        }
                    }
                }
            }
        }
        ws.gvec = &prob.gv0;
        ws.n_g = ng0;
    }

    // ── Self-consistent response, one perturbation μ = (J,β) at a time ──
    bool all_responses_converged = true;

    // Persistent Δψ storage: per (k, band) for the current μ.
    std::vector<std::vector<Cx>> dpsi(static_cast<std::size_t>(n_kpts)
                                      * static_cast<std::size_t>(n_bands));

    std::vector<AccReal> dn_mixed(static_cast<std::size_t>(nr));
    std::vector<AccReal> dn_new(static_cast<std::size_t>(nr));
    std::vector<AccReal> dv_hxc_r(static_cast<std::size_t>(nr));
    std::vector<AccReal> dv_scr_r(static_cast<std::size_t>(nr));
    std::vector<Cx> dn_g(static_cast<std::size_t>(ng0));
    std::vector<Cx> rhs, hpsi_tmp;

    for (int mu = 0; mu < dim; ++mu) {
        std::fill(dn_mixed.begin(), dn_mixed.end(), AccReal(0));
        std::fill(dv_hxc_r.begin(), dv_hxc_r.end(), AccReal(0));
        for (auto& v : dpsi)
            v.clear();
        const int atom_mu = mu / 3;
        const int axis_mu = mu % 3;

        bool converged = false;
        for (int rit = 0; rit < response_max_iter; ++rit) {
            // Screened perturbation on the grid.
            for (int i = 0; i < nr; ++i)
                dv_scr_r[static_cast<std::size_t>(i)] =
                    dvloc_r[static_cast<std::size_t>(mu)]
                           [static_cast<std::size_t>(i)]
                    + dv_hxc_r[static_cast<std::size_t>(i)];

            std::fill(dn_new.begin(), dn_new.end(), AccReal(0));

            for (int ik = 0; ik < n_kpts; ++ik) {
                bind_k(ik);
                const auto& ks = kstates[static_cast<std::size_t>(ik)];
                const int ngk = ks.ng;
                const AccReal w = prob.kg.weights[static_cast<std::size_t>(ik)];
                rhs.resize(static_cast<std::size_t>(ngk));
                hpsi_tmp.resize(static_cast<std::size_t>(ngk));

                for (int b = 0; b < n_bands; ++b) {
                    const Cx* psi =
                        ks.wfc.data() + static_cast<std::ptrdiff_t>(b) * ngk;
                    // rhs = −(ΔV_scr ψ + ΔV_nl^{μ} ψ)
                    apply_mult_potential(ws, dv_scr_r.data(), psi, rhs.data(),
                                         work_buf.data());
                    add_dnl_psi(ws, static_cast<std::size_t>(atom_mu), axis_mu,
                                psi, rhs.data());
                    for (int ig = 0; ig < ngk; ++ig)
                        rhs[static_cast<std::size_t>(ig)] =
                            -rhs[static_cast<std::size_t>(ig)];

                    auto& x = dpsi[static_cast<std::size_t>(ik) * n_bands
                                   + static_cast<std::size_t>(b)];
                    if (static_cast<int>(x.size()) != ngk)
                        x.assign(static_cast<std::size_t>(ngk), Cx(0, 0));
                    const int its = sternheimer_cg(
                        ws, ks.wfc.data(), n_bands, ks.eps[b], rhs.data(),
                        x.data(), cg_tol, cg_max_iter, work_buf2.data());
                    if (its > cg_max_iter)
                        ATV_WARN("[pw][dfpt] CG not converged (mu={} ik={} "
                                 "band={})",
                                 mu, ik, b);

                    // Δn(r) += w f 2 Re[ψ*(r) Δψ(r)] / Ω  (FFT-grid ψ).
                    pad_sphere_to_grid(ws, *ws.gvec, ngk, psi, work_buf.data());
                    ws.fft->c2c_bwd(work_buf.data());
                    pad_sphere_to_grid(ws, *ws.gvec, ngk, x.data(),
                                       work_buf2.data());
                    ws.fft->c2c_bwd(work_buf2.data());
                    const AccReal pref = w * occ_f * AccReal(2) / omega;
                    for (int i = 0; i < nr; ++i) {
                        const AccReal re =
                            static_cast<AccReal>(work_buf.data()[i][0])
                                * static_cast<AccReal>(work_buf2.data()[i][0])
                            + static_cast<AccReal>(work_buf.data()[i][1])
                                  * static_cast<AccReal>(
                                      work_buf2.data()[i][1]);
                        dn_new[static_cast<std::size_t>(i)] += pref * re;
                    }
                }
            }
            ws.gvec = &prob.gv0;
            ws.n_g = ng0;

            // Convergence on ‖Δn_new − Δn_mixed‖ (grid L2, volume weight).
            AccReal diff2 = 0;
            for (int i = 0; i < nr; ++i) {
                const AccReal d = dn_new[static_cast<std::size_t>(i)]
                                  - dn_mixed[static_cast<std::size_t>(i)];
                diff2 += d * d;
            }
            const AccReal delta = std::sqrt(diff2 * dV);
            const bool this_converged = (delta < response_tol);

            // Linear mix, then rebuild ΔV_HXC from the mixed Δn.
            for (int i = 0; i < nr; ++i)
                dn_mixed[static_cast<std::size_t>(i)] +=
                    mix_beta
                    * (dn_new[static_cast<std::size_t>(i)]
                       - dn_mixed[static_cast<std::size_t>(i)]);

            real_to_sphere(ws, gv0, dn_mixed.data(), dn_g.data(),
                           work_buf.data());
            // ΔV_H(G) = 4π Δρ(G)/G² (G = 0 → 0).
            for (int ig = 0; ig < ng0; ++ig) {
                const AccReal g2 = gv0.g[ig].norm_squared();
                if (g2 < AccReal(1e-16))
                    dn_g[static_cast<std::size_t>(ig)] = Cx(0, 0);
                else
                    dn_g[static_cast<std::size_t>(ig)] *=
                        AccReal(4) * AccReal(units::pi) / g2;
            }
            sphere_to_real(ws, gv0, dn_g.data(), dv_hxc_r.data(),
                           work_buf.data());
            // + f_xc Δn, sphere-filtered like the ground-state V_xc.
            {
                std::vector<AccReal> tmp(static_cast<std::size_t>(nr));
                for (int i = 0; i < nr; ++i)
                    tmp[static_cast<std::size_t>(i)] =
                        fxc_r[static_cast<std::size_t>(i)]
                        * dn_mixed[static_cast<std::size_t>(i)];
                real_to_sphere(ws, gv0, tmp.data(), dn_g.data(),
                               work_buf.data());
                std::vector<AccReal> tmp2(static_cast<std::size_t>(nr));
                sphere_to_real(ws, gv0, dn_g.data(), tmp2.data(),
                               work_buf.data());
                for (int i = 0; i < nr; ++i)
                    dv_hxc_r[static_cast<std::size_t>(i)] +=
                        tmp2[static_cast<std::size_t>(i)];
            }

            if (this_converged) {
                converged = true;
                ATV_DEBUG("[pw][dfpt] response mu={} converged iter={} "
                          "delta={:.3e}",
                          mu, rit + 1, delta);
                break;
            }
        }
        if (!converged) {
            all_responses_converged = false;
            ATV_WARN("[pw][dfpt] response mu={} not converged after {} iters",
                     mu, response_max_iter);
        }

        // ── Response contributions C[ν][μ] ──
        // Local: ∫ Δn^{(μ)}(r) ΔV_loc^{(ν)}(r) dr.
        for (int nu = 0; nu < dim; ++nu) {
            AccReal s = 0;
            const auto& dvl = dvloc_r[static_cast<std::size_t>(nu)];
            for (int i = 0; i < nr; ++i)
                s += dn_mixed[static_cast<std::size_t>(i)]
                     * dvl[static_cast<std::size_t>(i)];
            Cat(nu, mu) += s * dV;
        }

        // Nonlocal: Σ_{k,i} w f 2 Re⟨Δψ^{(μ)}|∂V_nl^{(ν)}|ψ⟩.
        {
            std::vector<KbMoments> mo_psi, mo_dpsi;
            for (int ik = 0; ik < n_kpts; ++ik) {
                bind_k(ik);
                const auto& ks = kstates[static_cast<std::size_t>(ik)];
                const AccReal w = prob.kg.weights[static_cast<std::size_t>(ik)];
                for (int b = 0; b < n_bands; ++b) {
                    const Cx* psi =
                        ks.wfc.data() + static_cast<std::ptrdiff_t>(b) * ks.ng;
                    const auto& x = dpsi[static_cast<std::size_t>(ik) * n_bands
                                         + static_cast<std::size_t>(b)];
                    for (std::size_t ja = 0; ja < ws.kb_proj.size(); ++ja) {
                        for (std::size_t ip = 0; ip < ws.kb_proj[ja].size();
                             ++ip) {
                            const AccReal Dc = ws.kb_proj[ja][ip].d;
                            if (std::abs(Dc) < AccReal(1e-16))
                                continue;
                            kb_moments(ws, ja, ip, psi, false, mo_psi);
                            kb_moments(ws, ja, ip, x.data(), false, mo_dpsi);
                            for (std::size_t m = 0; m < mo_psi.size(); ++m) {
                                const Cx cP = mo_psi[m].c;
                                const Cx cD = mo_dpsi[m].c;
                                for (int bb = 0; bb < 3; ++bb) {
                                    // i D (conj(cΔ) s_b − conj(sΔ_b) c)
                                    const Cx el =
                                        Cx(0, 1) * Dc
                                        * (std::conj(cD) * mo_psi[m].s[bb]
                                           - std::conj(mo_dpsi[m].s[bb]) * cP);
                                    Cat(3 * static_cast<int>(ja) + bb, mu) +=
                                        w * occ_f * AccReal(2) * el.real();
                                }
                            }
                        }
                    }
                }
            }
            ws.gvec = &prob.gv0;
            ws.n_g = ng0;
        }
    }
    out.response_converged = all_responses_converged;

    // ── Symmetrize; record the raw asymmetry as a diagnostic ──
    AccReal max_asym = 0;
    for (int i = 0; i < dim; ++i)
        for (int j = i + 1; j < dim; ++j) {
            const AccReal d = std::abs(Cat(i, j) - Cat(j, i)) * AccReal(0.5);
            if (d > max_asym)
                max_asym = d;
            const AccReal avg = AccReal(0.5) * (Cat(i, j) + Cat(j, i));
            Cat(i, j) = avg;
            Cat(j, i) = avg;
        }
    out.max_asymmetry = max_asym;
    out.force_constants = std::move(C);
    return out;
}

// ── Finite-q DFPT assembly ──────────────────────────────────────────

namespace {

/// First-order nonlocal KB at finite phonon q: writes into the k+q sphere.
/// Moments of ψ are taken on the k sphere; projectors for the ket on k+q.
/// Reduces to `add_dnl_psi` when kq == k (same sphere / same projectors).
void add_dnl_psi_q(const PWWorkspace& ws_k, const PWWorkspace& ws_kq,
                   std::size_t ia, int alpha, const Cx* psi_k, Cx* out_kq) {
    if (ws_k.kb_proj.empty() || ia >= ws_k.kb_proj.size())
        return;
    if (ws_kq.kb_proj.empty() || ia >= ws_kq.kb_proj.size())
        return;
    const auto& gv_kq = *ws_kq.gvec;
    const AccReal pref = AccReal(4.0 * M_PI) / std::sqrt(ws_kq.omega);
    const AccReal ax = ws_kq.atom_pos_bohr[ia][0];
    const AccReal ay = ws_kq.atom_pos_bohr[ia][1];
    const AccReal az = ws_kq.atom_pos_bohr[ia][2];

    std::vector<KbMoments> mo;
    const std::size_t nproj =
        std::min(ws_k.kb_proj[ia].size(), ws_kq.kb_proj[ia].size());
    for (std::size_t ip = 0; ip < nproj; ++ip) {
        const auto& proj_k = ws_k.kb_proj[ia][ip];
        const auto& proj_kq = ws_kq.kb_proj[ia][ip];
        const AccReal D = proj_k.d;
        if (std::abs(D) < AccReal(1e-16))
            continue;
        const int l = proj_k.l;
        const int nm = n_m_for_l(l);
        const auto im_l = im_l_factor(l);
        // Moments of ψ at k (bra projectors).
        {
            // Temporarily bind: kb_moments reads ws.gvec / kb_proj / n_g.
            // Call with ws_k.
            kb_moments(ws_k, ia, ip, psi_k, /*want_t=*/false, mo);
        }

        std::vector<AccReal> ylm(static_cast<std::size_t>(nm));
        for (int ig = 0; ig < ws_kq.n_g; ++ig) {
            const AccReal fl = proj_kq.f_G[ig];
            if (std::abs(fl) < AccReal(1e-30))
                continue;
            const AccReal qx = gv_kq.g[ig][0] + gv_kq.kvec[0];
            const AccReal qy = gv_kq.g[ig][1] + gv_kq.kvec[1];
            const AccReal qz = gv_kq.g[ig][2] + gv_kq.kvec[2];
            const AccReal qr = std::sqrt(qx * qx + qy * qy + qz * qz);
            const AccReal q[3] = {qx, qy, qz};
            const AccReal phase = qx * ax + qy * ay + qz * az;
            const AccReal cos_p = std::cos(phase), sin_p = std::sin(phase);
            compute_ylm(qx, qy, qz, qr, l, ylm.data());
            for (int m = 0; m < nm; ++m) {
                const Cx beta = pref * im_l * ylm[static_cast<std::size_t>(m)]
                                * fl * Cx(cos_p, -sin_p);
                const auto& mm = mo[static_cast<std::size_t>(m)];
                // i D β_kq (s_α^k − (k+q+G)_α c^k) — Γ-matching two-term form.
                out_kq[ig] +=
                    Cx(0, 1) * D * beta * (mm.s[alpha] - q[alpha] * mm.c);
            }
        }
    }
}

/// Complex field on the density sphere → complex real-space grid.
void sphere_to_complex_grid(const PWWorkspace& ws, const GVectorList& gv0,
                            const Cx* vg, Cx* vr, fftw_complex* work) {
    pad_sphere_to_grid(ws, gv0, gv0.count(), vg, work);
    ws.fft->c2c_bwd(work);
    for (int i = 0; i < ws.nr; ++i)
        vr[i] = Cx(static_cast<AccReal>(work[i][0]),
                   static_cast<AccReal>(work[i][1]));
}

/// Complex real-space field → density-sphere Fourier coefficients (1/N).
void complex_grid_to_sphere(const PWWorkspace& ws, const GVectorList& gv0,
                            const Cx* fr, Cx* fg, fftw_complex* work) {
    for (int i = 0; i < ws.nr; ++i) {
        work[i][0] = static_cast<double>(fr[i].real());
        work[i][1] = static_cast<double>(fr[i].imag());
    }
    ws.fft->c2c_fwd(work);
    extract_sphere_from_grid(ws, gv0, gv0.count(), work, fg);
    const AccReal inv_nr = AccReal(1) / AccReal(ws.nr);
    for (int ig = 0; ig < gv0.count(); ++ig)
        fg[ig] *= inv_nr;
}

} // namespace

DfptQOutput dfpt_q_force_constants(const cpu::Frame& frame, const PW& method,
                                   const Vec3& q_frac, AccReal response_tol,
                                   int response_max_iter, AccReal cg_tol,
                                   int cg_max_iter, AccReal mix_beta) {
    if (method.upf_path.empty())
        throw std::invalid_argument(
            "dfpt_q_force_constants: requires a UPF pseudopotential "
            "(toy V_loc path unsupported)");
    if (method.xc != cpu::scf::XCType::LSDA)
        throw std::invalid_argument(
            "dfpt_q_force_constants: response kernel implemented for "
            "XCType::LSDA only");
    if (method.smearing != Smearing::none)
        throw std::invalid_argument(
            "dfpt_q_force_constants: metallic smearing unsupported "
            "(fixed occupations only)");
    if (method.spin_polarized)
        throw std::invalid_argument(
            "dfpt_q_force_constants: collinear spin unsupported");
    if (!(response_tol > 0) || !(cg_tol > 0) || response_max_iter < 1
        || cg_max_iter < 1 || !(mix_beta > 0) || mix_beta > 1)
        throw std::invalid_argument(
            "dfpt_q_force_constants: invalid response/CG parameters");

    const int n = frame.nlocal();
    if (n <= 0)
        throw std::invalid_argument(
            "dfpt_q_force_constants: frame has no atoms");
    const int dim = 3 * n;

    int n_elec = 0;
    const int* types = frame.span(kAtom_type);
    for (int i = 0; i < n; ++i)
        n_elec += valence_from_z(types[i]);
    const int n_bands = std::max(n_elec / 2, 1);
    if (n_elec % 2 != 0)
        throw std::invalid_argument(
            "dfpt_q_force_constants: odd electron count needs "
            "spin/smearing (unsupported)");

    PoissonSolver poisson;
    XCBuilder xc;
    VeffBuilder veff;
    DensityBuilder density{AccReal(method.mix_beta), method.pulay_depth};
    DavidsonSolver davidson{n_bands, AccReal(method.tol) * AccReal(0.1), 40};

    DfptQOutput out;
    PWProblem prob;
    const PWResult scf = run_scf_fixed_occ(frame, method, prob, poisson, xc,
                                           veff, density, davidson, "dfpt_q");
    out.scf_converged = scf.converged;
    if (!scf.converged)
        ATV_WARN("[pw][dfpt_q] ground-state SCF not converged after {} iters; "
                 "response uses last density",
                 scf.iterations);

    PWWorkspace& ws = prob.ws;
    const GVectorList& gv0 = prob.gv0;
    ws.gvec = &prob.gv0;
    ws.n_g = gv0.count();
    const int ng0 = gv0.count();
    const int nr = ws.nr;
    const int n_kpts = static_cast<int>(prob.kg.kpoints.size());
    const AccReal omega = ws.omega;
    const AccReal dV = omega / AccReal(nr);
    const AccReal occ_f = AccReal(2);

    // Cartesian q [1/Bohr].
    cpu::Box cell_bohr;
    cell_bohr.set_cell(frame.simbox.cell() * AccReal(1.8897259886));
    const auto recip = cell_bohr.recip_cell();
    const Vec3 q_cart{recip(0, 0) * q_frac[0] + recip(1, 0) * q_frac[1]
                          + recip(2, 0) * q_frac[2],
                      recip(0, 1) * q_frac[0] + recip(1, 1) * q_frac[1]
                          + recip(2, 1) * q_frac[2],
                      recip(0, 2) * q_frac[0] + recip(1, 2) * q_frac[1]
                          + recip(2, 2) * q_frac[2]};

    FFTWBuffer<fftw_complex> work_buf(static_cast<std::size_t>(nr));
    FFTWBuffer<fftw_complex> work_buf2(static_cast<std::size_t>(nr));

    poisson.step(ws);
    xc.step(ws, method.xc);
    veff.step(ws);

    // Ground-state density + f_xc on the real grid.
    std::vector<AccReal> rho_real(static_cast<std::size_t>(nr));
    {
        std::vector<Cx> rho_g(static_cast<std::size_t>(ng0));
        for (int ig = 0; ig < ng0; ++ig)
            rho_g[static_cast<std::size_t>(ig)] = ws.rho_G[ig];
        sphere_to_real(ws, gv0, rho_g.data(), rho_real.data(), work_buf.data());
    }
    std::vector<AccReal> fxc_r(static_cast<std::size_t>(nr));
    for (int i = 0; i < nr; ++i)
        fxc_r[static_cast<std::size_t>(i)] = lda_fxc_at(
            std::max(rho_real[static_cast<std::size_t>(i)], AccReal(0)),
            method.xc);

    // Freeze V_eff_r for all subsequent Davidson / Sternheimer calls.
    std::vector<AccReal> veff_frozen(static_cast<std::size_t>(nr));
    for (int i = 0; i < nr; ++i)
        veff_frozen[static_cast<std::size_t>(i)] = ws.V_eff_r[i];

    struct KState {
        std::vector<Cx> wfc;
        std::vector<AccReal> eps;
        int ng = 0;
        std::vector<std::vector<PWWorkspace::KBProjG>> kb_proj;
    };

    auto bind_state = [&](KState& ks, const GVectorList& gvk) {
        ws.gvec = &gvk;
        ws.n_g = gvk.count();
        ws.T_G.resize({static_cast<std::size_t>(ws.n_g)});
        for (int ig = 0; ig < ws.n_g; ++ig)
            ws.T_G[ig] = gvk.ekin[ig];
        ws.kb_proj = ks.kb_proj;
        for (int i = 0; i < nr; ++i)
            ws.V_eff_r[i] = veff_frozen[static_cast<std::size_t>(i)];
    };

    // Eigenstates at every ground-state k; keep owned G-sphere copies.
    std::vector<GVectorList> gv_k0 = prob.gv_k;
    std::vector<KState> kstates(static_cast<std::size_t>(n_kpts));
    for (int ik = 0; ik < n_kpts; ++ik) {
        auto& gvk = gv_k0[static_cast<std::size_t>(ik)];
        auto& ks = kstates[static_cast<std::size_t>(ik)];
        if (!prob.kb_proj_k.empty())
            ks.kb_proj = prob.kb_proj_k[static_cast<std::size_t>(ik)];
        bind_state(ks, gvk);
        if (!prob.atom_wfc_k.empty())
            ws.atom_wfc = prob.atom_wfc_k[static_cast<std::size_t>(ik)];
        davidson.step(ws);
        ks.ng = ws.n_g;
        ks.eps.assign(ws.eigenvalues.begin(), ws.eigenvalues.begin() + n_bands);
        ks.wfc.resize(static_cast<std::size_t>(n_bands) * ws.n_g);
        for (int b = 0; b < n_bands; ++b)
            for (int ig = 0; ig < ws.n_g; ++ig)
                ks.wfc[static_cast<std::size_t>(b) * ws.n_g + ig] =
                    ws.wfc(b, ig);
    }
    const std::vector<AccReal> k_weights = prob.kg.weights;

    // Eigenstates at k+q (occupied manifold for P_c^{k+q}).
    std::vector<Vec3> kq_frac(static_cast<std::size_t>(n_kpts));
    for (int ik = 0; ik < n_kpts; ++ik) {
        const auto& kf = prob.kg.kpoints[static_cast<std::size_t>(ik)];
        kq_frac[static_cast<std::size_t>(ik)] = {
            kf[0] + q_frac[0], kf[1] + q_frac[1], kf[2] + q_frac[2]};
    }
    rebuild_k_tables(prob, frame, method.upf_path, method.ecut_ry, kq_frac);
    std::vector<GVectorList> gv_kq = prob.gv_k;
    std::vector<KState> kqstates(static_cast<std::size_t>(n_kpts));
    for (int ik = 0; ik < n_kpts; ++ik) {
        auto& gvk = gv_kq[static_cast<std::size_t>(ik)];
        auto& ks = kqstates[static_cast<std::size_t>(ik)];
        if (!prob.kb_proj_k.empty())
            ks.kb_proj = prob.kb_proj_k[static_cast<std::size_t>(ik)];
        bind_state(ks, gvk);
        if (!prob.atom_wfc_k.empty())
            ws.atom_wfc = prob.atom_wfc_k[static_cast<std::size_t>(ik)];
        davidson.step(ws);
        ks.ng = ws.n_g;
        ks.eps.assign(ws.eigenvalues.begin(), ws.eigenvalues.begin() + n_bands);
        ks.wfc.resize(static_cast<std::size_t>(n_bands) * ws.n_g);
        for (int b = 0; b < n_bands; ++b)
            for (int ig = 0; ig < ws.n_g; ++ig)
                ks.wfc[static_cast<std::size_t>(b) * ws.n_g + ig] =
                    ws.wfc(b, ig);
    }

    // Density sphere / V_eff stay on Γ; k-sphere tables are the owned copies.
    ws.gvec = &prob.gv0;
    ws.n_g = ng0;

    // v(|G+q|) on the density sphere: shift g → G+q for the radial transform.
    const UPFData upf = read_upf_v1(method.upf_path);
    GVectorList gv_qp = gv0;
    for (int ig = 0; ig < ng0; ++ig) {
        gv_qp.g[static_cast<std::size_t>(ig)] = {
            gv0.g[ig][0] + q_cart[0], gv0.g[ig][1] + q_cart[1],
            gv0.g[ig][2] + q_cart[2]};
    }
    const std::vector<AccReal> vg_qp = upf_vloc_to_gspace(upf, gv_qp, omega);

    std::vector<Vec3> pos_ang_v(static_cast<std::size_t>(n));
    {
        const StorageReal* sx = frame.span(kAtom_x);
        const StorageReal* sy = frame.span(kAtom_y);
        const StorageReal* sz = frame.span(kAtom_z);
        for (int i = 0; i < n; ++i)
            pos_ang_v[static_cast<std::size_t>(i)] = {
                static_cast<AccReal>(sx[i]), static_cast<AccReal>(sy[i]),
                static_cast<AccReal>(sz[i])};
    }

    // Bare local ΔV^{Iα}_q(r) — complex lattice-periodic part.
    std::vector<std::vector<Cx>> dvloc_r(
        static_cast<std::size_t>(dim),
        std::vector<Cx>(static_cast<std::size_t>(nr), Cx(0, 0)));
    {
        std::vector<Cx> dvg(static_cast<std::size_t>(ng0));
        for (int ia = 0; ia < n; ++ia) {
            const AccReal ax =
                ws.atom_pos_bohr[static_cast<std::size_t>(ia)][0];
            const AccReal ay =
                ws.atom_pos_bohr[static_cast<std::size_t>(ia)][1];
            const AccReal az =
                ws.atom_pos_bohr[static_cast<std::size_t>(ia)][2];
            for (int a = 0; a < 3; ++a) {
                for (int ig = 0; ig < ng0; ++ig) {
                    const AccReal gqa = gv0.g[ig][a] + q_cart[a];
                    const AccReal phase =
                        (gv0.g[ig][0] + q_cart[0]) * ax
                        + (gv0.g[ig][1] + q_cart[1]) * ay
                        + (gv0.g[ig][2] + q_cart[2]) * az;
                    const Cx eir(std::cos(phase), -std::sin(phase));
                    dvg[static_cast<std::size_t>(ig)] =
                        Cx(0, -1) * gqa * vg_qp[static_cast<std::size_t>(ig)]
                        * eir;
                }
                sphere_to_complex_grid(
                    ws, gv0, dvg.data(),
                    dvloc_r[static_cast<std::size_t>(3 * ia + a)].data(),
                    work_buf.data());
            }
        }
    }

    // Ion–ion Ewald C(q).
    std::vector<AccReal> valence(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        valence[static_cast<std::size_t>(i)] =
            static_cast<AccReal>(valence_from_z(types[i]));
    const AccReal ecut_ht = static_cast<AccReal>(method.ecut_ry) * AccReal(0.5);
    std::vector<Cx> C =
        ewald_dynamical_matrix_q(pos_ang_v, valence, frame.simbox, ecut_ht,
                                 q_frac);

    auto Cat = [&](int mu, int nu) -> Cx& {
        return C[static_cast<std::size_t>(mu) * dim + nu];
    };

    // Bare local second-derivative (same atom): −Ω (G+q)_a (G+q)_b Re/complex
    // contraction with ρ*(G) v(|G+q|) e^{-i(G+q)·R}.
    for (int ia = 0; ia < n; ++ia) {
        const AccReal ax = ws.atom_pos_bohr[static_cast<std::size_t>(ia)][0];
        const AccReal ay = ws.atom_pos_bohr[static_cast<std::size_t>(ia)][1];
        const AccReal az = ws.atom_pos_bohr[static_cast<std::size_t>(ia)][2];
        for (int ig = 0; ig < ng0; ++ig) {
            const AccReal gqx = gv0.g[ig][0] + q_cart[0];
            const AccReal gqy = gv0.g[ig][1] + q_cart[1];
            const AccReal gqz = gv0.g[ig][2] + q_cart[2];
            const AccReal gq[3] = {gqx, gqy, gqz};
            const AccReal phase = gqx * ax + gqy * ay + gqz * az;
            const Cx eir(std::cos(phase), -std::sin(phase));
            // Use ground-state ρ(G); contract as conj(ρ)*v*eir (complex).
            const Cx term = std::conj(Cx(ws.rho_G[ig]))
                            * (vg_qp[static_cast<std::size_t>(ig)] * eir);
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    Cat(3 * ia + a, 3 * ia + b) -=
                        omega * gq[a] * gq[b] * term;
        }
    }

    // KB second-derivative diagonal (evaluated at each k; q-independent body).
    {
        std::vector<KbMoments> mo;
        for (int ik = 0; ik < n_kpts; ++ik) {
            auto& ks = kstates[static_cast<std::size_t>(ik)];
            bind_state(ks, gv_k0[static_cast<std::size_t>(ik)]);
            const AccReal w = k_weights[static_cast<std::size_t>(ik)];
            for (int b = 0; b < n_bands; ++b) {
                const Cx* psi =
                    ks.wfc.data() + static_cast<std::ptrdiff_t>(b) * ks.ng;
                for (std::size_t ia = 0; ia < ws.kb_proj.size(); ++ia) {
                    for (std::size_t ip = 0; ip < ws.kb_proj[ia].size(); ++ip) {
                        const AccReal Dc = ws.kb_proj[ia][ip].d;
                        if (std::abs(Dc) < AccReal(1e-16))
                            continue;
                        kb_moments(ws, ia, ip, psi, /*want_t=*/true, mo);
                        for (const auto& mm : mo) {
                            for (int a = 0; a < 3; ++a)
                                for (int bb = 0; bb < 3; ++bb) {
                                    const AccReal term =
                                        w * occ_f * Dc
                                        * (AccReal(-2)
                                               * (mm.t[a][bb] * std::conj(mm.c))
                                                     .real()
                                           + AccReal(2)
                                                 * (mm.s[a]
                                                    * std::conj(mm.s[bb]))
                                                       .real());
                                    Cat(3 * static_cast<int>(ia) + a,
                                        3 * static_cast<int>(ia) + bb) += term;
                                }
                        }
                    }
                }
            }
        }
        ws.gvec = &prob.gv0;
        ws.n_g = ng0;
    }

    // Self-consistent response, one perturbation μ = (J,β) at a time.
    bool all_responses_converged = true;
    std::vector<std::vector<Cx>> dpsi(static_cast<std::size_t>(n_kpts)
                                      * static_cast<std::size_t>(n_bands));
    std::vector<Cx> dn_mixed(static_cast<std::size_t>(nr), Cx(0, 0));
    std::vector<Cx> dn_new(static_cast<std::size_t>(nr));
    std::vector<Cx> dv_hxc_r(static_cast<std::size_t>(nr), Cx(0, 0));
    std::vector<Cx> dv_scr_r(static_cast<std::size_t>(nr));
    std::vector<Cx> dn_g(static_cast<std::size_t>(ng0));
    std::vector<Cx> rhs;

    for (int mu = 0; mu < dim; ++mu) {
        std::fill(dn_mixed.begin(), dn_mixed.end(), Cx(0, 0));
        std::fill(dv_hxc_r.begin(), dv_hxc_r.end(), Cx(0, 0));
        for (auto& v : dpsi)
            v.clear();
        const int atom_mu = mu / 3;
        const int axis_mu = mu % 3;
        bool converged = false;

        for (int rit = 0; rit < response_max_iter; ++rit) {
            for (int i = 0; i < nr; ++i)
                dv_scr_r[static_cast<std::size_t>(i)] =
                    dvloc_r[static_cast<std::size_t>(mu)]
                           [static_cast<std::size_t>(i)]
                    + dv_hxc_r[static_cast<std::size_t>(i)];
            std::fill(dn_new.begin(), dn_new.end(), Cx(0, 0));

            for (int ik = 0; ik < n_kpts; ++ik) {
                auto& ks = kstates[static_cast<std::size_t>(ik)];
                auto& kqs = kqstates[static_cast<std::size_t>(ik)];
                const AccReal w = k_weights[static_cast<std::size_t>(ik)];
                const int ngk = ks.ng;
                const int ngkq = kqs.ng;

                for (int b = 0; b < n_bands; ++b) {
                    const Cx* psi =
                        ks.wfc.data() + static_cast<std::ptrdiff_t>(b) * ngk;

                    // ΔV_scr ψ on real grid from ψ_k, extract on k+q sphere.
                    rhs.assign(static_cast<std::size_t>(ngkq), Cx(0, 0));
                    {
                        ws.gvec = &gv_k0[static_cast<std::size_t>(ik)];
                        ws.n_g = ngk;
                        pad_sphere_to_grid(ws, *ws.gvec, ngk, psi,
                                           work_buf.data());
                        ws.fft->c2c_bwd(work_buf.data());
                        for (int i = 0; i < nr; ++i) {
                            const Cx p(
                                static_cast<AccReal>(work_buf.data()[i][0]),
                                static_cast<AccReal>(work_buf.data()[i][1]));
                            const Cx vp =
                                dv_scr_r[static_cast<std::size_t>(i)] * p;
                            work_buf2.data()[i][0] =
                                static_cast<double>(vp.real());
                            work_buf2.data()[i][1] =
                                static_cast<double>(vp.imag());
                        }
                        ws.fft->c2c_fwd(work_buf2.data());
                        ws.gvec = &gv_kq[static_cast<std::size_t>(ik)];
                        ws.n_g = ngkq;
                        extract_sphere_from_grid(ws, *ws.gvec, ngkq,
                                                 work_buf2.data(), rhs.data());
                        const AccReal inv_nr = AccReal(1) / AccReal(nr);
                        for (int ig = 0; ig < ngkq; ++ig)
                            rhs[static_cast<std::size_t>(ig)] *= inv_nr;
                    }

                    // + ΔV_nl (cross k → k+q).
                    {
                        ws.gvec = &gv_k0[static_cast<std::size_t>(ik)];
                        ws.n_g = ngk;
                        ws.T_G.resize({static_cast<std::size_t>(ngk)});
                        for (int ig = 0; ig < ngk; ++ig)
                            ws.T_G[ig] =
                                gv_k0[static_cast<std::size_t>(ik)].ekin[ig];
                        ws.kb_proj = ks.kb_proj;
                        PWWorkspace ws_kq;
                        ws_kq.gvec = &gv_kq[static_cast<std::size_t>(ik)];
                        ws_kq.n_g = ngkq;
                        ws_kq.omega = omega;
                        ws_kq.nr = nr;
                        ws_kq.fft_nx = ws.fft_nx;
                        ws_kq.fft_ny = ws.fft_ny;
                        ws_kq.fft_nz = ws.fft_nz;
                        ws_kq.atom_pos_bohr = ws.atom_pos_bohr;
                        ws_kq.kb_proj = kqs.kb_proj;
                        add_dnl_psi_q(ws, ws_kq,
                                     static_cast<std::size_t>(atom_mu), axis_mu,
                                     psi, rhs.data());
                    }

                    for (int ig = 0; ig < ngkq; ++ig)
                        rhs[static_cast<std::size_t>(ig)] =
                            -rhs[static_cast<std::size_t>(ig)];

                    // Sternheimer at k+q.
                    ws.gvec = &gv_kq[static_cast<std::size_t>(ik)];
                    ws.n_g = ngkq;
                    ws.T_G.resize({static_cast<std::size_t>(ngkq)});
                    for (int ig = 0; ig < ngkq; ++ig)
                        ws.T_G[ig] =
                            gv_kq[static_cast<std::size_t>(ik)].ekin[ig];
                    ws.kb_proj = kqs.kb_proj;
                    for (int i = 0; i < nr; ++i)
                        ws.V_eff_r[i] = veff_frozen[static_cast<std::size_t>(i)];

                    auto& x = dpsi[static_cast<std::size_t>(ik) * n_bands
                                   + static_cast<std::size_t>(b)];
                    if (static_cast<int>(x.size()) != ngkq)
                        x.assign(static_cast<std::size_t>(ngkq), Cx(0, 0));
                    const int its = sternheimer_cg(
                        ws, kqs.wfc.data(), n_bands, ks.eps[b], rhs.data(),
                        x.data(), cg_tol, cg_max_iter, work_buf2.data());
                    if (its > cg_max_iter)
                        ATV_WARN("[pw][dfpt_q] CG not converged (mu={} ik={} "
                                 "band={})",
                                 mu, ik, b);

                    // Δn(r) += 2 w f ψ_k*(r) Δψ_{k+q}(r) / Ω
                    // (factor 2 matches Γ path / QE TR convention for the
                    // conjugate partner; lattice-periodic part only).
                    {
                        ws.gvec = &gv_k0[static_cast<std::size_t>(ik)];
                        ws.n_g = ngk;
                        pad_sphere_to_grid(ws, *ws.gvec, ngk, psi,
                                           work_buf.data());
                        ws.fft->c2c_bwd(work_buf.data());
                        ws.gvec = &gv_kq[static_cast<std::size_t>(ik)];
                        ws.n_g = ngkq;
                        pad_sphere_to_grid(ws, *ws.gvec, ngkq, x.data(),
                                           work_buf2.data());
                        ws.fft->c2c_bwd(work_buf2.data());
                        const AccReal pref =
                            w * occ_f * AccReal(2) / omega;
                        for (int i = 0; i < nr; ++i) {
                            const Cx pk(
                                static_cast<AccReal>(work_buf.data()[i][0]),
                                static_cast<AccReal>(work_buf.data()[i][1]));
                            const Cx dk(
                                static_cast<AccReal>(work_buf2.data()[i][0]),
                                static_cast<AccReal>(work_buf2.data()[i][1]));
                            dn_new[static_cast<std::size_t>(i)] +=
                                pref * std::conj(pk) * dk;
                        }
                    }
                }
            }
            ws.gvec = &prob.gv0;
            ws.n_g = ng0;

            AccReal diff2 = 0;
            for (int i = 0; i < nr; ++i) {
                const Cx d = dn_new[static_cast<std::size_t>(i)]
                             - dn_mixed[static_cast<std::size_t>(i)];
                diff2 += (d * std::conj(d)).real();
            }
            const AccReal delta = std::sqrt(diff2 * dV);
            const bool this_converged = (delta < response_tol);

            for (int i = 0; i < nr; ++i)
                dn_mixed[static_cast<std::size_t>(i)] +=
                    mix_beta
                    * (dn_new[static_cast<std::size_t>(i)]
                       - dn_mixed[static_cast<std::size_t>(i)]);

            complex_grid_to_sphere(ws, gv0, dn_mixed.data(), dn_g.data(),
                                   work_buf.data());
            for (int ig = 0; ig < ng0; ++ig) {
                const AccReal gqx = gv0.g[ig][0] + q_cart[0];
                const AccReal gqy = gv0.g[ig][1] + q_cart[1];
                const AccReal gqz = gv0.g[ig][2] + q_cart[2];
                const AccReal gq2 = gqx * gqx + gqy * gqy + gqz * gqz;
                if (gq2 < AccReal(1e-16))
                    dn_g[static_cast<std::size_t>(ig)] = Cx(0, 0);
                else
                    dn_g[static_cast<std::size_t>(ig)] *=
                        AccReal(4) * AccReal(units::pi) / gq2;
            }
            sphere_to_complex_grid(ws, gv0, dn_g.data(), dv_hxc_r.data(),
                                   work_buf.data());
            {
                std::vector<Cx> tmp(static_cast<std::size_t>(nr));
                for (int i = 0; i < nr; ++i)
                    tmp[static_cast<std::size_t>(i)] =
                        fxc_r[static_cast<std::size_t>(i)]
                        * dn_mixed[static_cast<std::size_t>(i)];
                complex_grid_to_sphere(ws, gv0, tmp.data(), dn_g.data(),
                                       work_buf.data());
                std::vector<Cx> tmp2(static_cast<std::size_t>(nr));
                sphere_to_complex_grid(ws, gv0, dn_g.data(), tmp2.data(),
                                       work_buf.data());
                for (int i = 0; i < nr; ++i)
                    dv_hxc_r[static_cast<std::size_t>(i)] +=
                        tmp2[static_cast<std::size_t>(i)];
            }

            if (this_converged) {
                converged = true;
                ATV_DEBUG("[pw][dfpt_q] response mu={} converged iter={} "
                          "delta={:.3e}",
                          mu, rit + 1, delta);
                break;
            }
        }
        if (!converged) {
            all_responses_converged = false;
            ATV_WARN("[pw][dfpt_q] response mu={} not converged after {} iters",
                     mu, response_max_iter);
        }

        // Local response contribution C[ν][μ] += ∫ Δn^{(μ)} (ΔV_loc^{(ν)})*
        for (int nu = 0; nu < dim; ++nu) {
            Cx s(0, 0);
            const auto& dvl = dvloc_r[static_cast<std::size_t>(nu)];
            for (int i = 0; i < nr; ++i)
                s += dn_mixed[static_cast<std::size_t>(i)]
                     * std::conj(dvl[static_cast<std::size_t>(i)]);
            Cat(nu, mu) += s * dV;
        }

        // Nonlocal response.
        {
            std::vector<KbMoments> mo_psi, mo_dpsi;
            for (int ik = 0; ik < n_kpts; ++ik) {
                auto& ks = kstates[static_cast<std::size_t>(ik)];
                auto& kqs = kqstates[static_cast<std::size_t>(ik)];
                const AccReal w = k_weights[static_cast<std::size_t>(ik)];
                const int ngk = ks.ng;
                for (int b = 0; b < n_bands; ++b) {
                    const Cx* psi =
                        ks.wfc.data() + static_cast<std::ptrdiff_t>(b) * ngk;
                    const auto& x = dpsi[static_cast<std::size_t>(ik) * n_bands
                                         + static_cast<std::size_t>(b)];
                    ws.gvec = &gv_k0[static_cast<std::size_t>(ik)];
                    ws.n_g = ngk;
                    ws.kb_proj = ks.kb_proj;
                    PWWorkspace ws_kq;
                    ws_kq.gvec = &gv_kq[static_cast<std::size_t>(ik)];
                    ws_kq.n_g = kqs.ng;
                    ws_kq.omega = omega;
                    ws_kq.atom_pos_bohr = ws.atom_pos_bohr;
                    ws_kq.kb_proj = kqs.kb_proj;
                    for (std::size_t ja = 0;
                         ja < std::min(ks.kb_proj.size(), kqs.kb_proj.size());
                         ++ja) {
                        const std::size_t nproj = std::min(
                            ks.kb_proj[ja].size(), kqs.kb_proj[ja].size());
                        for (std::size_t ip = 0; ip < nproj; ++ip) {
                            const AccReal Dc = ks.kb_proj[ja][ip].d;
                            if (std::abs(Dc) < AccReal(1e-16))
                                continue;
                            kb_moments(ws, ja, ip, psi, false, mo_psi);
                            kb_moments(ws_kq, ja, ip, x.data(), false, mo_dpsi);
                            for (std::size_t m = 0; m < mo_psi.size(); ++m) {
                                const Cx cP = mo_psi[m].c;
                                const Cx cD = mo_dpsi[m].c;
                                for (int bb = 0; bb < 3; ++bb) {
                                    const Cx el =
                                        Cx(0, 1) * Dc
                                        * (std::conj(cD) * mo_psi[m].s[bb]
                                           - std::conj(mo_dpsi[m].s[bb]) * cP);
                                    Cat(3 * static_cast<int>(ja) + bb, mu) +=
                                        w * occ_f * AccReal(2) * el;
                                }
                            }
                        }
                    }
                }
            }
            ws.gvec = &prob.gv0;
            ws.n_g = ng0;
        }
    }
    out.response_converged = all_responses_converged;

    // Hermitian projection.
    AccReal max_asym = 0;
    for (int i = 0; i < dim; ++i)
        for (int j = i; j < dim; ++j) {
            const Cx cij = Cat(i, j);
            const Cx cji = Cat(j, i);
            const Cx avg = AccReal(0.5) * (cij + std::conj(cji));
            const AccReal d = std::abs(cij - std::conj(cji)) * AccReal(0.5);
            if (d > max_asym)
                max_asym = d;
            Cat(i, j) = avg;
            Cat(j, i) = std::conj(avg);
        }
    out.max_asymmetry = max_asym;
    out.force_constants = std::move(C);
    return out;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
