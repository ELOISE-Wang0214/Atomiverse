/// @file cpu/pw/internal/pw_subspace.cpp
/// @brief Davidson subspace utilities (moved verbatim from pw_davidson.cpp).

#include "atomiverse/cpu/pw/internal/pw_subspace.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>

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

} // namespace

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

void jacobi_eigh(int n, Real* A, Real* eig, Real* vecs, int max_sweep) {
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

int init_subspace(const PWWorkspace& ws, int n_bands,
                  std::vector<std::complex<Real>>& sub_vectors) {
    int ng = ws.n_g;
    int nb = n_bands;

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
    sub_vectors.assign(static_cast<size_t>(n_sub) * ng,
                       std::complex<Real>(0, 0));

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
                        std::complex<Real> S(std::cos(phase),
                                             -std::sin(phase));
                        sub_vectors[j * ng + ig] = il * ylm[m] * chi * S;
                    }
                    // QE startingwfc='atomic+random': add random perturbation
                    for (int ig = 0; ig < ng; ++ig) {
                        Real rnd =
                            Real(std::rand()) / Real(RAND_MAX) - Real(0.5);
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
    return n_sub;
}

} // namespace internal
} // namespace pw
} // namespace cpu
} // namespace ATV
