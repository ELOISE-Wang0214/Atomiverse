/// @file src/cpu/gpw/internal/mp2_kernels.cpp
/// @brief GPW MP2 kernels — MO-on-grid, (ia|jb) via per-pair Poisson or
///        RI-GPW G-space Coulomb intermediates, and the closed-shell
///        doubles sum. See mp2_kernels.h.

#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"

// solve_poisson_q is used by multi-k RPA finite-q pair densities.

namespace ATV {
namespace cpu {
namespace gpw {
namespace internal {

std::vector<double> mo_on_grid(const CompactBasisGrid& cbg,
                               const linalg::RealMatrix& C, int p0, int np) {
    const int n_mo = static_cast<int>(C.shape(1));
    if (p0 < 0 || np < 0 || p0 + np > n_mo)
        throw std::invalid_argument(
            "gpw mo_on_grid: MO range [" + std::to_string(p0) + ", "
            + std::to_string(p0 + np) + ") outside 0.." + std::to_string(n_mo));
    const long npts = cbg.npts;
    std::vector<double> psi(
        static_cast<std::size_t>(np) * static_cast<std::size_t>(npts), 0.0);
    const std::size_t na = cbg.gidx.size();
    for (std::size_t k = 0; k < na; ++k) {
        const long g = cbg.gidx[k];
        for (std::size_t nz = cbg.off[k]; nz < cbg.off[k + 1]; ++nz) {
            const int u = cbg.fu[nz];
            const double val = cbg.fv[nz];
            for (int p = 0; p < np; ++p)
                psi[static_cast<std::size_t>(p) * static_cast<std::size_t>(npts)
                    + static_cast<std::size_t>(g)] += C(u, p0 + p) * val;
        }
    }
    return psi;
}

std::vector<double> mo_eri_ovov(const UniformGrid& grid,
                                const std::vector<double>& psi_occ, int n_occ,
                                const std::vector<double>& psi_virt, int n_virt,
                                bool nonperiodic) {
    const long npts = grid.npts();
    const double w = grid.weight();
    if (psi_occ.size()
            != static_cast<std::size_t>(n_occ) * static_cast<std::size_t>(npts)
        || psi_virt.size()
               != static_cast<std::size_t>(n_virt)
                      * static_cast<std::size_t>(npts))
        throw std::invalid_argument(
            "gpw mo_eri_ovov: MO grid arrays do not match n_occ/n_virt/npts");
    const int nov = n_occ * n_virt;
    std::vector<double> v(static_cast<std::size_t>(nov)
                          * static_cast<std::size_t>(nov));
    std::vector<double> rho(static_cast<std::size_t>(npts));
    for (int ia = 0; ia < nov; ++ia) {
        const double* pi = psi_occ.data()
                           + static_cast<std::size_t>(ia / n_virt)
                                 * static_cast<std::size_t>(npts);
        const double* pa = psi_virt.data()
                           + static_cast<std::size_t>(ia % n_virt)
                                 * static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            rho[static_cast<std::size_t>(g)] = pi[g] * pa[g];
        const PoissonResult pr = nonperiodic ? solve_poisson_mt(grid, rho)
                                             : solve_poisson(grid, rho);
        // (ia|jb) = (jb|ia): fill the upper triangle jb >= ia and mirror.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int jb = ia; jb < nov; ++jb) {
            const double* pj = psi_occ.data()
                               + static_cast<std::size_t>(jb / n_virt)
                                     * static_cast<std::size_t>(npts);
            const double* pb = psi_virt.data()
                               + static_cast<std::size_t>(jb % n_virt)
                                     * static_cast<std::size_t>(npts);
            double s = 0.0;
            for (long g = 0; g < npts; ++g)
                s += pr.v_r[static_cast<std::size_t>(g)] * pj[g] * pb[g];
            const double val = w * s;
            v[static_cast<std::size_t>(ia) * static_cast<std::size_t>(nov)
              + static_cast<std::size_t>(jb)] = val;
            v[static_cast<std::size_t>(jb) * static_cast<std::size_t>(nov)
              + static_cast<std::size_t>(ia)] = val;
        }
    }
    return v;
}

std::vector<double> mo_eri_ovov_ri(const UniformGrid& grid,
                                   const std::vector<double>& psi_occ,
                                   int n_occ,
                                   const std::vector<double>& psi_virt,
                                   int n_virt) {
    // Product RI-GPW path (Del Ben 2012 / CP2K MP2_GPW): form the Coulomb
    // potential of every pair density once (G-space 4π/G² intermediate via
    // the shared FFT Poisson), store the real-space V_ia intermediates
    // (three-index product structure, memory O(n_ov · N_g)), then contract
    // (ia|jb) = Σ_g w V_ia(g) ρ_jb(g). Numerically identical to
    // mo_eri_ovov(..., nonperiodic=false) — the Poisson kernel is the PW
    // resolution of the identity for the pair density on this grid.
    const long npts = grid.npts();
    const double w = grid.weight();
    if (psi_occ.size()
            != static_cast<std::size_t>(n_occ) * static_cast<std::size_t>(npts)
        || psi_virt.size()
               != static_cast<std::size_t>(n_virt)
                      * static_cast<std::size_t>(npts))
        throw std::invalid_argument(
            "gpw mo_eri_ovov_ri: MO grid arrays do not match n_occ/n_virt/npts");
    if (n_occ < 0 || n_virt < 0)
        throw std::invalid_argument("gpw mo_eri_ovov_ri: negative orbital count");

    const int nov = n_occ * n_virt;
    std::vector<double> V(static_cast<std::size_t>(nov)
                          * static_cast<std::size_t>(npts));
    std::vector<double> rho(static_cast<std::size_t>(npts));
    for (int ia = 0; ia < nov; ++ia) {
        const double* pi = psi_occ.data()
                           + static_cast<std::size_t>(ia / n_virt)
                                 * static_cast<std::size_t>(npts);
        const double* pa = psi_virt.data()
                           + static_cast<std::size_t>(ia % n_virt)
                                 * static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            rho[static_cast<std::size_t>(g)] = pi[g] * pa[g];
        // FFT → 4π/G² → IFFT (G=0 dropped); owns the RI-GPW intermediate.
        const PoissonResult pr = solve_poisson(grid, rho);
        double* Via = V.data() + static_cast<std::size_t>(ia)
                                      * static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            Via[g] = pr.v_r[static_cast<std::size_t>(g)];
    }

    std::vector<double> v(static_cast<std::size_t>(nov)
                          * static_cast<std::size_t>(nov));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int ia = 0; ia < nov; ++ia) {
        const double* Via = V.data() + static_cast<std::size_t>(ia)
                                           * static_cast<std::size_t>(npts);
        for (int jb = ia; jb < nov; ++jb) {
            const double* pj = psi_occ.data()
                               + static_cast<std::size_t>(jb / n_virt)
                                     * static_cast<std::size_t>(npts);
            const double* pb = psi_virt.data()
                               + static_cast<std::size_t>(jb % n_virt)
                                     * static_cast<std::size_t>(npts);
            double s = 0.0;
            for (long g = 0; g < npts; ++g)
                s += Via[g] * pj[g] * pb[g];
            const double val = w * s;
            v[static_cast<std::size_t>(ia) * static_cast<std::size_t>(nov)
              + static_cast<std::size_t>(jb)] = val;
            v[static_cast<std::size_t>(jb) * static_cast<std::size_t>(nov)
              + static_cast<std::size_t>(ia)] = val;
        }
    }
    return v;
}

void mp2_rhf_energy(const std::vector<double>& v_ovov,
                    const std::vector<double>& eps_occ,
                    const std::vector<double>& eps_virt, double& e_os,
                    double& e_ss) {
    const int no = static_cast<int>(eps_occ.size());
    const int nv = static_cast<int>(eps_virt.size());
    const std::size_t nov =
        static_cast<std::size_t>(no) * static_cast<std::size_t>(nv);
    if (v_ovov.size() != nov * nov)
        throw std::invalid_argument(
            "gpw mp2_rhf_energy: tensor size does not match eps_occ/eps_virt");
    e_os = 0.0;
    e_ss = 0.0;
    for (int i = 0; i < no; ++i)
        for (int j = 0; j < no; ++j)
            for (int a = 0; a < nv; ++a)
                for (int b = 0; b < nv; ++b) {
                    const double v_iajb =
                        v_ovov[(static_cast<std::size_t>(i) * nv + a) * nov
                               + static_cast<std::size_t>(j) * nv + b];
                    const double v_ibja =
                        v_ovov[(static_cast<std::size_t>(i) * nv + b) * nov
                               + static_cast<std::size_t>(j) * nv + a];
                    const double denom =
                        eps_occ[static_cast<std::size_t>(i)]
                        + eps_occ[static_cast<std::size_t>(j)]
                        - eps_virt[static_cast<std::size_t>(a)]
                        - eps_virt[static_cast<std::size_t>(b)];
                    e_os += v_iajb * v_iajb / denom;
                    e_ss += v_iajb * (v_iajb - v_ibja) / denom;
                }
}

std::vector<Cplx> mo_on_grid_k(const KpointBasisGrid& kbg, int ik,
                               const linalg::CplxMatrix& C, int p0, int np) {
    const int n_mo = static_cast<int>(C.shape(1));
    const int nb = static_cast<int>(C.shape(0));
    if (p0 < 0 || np < 0 || p0 + np > n_mo)
        throw std::invalid_argument(
            "gpw mo_on_grid_k: MO range [" + std::to_string(p0) + ", "
            + std::to_string(p0 + np) + ") outside 0.." + std::to_string(n_mo));
    if (ik < 0 || ik >= static_cast<int>(kbg.chi.size()))
        throw std::invalid_argument("gpw mo_on_grid_k: k index out of range");
    if (nb != kbg.nb)
        throw std::invalid_argument(
            "gpw mo_on_grid_k: C.nb does not match KpointBasisGrid.nb");
    const long npts = kbg.npts;
    std::vector<Cplx> psi(
        static_cast<std::size_t>(np) * static_cast<std::size_t>(npts),
        Cplx{0.0, 0.0});
    // KpointBasisGrid::chi[k] is dense over all npts: chi[g * nb + u].
    const auto& chi_k = kbg.chi[static_cast<std::size_t>(ik)];
    if (chi_k.size()
        != static_cast<std::size_t>(npts) * static_cast<std::size_t>(nb))
        throw std::logic_error(
            "gpw mo_on_grid_k: chi_k size is not npts * nb (dense layout)");
    for (long g = 0; g < npts; ++g) {
        const Cplx* chi =
            chi_k.data() + static_cast<std::size_t>(g) * static_cast<std::size_t>(nb);
        for (int p = 0; p < np; ++p) {
            Cplx s{0.0, 0.0};
            for (int u = 0; u < nb; ++u)
                s += C(u, p0 + p) * chi[u];
            psi[static_cast<std::size_t>(p) * static_cast<std::size_t>(npts)
                + static_cast<std::size_t>(g)] = s;
        }
    }
    return psi;
}

namespace {

// Wrap fractional κ into (−1/2, 1/2].
double wrap_frac(double x) {
    x -= std::floor(x + 0.5);
    if (x <= -0.5)
        x += 1.0;
    if (x > 0.5)
        x -= 1.0;
    return x;
}

int find_k_index(const std::vector<std::array<double, 3>>& k_frac, double fx,
                 double fy, double fz, double tol = 1e-8) {
    const double tx = wrap_frac(fx), ty = wrap_frac(fy), tz = wrap_frac(fz);
    for (std::size_t i = 0; i < k_frac.size(); ++i) {
        const double dx = wrap_frac(k_frac[i][0] - tx);
        const double dy = wrap_frac(k_frac[i][1] - ty);
        const double dz = wrap_frac(k_frac[i][2] - tz);
        if (std::abs(dx) < tol && std::abs(dy) < tol && std::abs(dz) < tol)
            return static_cast<int>(i);
    }
    return -1;
}

// Coulomb potential of a complex pair density via linear real Poisson.
// Optional crystal-momentum transfer q_frac (nullptr ⇒ q = 0, G=0 dropped).
void poisson_complex(const UniformGrid& grid, const Cplx* rho,
                     std::vector<Cplx>& V, const double* q_frac = nullptr) {
    const long npts = grid.npts();
    std::vector<double> re(static_cast<std::size_t>(npts));
    std::vector<double> im(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g) {
        re[static_cast<std::size_t>(g)] = rho[g].real();
        im[static_cast<std::size_t>(g)] = rho[g].imag();
    }
    const double q0[3] = {0.0, 0.0, 0.0};
    const double* q = q_frac ? q_frac : q0;
    const PoissonResult vr = solve_poisson_q(grid, re, q);
    const PoissonResult vi = solve_poisson_q(grid, im, q);
    V.resize(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g)
        V[static_cast<std::size_t>(g)] =
            Cplx{vr.v_r[static_cast<std::size_t>(g)],
                 vi.v_r[static_cast<std::size_t>(g)]};
}

} // namespace

void mp2_rhf_energy_kpoints(
    const UniformGrid& grid,
    const std::vector<std::vector<Cplx>>& psi_occ_k,
    const std::vector<std::vector<Cplx>>& psi_virt_k,
    const std::vector<std::vector<double>>& eps_occ_k,
    const std::vector<std::vector<double>>& eps_virt_k,
    const std::vector<std::array<double, 3>>& k_frac,
    const std::vector<double>& k_weights, double& e_os, double& e_ss) {
    const int nk = static_cast<int>(k_frac.size());
    if (nk < 1)
        throw std::invalid_argument("gpw mp2_rhf_energy_kpoints: empty k mesh");
    if (static_cast<int>(psi_occ_k.size()) != nk
        || static_cast<int>(psi_virt_k.size()) != nk
        || static_cast<int>(eps_occ_k.size()) != nk
        || static_cast<int>(eps_virt_k.size()) != nk
        || static_cast<int>(k_weights.size()) != nk)
        throw std::invalid_argument(
            "gpw mp2_rhf_energy_kpoints: per-k array length mismatch");

    const int n_occ = static_cast<int>(eps_occ_k[0].size());
    const int n_virt = static_cast<int>(eps_virt_k[0].size());
    if (n_occ < 1 || n_virt < 1)
        throw std::invalid_argument(
            "gpw mp2_rhf_energy_kpoints: need n_occ>=1 and n_virt>=1");
    const long npts = grid.npts();
    const double w = grid.weight();
    const int nov = n_occ * n_virt;

    for (int ik = 0; ik < nk; ++ik) {
        if (static_cast<int>(eps_occ_k[static_cast<std::size_t>(ik)].size())
                != n_occ
            || static_cast<int>(
                   eps_virt_k[static_cast<std::size_t>(ik)].size())
                   != n_virt)
            throw std::invalid_argument(
                "gpw mp2_rhf_energy_kpoints: inconsistent n_occ/n_virt across k");
        const std::size_t need_o =
            static_cast<std::size_t>(n_occ) * static_cast<std::size_t>(npts);
        const std::size_t need_v =
            static_cast<std::size_t>(n_virt) * static_cast<std::size_t>(npts);
        if (psi_occ_k[static_cast<std::size_t>(ik)].size() != need_o
            || psi_virt_k[static_cast<std::size_t>(ik)].size() != need_v)
            throw std::invalid_argument(
                "gpw mp2_rhf_energy_kpoints: MO grid size mismatch");
    }

    // Equal-weight factor: E = (1/Nk) Σ_{ki,kj,ka} ... when w_k = 1/Nk.
    // General equal-ish: use Nk * w_ki * w_kj * Nk * w_ka… → for equal w:
    // weight_triple = Nk * w_i * w_j  (kb fixed; one free sum over ka with
    // w absorbed as 1 when all w=1/Nk ⇒ weight = 1/Nk).
    // Use: weight = nk * w[ki] * w[kj] * nk * w[ka] / (nk * w[kb]) simplified
    // for equal w to 1/nk. For robustness with equal weights only:
    double wsum = 0.0;
    for (double wk : k_weights)
        wsum += wk;
    if (std::abs(wsum - 1.0) > 1e-8)
        throw std::invalid_argument(
            "gpw mp2_rhf_energy_kpoints: k_weights must sum to 1");

    // Precompute V_ia for every (ki, ka, i, a): complex Coulomb potential of
    // ρ_ia = conj(ψ_i^{ki}) ψ_a^{ka}. Memory O(nk² · nov · npts).
    // Layout: Via_cache[ki * nk + ka][ia * npts + g]
    std::vector<std::vector<Cplx>> Via_cache(
        static_cast<std::size_t>(nk) * static_cast<std::size_t>(nk));
    std::vector<Cplx> rho(static_cast<std::size_t>(npts));
    for (int ki = 0; ki < nk; ++ki) {
        for (int ka = 0; ka < nk; ++ka) {
            auto& block =
                Via_cache[static_cast<std::size_t>(ki) * nk
                          + static_cast<std::size_t>(ka)];
            block.resize(static_cast<std::size_t>(nov)
                         * static_cast<std::size_t>(npts));
            const auto& po = psi_occ_k[static_cast<std::size_t>(ki)];
            const auto& pv = psi_virt_k[static_cast<std::size_t>(ka)];
            for (int i = 0; i < n_occ; ++i) {
                for (int a = 0; a < n_virt; ++a) {
                    const int ia = i * n_virt + a;
                    const Cplx* pi =
                        po.data()
                        + static_cast<std::size_t>(i)
                              * static_cast<std::size_t>(npts);
                    const Cplx* pa =
                        pv.data()
                        + static_cast<std::size_t>(a)
                              * static_cast<std::size_t>(npts);
                    for (long g = 0; g < npts; ++g)
                        rho[static_cast<std::size_t>(g)] =
                            std::conj(pi[g]) * pa[g];
                    std::vector<Cplx> V;
                    poisson_complex(grid, rho.data(), V);
                    Cplx* dest =
                        block.data()
                        + static_cast<std::size_t>(ia)
                              * static_cast<std::size_t>(npts);
                    for (long g = 0; g < npts; ++g)
                        dest[g] = V[static_cast<std::size_t>(g)];
                }
            }
        }
    }

    auto contract = [&](int ki, int ka, int i, int a, int kj, int kb, int j,
                        int b) -> Cplx {
        // (ia|jb) = Σ_g w V_ia(g) ρ_jb(g), ρ_jb = conj(ψ_j) ψ_b
        const auto& Vblock =
            Via_cache[static_cast<std::size_t>(ki) * nk
                      + static_cast<std::size_t>(ka)];
        const int ia = i * n_virt + a;
        const Cplx* Via =
            Vblock.data()
            + static_cast<std::size_t>(ia) * static_cast<std::size_t>(npts);
        const Cplx* pj =
            psi_occ_k[static_cast<std::size_t>(kj)].data()
            + static_cast<std::size_t>(j) * static_cast<std::size_t>(npts);
        const Cplx* pb =
            psi_virt_k[static_cast<std::size_t>(kb)].data()
            + static_cast<std::size_t>(b) * static_cast<std::size_t>(npts);
        Cplx s{0.0, 0.0};
        for (long g = 0; g < npts; ++g)
            s += Via[g] * (std::conj(pj[g]) * pb[g]);
        return w * s;
    };

    e_os = 0.0;
    e_ss = 0.0;
    // weight_triple = 1/Nk for equal weights (w_k = 1/Nk).
    // Generalization used: Nk * w[ki] * w[kj]  (matches 1/Nk when equal;
    // free sum over ka with kb determined — each (ki,kj,ka) counted once).
    for (int ki = 0; ki < nk; ++ki) {
        for (int kj = 0; kj < nk; ++kj) {
            const double w_ij =
                static_cast<double>(nk) * k_weights[static_cast<std::size_t>(ki)]
                * k_weights[static_cast<std::size_t>(kj)];
            for (int ka = 0; ka < nk; ++ka) {
                const double fbx = k_frac[static_cast<std::size_t>(ki)][0]
                                   + k_frac[static_cast<std::size_t>(kj)][0]
                                   - k_frac[static_cast<std::size_t>(ka)][0];
                const double fby = k_frac[static_cast<std::size_t>(ki)][1]
                                   + k_frac[static_cast<std::size_t>(kj)][1]
                                   - k_frac[static_cast<std::size_t>(ka)][1];
                const double fbz = k_frac[static_cast<std::size_t>(ki)][2]
                                   + k_frac[static_cast<std::size_t>(kj)][2]
                                   - k_frac[static_cast<std::size_t>(ka)][2];
                const int kb = find_k_index(k_frac, fbx, fby, fbz);
                if (kb < 0)
                    throw std::runtime_error(
                        "gpw mp2_rhf_energy_kpoints: momentum-conserving kb "
                        "not in the sampled mesh — use a closed (uniform) "
                        "k-mesh");
                for (int i = 0; i < n_occ; ++i)
                    for (int j = 0; j < n_occ; ++j)
                        for (int a = 0; a < n_virt; ++a)
                            for (int b = 0; b < n_virt; ++b) {
                                const Cplx v_iajb =
                                    contract(ki, ka, i, a, kj, kb, j, b);
                                const Cplx v_ibja =
                                    contract(ki, kb, i, b, kj, ka, j, a);
                                const double denom =
                                    eps_occ_k[static_cast<std::size_t>(ki)]
                                             [static_cast<std::size_t>(i)]
                                    + eps_occ_k[static_cast<std::size_t>(kj)]
                                               [static_cast<std::size_t>(j)]
                                    - eps_virt_k[static_cast<std::size_t>(ka)]
                                                [static_cast<std::size_t>(a)]
                                    - eps_virt_k[static_cast<std::size_t>(kb)]
                                                [static_cast<std::size_t>(b)];
                                if (std::abs(denom) < 1e-12)
                                    continue;
                                const double inv = w_ij / denom;
                                e_os += std::norm(v_iajb) * inv;
                                e_ss += (v_iajb * std::conj(v_iajb - v_ibja))
                                            .real()
                                        * inv;
                            }
            }
        }
    }
}

void direct_rpa_rhf(const std::vector<double>& v_ovov,
                    const std::vector<double>& eps_occ,
                    const std::vector<double>& eps_virt,
                    DirectRpaAmplitudes& out) {
    using linalg::RealMatrix;
    using linalg::sym_power;

    const int no = static_cast<int>(eps_occ.size());
    const int nv = static_cast<int>(eps_virt.size());
    if (no < 1 || nv < 1)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf: need n_occ >= 1 and n_virt >= 1");
    const std::size_t nov =
        static_cast<std::size_t>(no) * static_cast<std::size_t>(nv);
    if (v_ovov.size() != nov * nov)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf: tensor size does not match eps_occ/eps_virt");

    // Singlet dRPA (α_hf = 0): A = Δε δ + 2 V, B = 2 V, index p = i*nv + a.
    RealMatrix A =
        xt::zeros<double>({static_cast<int>(nov), static_cast<int>(nov)});
    RealMatrix B =
        xt::zeros<double>({static_cast<int>(nov), static_cast<int>(nov)});
    out.tr_A = 0.0;
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a) {
            const int p = i * nv + a;
            const double de = eps_virt[static_cast<std::size_t>(a)]
                              - eps_occ[static_cast<std::size_t>(i)];
            for (int j = 0; j < no; ++j)
                for (int b = 0; b < nv; ++b) {
                    const int q = j * nv + b;
                    const double v =
                        v_ovov[static_cast<std::size_t>(p) * nov
                               + static_cast<std::size_t>(q)];
                    const double two_v = 2.0 * v;
                    double Aval = two_v;
                    if (p == q)
                        Aval += de;
                    A(p, q) = Aval;
                    B(p, q) = two_v;
                }
            out.tr_A += A(p, p);
        }

    constexpr double kAmBFloor = 1e-10;
    constexpr double kOmega2Tol = 1e-10;
    const RealMatrix Mp = A + B;
    const RealMatrix Mm = A - B;
    {
        auto [s, U] = xt::linalg::eigh(Mm);
        (void)U;
        if (s(0) <= kAmBFloor)
            throw std::runtime_error(
                "gpw direct_rpa_rhf: (A - B) is not positive definite "
                "(lowest eigenvalue = "
                + std::to_string(s(0))
                + " Ha) — dRPA reduction undefined (no fallback)");
    }
    const RealMatrix T = sym_power(Mm, 0.5, kAmBFloor, "gpw direct_rpa_rhf");
    const RealMatrix Ti =
        sym_power(Mm, -0.5, kAmBFloor, "gpw direct_rpa_rhf");
    RealMatrix M = xt::linalg::dot(T, RealMatrix(xt::linalg::dot(Mp, T)));
    for (std::size_t r = 0; r < nov; ++r)
        for (std::size_t c = r + 1; c < nov; ++c) {
            const double v = 0.5 * (M(r, c) + M(c, r));
            M(r, c) = v;
            M(c, r) = v;
        }
    auto [w2, Z] = xt::linalg::eigh(M);

    out.omega.assign(nov, 0.0);
    out.X.assign(nov * nov, 0.0);
    out.Y.assign(nov * nov, 0.0);
    double sum_om = 0.0;
    for (std::size_t s = 0; s < nov; ++s) {
        if (w2(s) < -kOmega2Tol)
            throw std::runtime_error(
                "gpw direct_rpa_rhf: dRPA omega^2 = " + std::to_string(w2(s))
                + " Ha^2 (root " + std::to_string(s)
                + ") is negative — imaginary excitation (no fallback)");
        if (w2(s) <= kOmega2Tol)
            throw std::runtime_error(
                "gpw direct_rpa_rhf: dRPA omega^2 = " + std::to_string(w2(s))
                + " Ha^2 (root " + std::to_string(s)
                + ") is numerically zero (no fallback)");
        const double w = std::sqrt(w2(s));
        out.omega[s] = w;
        sum_om += w;
        // (X+Y) = ω^{-1/2} T z, (X−Y) = ω^{1/2} T^{-1} z
        for (std::size_t r = 0; r < nov; ++r) {
            double xpy = 0.0, xmy = 0.0;
            for (std::size_t c = 0; c < nov; ++c) {
                xpy += T(static_cast<int>(r), static_cast<int>(c))
                       * Z(static_cast<int>(c), static_cast<int>(s));
                xmy += Ti(static_cast<int>(r), static_cast<int>(c))
                       * Z(static_cast<int>(c), static_cast<int>(s));
            }
            xpy /= std::sqrt(w);
            xmy *= std::sqrt(w);
            out.X[s * nov + r] = 0.5 * (xpy + xmy);
            out.Y[s * nov + r] = 0.5 * (xpy - xmy);
        }
    }
    out.e_corr = 0.5 * (sum_om - out.tr_A);
}

void direct_rpa_rhf_energy(const std::vector<double>& v_ovov,
                           const std::vector<double>& eps_occ,
                           const std::vector<double>& eps_virt, double& e_corr,
                           double& tr_A, std::vector<double>& omega) {
    DirectRpaAmplitudes full;
    direct_rpa_rhf(v_ovov, eps_occ, eps_virt, full);
    e_corr = full.e_corr;
    tr_A = full.tr_A;
    omega = std::move(full.omega);
}

std::vector<double> mo_pair_potentials_ov(const UniformGrid& grid,
                                          const std::vector<double>& psi_occ,
                                          int n_occ,
                                          const std::vector<double>& psi_virt,
                                          int n_virt) {
    const long npts = grid.npts();
    if (psi_occ.size()
            != static_cast<std::size_t>(n_occ) * static_cast<std::size_t>(npts)
        || psi_virt.size()
               != static_cast<std::size_t>(n_virt)
                      * static_cast<std::size_t>(npts))
        throw std::invalid_argument(
            "gpw mo_pair_potentials_ov: MO grid arrays do not match "
            "n_occ/n_virt/npts");
    if (n_occ < 0 || n_virt < 0)
        throw std::invalid_argument(
            "gpw mo_pair_potentials_ov: negative orbital count");

    const int nov = n_occ * n_virt;
    std::vector<double> V(static_cast<std::size_t>(nov)
                          * static_cast<std::size_t>(npts));
    std::vector<double> rho(static_cast<std::size_t>(npts));
    for (int ia = 0; ia < nov; ++ia) {
        const double* pi = psi_occ.data()
                           + static_cast<std::size_t>(ia / n_virt)
                                 * static_cast<std::size_t>(npts);
        const double* pa = psi_virt.data()
                           + static_cast<std::size_t>(ia % n_virt)
                                 * static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            rho[static_cast<std::size_t>(g)] = pi[g] * pa[g];
        const PoissonResult pr = solve_poisson(grid, rho);
        double* Via = V.data() + static_cast<std::size_t>(ia)
                                      * static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            Via[g] = pr.v_r[static_cast<std::size_t>(g)];
    }
    return V;
}

std::vector<double> exchange_diag_mt(const UniformGrid& grid,
                                     const std::vector<double>& psi_mo, int n_mo,
                                     int n_occ) {
    const long npts = grid.npts();
    const double w = grid.weight();
    if (n_mo < 1 || n_occ < 1 || n_occ > n_mo)
        throw std::invalid_argument(
            "gpw exchange_diag_mt: bad n_mo/n_occ");
    if (psi_mo.size()
        != static_cast<std::size_t>(n_mo) * static_cast<std::size_t>(npts))
        throw std::invalid_argument(
            "gpw exchange_diag_mt: MO grid size mismatch");

    std::vector<double> sx(static_cast<std::size_t>(n_mo), 0.0);
    std::vector<double> rho(static_cast<std::size_t>(npts));
    for (int p = 0; p < n_mo; ++p) {
        const double* pp = psi_mo.data()
                           + static_cast<std::size_t>(p)
                                 * static_cast<std::size_t>(npts);
        double sum = 0.0;
        for (int i = 0; i < n_occ; ++i) {
            const double* pi = psi_mo.data()
                               + static_cast<std::size_t>(i)
                                     * static_cast<std::size_t>(npts);
            for (long g = 0; g < npts; ++g)
                rho[static_cast<std::size_t>(g)] = pp[g] * pi[g];
            // MT truncated Coulomb — G=0-safe for charged (ii|ii) pairs.
            const PoissonResult pr = solve_poisson_mt(grid, rho);
            double s = 0.0;
            for (long g = 0; g < npts; ++g)
                s += pr.v_r[static_cast<std::size_t>(g)]
                     * rho[static_cast<std::size_t>(g)];
            sum -= w * s; // Σ_x,p = −Σ_i (pi|pi)
        }
        sx[static_cast<std::size_t>(p)] = sum;
    }
    return sx;
}

void direct_rpa_rhf_energy_kpoints(
    const UniformGrid& grid,
    const std::vector<std::vector<Cplx>>& psi_occ_k,
    const std::vector<std::vector<Cplx>>& psi_virt_k,
    const std::vector<std::vector<double>>& eps_occ_k,
    const std::vector<std::vector<double>>& eps_virt_k,
    const std::vector<std::array<double, 3>>& k_frac,
    const std::vector<double>& k_weights, double& e_corr, double& tr_A,
    std::vector<double>& omega) {
    using linalg::CplxMatrix;

    const int nk = static_cast<int>(k_frac.size());
    if (nk < 1)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf_energy_kpoints: empty k mesh");
    if (static_cast<int>(psi_occ_k.size()) != nk
        || static_cast<int>(psi_virt_k.size()) != nk
        || static_cast<int>(eps_occ_k.size()) != nk
        || static_cast<int>(eps_virt_k.size()) != nk
        || static_cast<int>(k_weights.size()) != nk)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf_energy_kpoints: per-k array length mismatch");

    const int n_occ = static_cast<int>(eps_occ_k[0].size());
    const int n_virt = static_cast<int>(eps_virt_k[0].size());
    if (n_occ < 1 || n_virt < 1)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf_energy_kpoints: need n_occ>=1 and n_virt>=1");
    const long npts = grid.npts();
    const double w_grid = grid.weight();
    const int nov = n_occ * n_virt;

    for (int ik = 0; ik < nk; ++ik) {
        if (static_cast<int>(eps_occ_k[static_cast<std::size_t>(ik)].size())
                != n_occ
            || static_cast<int>(
                   eps_virt_k[static_cast<std::size_t>(ik)].size())
                   != n_virt)
            throw std::invalid_argument(
                "gpw direct_rpa_rhf_energy_kpoints: inconsistent n_occ/n_virt "
                "across k");
        const std::size_t need_o =
            static_cast<std::size_t>(n_occ) * static_cast<std::size_t>(npts);
        const std::size_t need_v =
            static_cast<std::size_t>(n_virt) * static_cast<std::size_t>(npts);
        if (psi_occ_k[static_cast<std::size_t>(ik)].size() != need_o
            || psi_virt_k[static_cast<std::size_t>(ik)].size() != need_v)
            throw std::invalid_argument(
                "gpw direct_rpa_rhf_energy_kpoints: MO grid size mismatch");
    }
    double wsum = 0.0;
    for (double wk : k_weights)
        wsum += wk;
    if (std::abs(wsum - 1.0) > 1e-8)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf_energy_kpoints: k_weights must sum to 1");

    // Product multi-k dRPA: finite-q inter-k dielectric (Furche plasmon).
    // q-points are unique wrapped differences ka − ki over the sampled k mesh
    // (Γ-avoiding MP meshes are not groups under addition, so q ≠ k in
    // general). For each transfer q, Casida from transitions i_k → a_{k+q}
    // for every ik with k+q on the mesh, Coulomb at q (RI-GPW Poisson),
    // Furche-averaged
    //   E_c = Σ_q w_q · ½(Σ_s Ω_{q,s} − Tr A_q),
    // with w_q = n_src(q) / N_k² * N_k = n_src(q)/N_k  so Σ_q w_q = 1 when
    // every ordered pair (ki,ka) contributes once (n_src summed over q = N_k).
    // Recovers the real Γ formula when N_k = 1 and orbitals are real.
    constexpr double kAmBFloor = 1e-10;
    constexpr double kOmega2Tol = 1e-10;
    constexpr double kQTol = 1e-10;

    e_corr = 0.0;
    tr_A = 0.0;
    omega.clear();
    std::vector<Cplx> rho(static_cast<std::size_t>(npts));

    // Unique q = wrap(ka − ki) over all ordered pairs.
    std::vector<std::array<double, 3>> q_list;
    auto q_index_of = [&](double qx, double qy, double qz) -> int {
        const double tx = wrap_frac(qx), ty = wrap_frac(qy),
                     tz = wrap_frac(qz);
        for (std::size_t iq = 0; iq < q_list.size(); ++iq) {
            const double dx = wrap_frac(q_list[iq][0] - tx);
            const double dy = wrap_frac(q_list[iq][1] - ty);
            const double dz = wrap_frac(q_list[iq][2] - tz);
            if (std::abs(dx) < kQTol && std::abs(dy) < kQTol
                && std::abs(dz) < kQTol)
                return static_cast<int>(iq);
        }
        return -1;
    };
    for (int ik = 0; ik < nk; ++ik) {
        for (int ka = 0; ka < nk; ++ka) {
            const double qx = k_frac[static_cast<std::size_t>(ka)][0]
                              - k_frac[static_cast<std::size_t>(ik)][0];
            const double qy = k_frac[static_cast<std::size_t>(ka)][1]
                              - k_frac[static_cast<std::size_t>(ik)][1];
            const double qz = k_frac[static_cast<std::size_t>(ka)][2]
                              - k_frac[static_cast<std::size_t>(ik)][2];
            if (q_index_of(qx, qy, qz) < 0)
                q_list.push_back(
                    {{wrap_frac(qx), wrap_frac(qy), wrap_frac(qz)}});
        }
    }
    const int nq = static_cast<int>(q_list.size());
    if (nq < 1)
        throw std::runtime_error(
            "gpw direct_rpa_rhf_energy_kpoints: empty q mesh from k "
            "differences");

    for (int iq = 0; iq < nq; ++iq) {
        const double q_frac[3] = {q_list[static_cast<std::size_t>(iq)][0],
                                  q_list[static_cast<std::size_t>(iq)][1],
                                  q_list[static_cast<std::size_t>(iq)][2]};
        (void)q_frac; // diagnostics; Bloch ρ already carries e^{iq·r}

        // Sources ik with ka = k[ik] + q on the mesh.
        std::vector<int> src_k;
        std::vector<int> tgt_k;
        src_k.reserve(static_cast<std::size_t>(nk));
        tgt_k.reserve(static_cast<std::size_t>(nk));
        for (int ik = 0; ik < nk; ++ik) {
            const double fx = k_frac[static_cast<std::size_t>(ik)][0]
                              + q_frac[0];
            const double fy = k_frac[static_cast<std::size_t>(ik)][1]
                              + q_frac[1];
            const double fz = k_frac[static_cast<std::size_t>(ik)][2]
                              + q_frac[2];
            const int ka = find_k_index(k_frac, fx, fy, fz);
            if (ka >= 0) {
                src_k.push_back(ik);
                tgt_k.push_back(ka);
            }
        }
        const int n_src = static_cast<int>(src_k.size());
        if (n_src < 1)
            continue;
        // Equal-pair weight: each ordered (ki,ka) contributes 1/N_k once
        // when summed over q of n_src(q)/N_k  (Σ_q n_src = N_k when every
        // ki has exactly one ka per… wait: Σ_q n_src(q) = N_k² ordered
        // pairs that land on mesh = N_k² always since every (ki,ka) defines
        // a q). So Σ_q w_q = Σ_q n_src/N_k = N_k. Need /N_k more → w_q =
        // n_src / N_k².
        const double wq =
            static_cast<double>(n_src)
            / (static_cast<double>(nk) * static_cast<double>(nk));

        const int n_block = n_src * nov;
        std::vector<std::vector<Cplx>> Via(
            static_cast<std::size_t>(n_block));
        std::vector<double> de(static_cast<std::size_t>(n_block));
        std::vector<std::vector<Cplx>> rho_alpha(
            static_cast<std::size_t>(n_block));

        for (int is = 0; is < n_src; ++is) {
            const int ik = src_k[static_cast<std::size_t>(is)];
            const int ka = tgt_k[static_cast<std::size_t>(is)];
            const auto& po = psi_occ_k[static_cast<std::size_t>(ik)];
            const auto& pv = psi_virt_k[static_cast<std::size_t>(ka)];
            const auto& eo = eps_occ_k[static_cast<std::size_t>(ik)];
            const auto& ev = eps_virt_k[static_cast<std::size_t>(ka)];

            for (int i = 0; i < n_occ; ++i) {
                for (int a = 0; a < n_virt; ++a) {
                    const int ia = i * n_virt + a;
                    const int alpha = is * nov + ia;
                    const Cplx* pi =
                        po.data()
                        + static_cast<std::size_t>(i)
                              * static_cast<std::size_t>(npts);
                    const Cplx* pa =
                        pv.data()
                        + static_cast<std::size_t>(a)
                              * static_cast<std::size_t>(npts);
                    rho_alpha[static_cast<std::size_t>(alpha)].resize(
                        static_cast<std::size_t>(npts));
                    for (long g = 0; g < npts; ++g) {
                        const Cplx r = std::conj(pi[g]) * pa[g];
                        rho[static_cast<std::size_t>(g)] = r;
                        rho_alpha[static_cast<std::size_t>(alpha)]
                                 [static_cast<std::size_t>(g)] = r;
                    }
                    // Full Bloch product ρ = ψ_i^{k*} ψ_a^{k+q} already
                    // carries the e^{iq·r} envelope from the AO Bloch sum —
                    // use plain Poisson (G=0 dropped; ρ is charge-neutral by
                    // k-orthogonality). Do NOT also shift by q in G-space
                    // (that double-counts the transfer; same as multi-k MP2).
                    poisson_complex(grid, rho.data(),
                                    Via[static_cast<std::size_t>(alpha)]);

                    de[static_cast<std::size_t>(alpha)] =
                        ev[static_cast<std::size_t>(a)]
                        - eo[static_cast<std::size_t>(i)];
                    if (de[static_cast<std::size_t>(alpha)] <= kAmBFloor)
                        throw std::runtime_error(
                            "gpw direct_rpa_rhf_energy_kpoints: non-positive "
                            "Δε at a (q,k,i,a) transition — dRPA (A−B) not PD");
                }
            }
        }

        // Hermitian Coulomb form (ρ_p | ρ_q) = ∫ V[ρ_p]* ρ_q  (complex PD).
        CplxMatrix Vmat = xt::zeros<Cplx>({n_block, n_block});
        for (int p = 0; p < n_block; ++p) {
            const Cplx* Vp = Via[static_cast<std::size_t>(p)].data();
            for (int qq = 0; qq < n_block; ++qq) {
                const Cplx* rq =
                    rho_alpha[static_cast<std::size_t>(qq)].data();
                Cplx s{0.0, 0.0};
                for (long g = 0; g < npts; ++g)
                    s += std::conj(Vp[g]) * rq[g];
                Vmat(p, qq) = w_grid * s;
            }
        }
        for (int p = 0; p < n_block; ++p) {
            Vmat(p, p) = Cplx{Vmat(p, p).real(), 0.0};
            for (int qq = p + 1; qq < n_block; ++qq) {
                const Cplx avg =
                    0.5 * (Vmat(p, qq) + std::conj(Vmat(qq, p)));
                Vmat(p, qq) = avg;
                Vmat(qq, p) = std::conj(avg);
            }
        }

        CplxMatrix M = xt::zeros<Cplx>({n_block, n_block});
        double tr_Aq = 0.0;
        for (int p = 0; p < n_block; ++p) {
            const double dep = de[static_cast<std::size_t>(p)];
            const double sp = std::sqrt(dep);
            tr_Aq += dep + 2.0 * Vmat(p, p).real();
            for (int qq = 0; qq < n_block; ++qq) {
                const double sq = std::sqrt(de[static_cast<std::size_t>(qq)]);
                Cplx apb = 4.0 * Vmat(p, qq);
                if (p == qq)
                    apb += dep;
                M(p, qq) = sp * apb * sq;
            }
        }
        for (int p = 0; p < n_block; ++p) {
            M(p, p) = Cplx{M(p, p).real(), 0.0};
            for (int qq = p + 1; qq < n_block; ++qq) {
                const Cplx avg = 0.5 * (M(p, qq) + std::conj(M(qq, p)));
                M(p, qq) = avg;
                M(qq, p) = std::conj(avg);
            }
        }

        std::vector<double> w2;
        CplxMatrix Z;
        eigh_hermitian(M, w2, Z);
        (void)Z;

        double sum_om = 0.0;
        for (int s = 0; s < n_block; ++s) {
            double lam = w2[static_cast<std::size_t>(s)];
            if (lam < -kOmega2Tol)
                throw std::runtime_error(
                    "gpw direct_rpa_rhf_energy_kpoints: dRPA omega^2 = "
                    + std::to_string(lam) + " Ha^2 (q=" + std::to_string(iq)
                    + " root " + std::to_string(s)
                    + ") is negative — imaginary excitation (no fallback)");
            if (lam <= kOmega2Tol)
                throw std::runtime_error(
                    "gpw direct_rpa_rhf_energy_kpoints: dRPA omega^2 = "
                    + std::to_string(lam) + " Ha^2 (q=" + std::to_string(iq)
                    + " root " + std::to_string(s)
                    + ") is numerically zero (no fallback)");
            const double om = std::sqrt(lam);
            omega.push_back(om);
            sum_om += om;
        }
        e_corr += wq * 0.5 * (sum_om - tr_Aq);
        tr_A += wq * tr_Aq;
    }
    std::sort(omega.begin(), omega.end());
}

namespace {

constexpr double kGwEta = 1e-8;
constexpr double kGwQpeTol = 1e-8;
constexpr int kGwQpeMaxIter = 100;

/// Diagonal Re Σ_c,p(ω) with complex transition densities stored as |t|².
double sigma_c_diag_abs2(double omega, int p, int n_occ, int n_mo,
                         const std::vector<double>& eps,
                         const std::vector<double>& Om,
                         const std::vector<double>& tdm_abs2) {
    const int n_exc = static_cast<int>(Om.size());
    const std::size_t nmo = static_cast<std::size_t>(n_mo);
    double re = 0.0;
    for (int s = 0; s < n_exc; ++s) {
        const double Os = Om[static_cast<std::size_t>(s)];
        const std::size_t s0 =
            static_cast<std::size_t>(s) * nmo * nmo
            + static_cast<std::size_t>(p);
        for (int m = 0; m < n_occ; ++m) {
            const double t2 = tdm_abs2[s0 + static_cast<std::size_t>(m) * nmo];
            const double den_r = omega + Os - eps[static_cast<std::size_t>(m)];
            const double den2 = den_r * den_r + kGwEta * kGwEta;
            re += t2 * den_r / den2;
        }
        for (int m = n_occ; m < n_mo; ++m) {
            const double t2 = tdm_abs2[s0 + static_cast<std::size_t>(m) * nmo];
            const double den_r = omega - Os - eps[static_cast<std::size_t>(m)];
            const double den2 = den_r * den_r + kGwEta * kGwEta;
            re += t2 * den_r / den2;
        }
    }
    return re;
}

bool solve_qp_abs2(double eps_p, double vx_m_vxc, int p, int n_occ, int n_mo,
                   const std::vector<double>& eps, const std::vector<double>& Om,
                   const std::vector<double>& tdm_abs2, double& omega_out) {
    double w = eps_p;
    constexpr double kH = 1e-6;
    for (int it = 0; it < kGwQpeMaxIter; ++it) {
        const double sig =
            sigma_c_diag_abs2(w, p, n_occ, n_mo, eps, Om, tdm_abs2);
        const double f = w - eps_p - (sig + vx_m_vxc);
        if (std::abs(f) < kGwQpeTol) {
            omega_out = w;
            return true;
        }
        const double sig_p =
            sigma_c_diag_abs2(w + kH, p, n_occ, n_mo, eps, Om, tdm_abs2);
        const double dsig = (sig_p - sig) / kH;
        const double df = 1.0 - dsig;
        if (!(std::abs(df) > 1e-14))
            return false;
        w -= f / df;
    }
    omega_out = w;
    return std::abs(w - eps_p
                    - (sigma_c_diag_abs2(w, p, n_occ, n_mo, eps, Om, tdm_abs2)
                       + vx_m_vxc))
           < 10.0 * kGwQpeTol;
}

/// MT Poisson on a complex density via independent real/imag solves.
void poisson_mt_complex(const UniformGrid& grid, const Cplx* rho,
                        std::vector<Cplx>& V) {
    const long npts = grid.npts();
    std::vector<double> re(static_cast<std::size_t>(npts));
    std::vector<double> im(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g) {
        re[static_cast<std::size_t>(g)] = rho[g].real();
        im[static_cast<std::size_t>(g)] = rho[g].imag();
    }
    const PoissonResult vr = solve_poisson_mt(grid, re);
    const PoissonResult vi = solve_poisson_mt(grid, im);
    V.resize(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g)
        V[static_cast<std::size_t>(g)] =
            Cplx{vr.v_r[static_cast<std::size_t>(g)],
                 vi.v_r[static_cast<std::size_t>(g)]};
}

} // namespace

void spectral_g0w0_rhf_kpoints(
    const UniformGrid& grid,
    const std::vector<std::vector<Cplx>>& psi_mo_k,
    const std::vector<std::vector<double>>& eps_k,
    const std::vector<double>& k_weights, int n_occ, int n_virt,
    const std::vector<double>& vxc_r,
    std::vector<std::vector<double>>& qp_energy_k,
    std::vector<std::vector<double>>& sigma_x_k,
    std::vector<std::vector<double>>& vxc_diag_k, bool& all_converged) {
    using linalg::CplxMatrix;

    const int nk = static_cast<int>(psi_mo_k.size());
    if (nk < 1)
        throw std::invalid_argument(
            "gpw spectral_g0w0_rhf_kpoints: empty k mesh");
    if (static_cast<int>(eps_k.size()) != nk
        || static_cast<int>(k_weights.size()) != nk)
        throw std::invalid_argument(
            "gpw spectral_g0w0_rhf_kpoints: per-k array length mismatch");
    if (n_occ < 1 || n_virt < 1)
        throw std::invalid_argument(
            "gpw spectral_g0w0_rhf_kpoints: need n_occ>=1 and n_virt>=1");

    const long npts = grid.npts();
    const double w_grid = grid.weight();
    if (static_cast<int>(vxc_r.size()) != static_cast<int>(npts))
        throw std::invalid_argument(
            "gpw spectral_g0w0_rhf_kpoints: vxc_r size != npts");

    const int n_mo = static_cast<int>(eps_k[0].size());
    if (n_mo < n_occ + n_virt)
        throw std::invalid_argument(
            "gpw spectral_g0w0_rhf_kpoints: n_mo < n_occ + n_virt");
    const int nov = n_occ * n_virt;
    const std::size_t need_mo =
        static_cast<std::size_t>(n_mo) * static_cast<std::size_t>(npts);

    for (int ik = 0; ik < nk; ++ik) {
        if (static_cast<int>(eps_k[static_cast<std::size_t>(ik)].size())
            != n_mo)
            throw std::invalid_argument(
                "gpw spectral_g0w0_rhf_kpoints: inconsistent n_mo across k");
        if (psi_mo_k[static_cast<std::size_t>(ik)].size() != need_mo)
            throw std::invalid_argument(
                "gpw spectral_g0w0_rhf_kpoints: MO grid size mismatch");
    }
    double wsum = 0.0;
    for (double wk : k_weights)
        wsum += wk;
    if (std::abs(wsum - 1.0) > 1e-8)
        throw std::invalid_argument(
            "gpw spectral_g0w0_rhf_kpoints: k_weights must sum to 1");

    constexpr double kAmBFloor = 1e-10;
    constexpr double kOmega2Tol = 1e-10;
    const double sqrt2 = std::sqrt(2.0);

    qp_energy_k.assign(static_cast<std::size_t>(nk),
                       std::vector<double>(static_cast<std::size_t>(n_mo), 0.0));
    sigma_x_k.assign(static_cast<std::size_t>(nk),
                     std::vector<double>(static_cast<std::size_t>(n_mo), 0.0));
    vxc_diag_k.assign(static_cast<std::size_t>(nk),
                      std::vector<double>(static_cast<std::size_t>(n_mo), 0.0));
    all_converged = true;

    std::vector<Cplx> rho(static_cast<std::size_t>(npts));

    for (int ik = 0; ik < nk; ++ik) {
        const auto& psi = psi_mo_k[static_cast<std::size_t>(ik)];
        const auto& eps = eps_k[static_cast<std::size_t>(ik)];

        // ── OV Coulomb intermediates (q = 0 Poisson) ──────────────────
        std::vector<std::vector<Cplx>> Via(static_cast<std::size_t>(nov));
        for (int i = 0; i < n_occ; ++i) {
            for (int a = 0; a < n_virt; ++a) {
                const int ia = i * n_virt + a;
                const Cplx* pi =
                    psi.data()
                    + static_cast<std::size_t>(i) * static_cast<std::size_t>(npts);
                const Cplx* pa =
                    psi.data()
                    + static_cast<std::size_t>(n_occ + a)
                          * static_cast<std::size_t>(npts);
                for (long g = 0; g < npts; ++g)
                    rho[static_cast<std::size_t>(g)] = std::conj(pi[g]) * pa[g];
                poisson_complex(grid, rho.data(),
                                Via[static_cast<std::size_t>(ia)]);
            }
        }

        auto contract_ov = [&](int i, int a, int j, int b) -> Cplx {
            const int ia = i * n_virt + a;
            const Cplx* Via_g = Via[static_cast<std::size_t>(ia)].data();
            const Cplx* pj =
                psi.data()
                + static_cast<std::size_t>(j) * static_cast<std::size_t>(npts);
            const Cplx* pb =
                psi.data()
                + static_cast<std::size_t>(n_occ + b)
                      * static_cast<std::size_t>(npts);
            Cplx s{0.0, 0.0};
            for (long g = 0; g < npts; ++g)
                s += Via_g[g] * (std::conj(pj[g]) * pb[g]);
            return w_grid * s;
        };

        std::vector<double> de(static_cast<std::size_t>(nov));
        for (int i = 0; i < n_occ; ++i)
            for (int a = 0; a < n_virt; ++a) {
                const int p = i * n_virt + a;
                de[static_cast<std::size_t>(p)] =
                    eps[static_cast<std::size_t>(n_occ + a)]
                    - eps[static_cast<std::size_t>(i)];
                if (de[static_cast<std::size_t>(p)] <= kAmBFloor)
                    throw std::runtime_error(
                        "gpw spectral_g0w0_rhf_kpoints: non-positive Δε at a "
                        "(k,i,a) transition — dRPA (A−B) not PD");
            }

        CplxMatrix Vmat = xt::zeros<Cplx>({nov, nov});
        for (int i = 0; i < n_occ; ++i)
            for (int a = 0; a < n_virt; ++a) {
                const int p = i * n_virt + a;
                for (int j = 0; j < n_occ; ++j)
                    for (int b = 0; b < n_virt; ++b) {
                        const int q = j * n_virt + b;
                        Vmat(p, q) = contract_ov(i, a, j, b);
                    }
            }
        for (int p = 0; p < nov; ++p) {
            Vmat(p, p) = Cplx{Vmat(p, p).real(), 0.0};
            for (int q = p + 1; q < nov; ++q) {
                const Cplx avg = 0.5 * (Vmat(p, q) + std::conj(Vmat(q, p)));
                Vmat(p, q) = avg;
                Vmat(q, p) = std::conj(avg);
            }
        }

        // Hermitian reduction M = (A−B)^{1/2}(A+B)(A−B)^{1/2} with
        // A−B = diag(Δε), A+B = diag(Δε) + 4 V.
        CplxMatrix M = xt::zeros<Cplx>({nov, nov});
        for (int p = 0; p < nov; ++p) {
            const double dep = de[static_cast<std::size_t>(p)];
            const double sp = std::sqrt(dep);
            for (int q = 0; q < nov; ++q) {
                const double sq = std::sqrt(de[static_cast<std::size_t>(q)]);
                Cplx apb = 4.0 * Vmat(p, q);
                if (p == q)
                    apb += dep;
                M(p, q) = sp * apb * sq;
            }
        }
        for (int p = 0; p < nov; ++p) {
            M(p, p) = Cplx{M(p, p).real(), 0.0};
            for (int q = p + 1; q < nov; ++q) {
                const Cplx avg = 0.5 * (M(p, q) + std::conj(M(q, p)));
                M(p, q) = avg;
                M(q, p) = std::conj(avg);
            }
        }

        std::vector<double> w2;
        CplxMatrix Z;
        eigh_hermitian(M, w2, Z);

        std::vector<double> Om(static_cast<std::size_t>(nov));
        // X+Y [s, ia] complex row-major
        std::vector<Cplx> Xpy(static_cast<std::size_t>(nov)
                              * static_cast<std::size_t>(nov));
        for (int s = 0; s < nov; ++s) {
            double lam = w2[static_cast<std::size_t>(s)];
            if (lam < -kOmega2Tol)
                throw std::runtime_error(
                    "gpw spectral_g0w0_rhf_kpoints: dRPA omega^2 = "
                    + std::to_string(lam) + " Ha^2 (k=" + std::to_string(ik)
                    + " root " + std::to_string(s)
                    + ") is negative — imaginary excitation (no fallback)");
            if (lam <= kOmega2Tol)
                throw std::runtime_error(
                    "gpw spectral_g0w0_rhf_kpoints: dRPA omega^2 = "
                    + std::to_string(lam) + " Ha^2 (k=" + std::to_string(ik)
                    + " root " + std::to_string(s)
                    + ") is numerically zero (no fallback)");
            const double om = std::sqrt(lam);
            Om[static_cast<std::size_t>(s)] = om;
            // (X+Y) = Ω^{-1/2} (A−B)^{1/2} z = Ω^{-1/2} sqrt(Δε) z
            const double inv_sqrt_om = 1.0 / std::sqrt(om);
            for (int p = 0; p < nov; ++p) {
                const double sp = std::sqrt(de[static_cast<std::size_t>(p)]);
                Xpy[static_cast<std::size_t>(s)
                        * static_cast<std::size_t>(nov)
                    + static_cast<std::size_t>(p)] =
                    inv_sqrt_om * sp * Z(p, s);
            }
        }

        // Transition |t[s,m,p]|²: t = Σ_ia √2 (X+Y)_ia^s (ia|mp)
        // (ia|mp) = ∫ V_ia · conj(ψ_m) ψ_p
        std::vector<double> tdm_abs2(
            static_cast<std::size_t>(nov) * static_cast<std::size_t>(n_mo)
                * static_cast<std::size_t>(n_mo),
            0.0);
        std::vector<Cplx> Ws(static_cast<std::size_t>(npts));
        for (int s = 0; s < nov; ++s) {
            std::fill(Ws.begin(), Ws.end(), Cplx{0.0, 0.0});
            for (int ia = 0; ia < nov; ++ia) {
                const Cplx z =
                    sqrt2
                    * Xpy[static_cast<std::size_t>(s)
                              * static_cast<std::size_t>(nov)
                          + static_cast<std::size_t>(ia)];
                if (z.real() == 0.0 && z.imag() == 0.0)
                    continue;
                const Cplx* Via_g = Via[static_cast<std::size_t>(ia)].data();
                for (long g = 0; g < npts; ++g)
                    Ws[static_cast<std::size_t>(g)] += z * Via_g[g];
            }
            const std::size_t s0 = static_cast<std::size_t>(s)
                                   * static_cast<std::size_t>(n_mo)
                                   * static_cast<std::size_t>(n_mo);
            for (int m = 0; m < n_mo; ++m) {
                const Cplx* pm =
                    psi.data()
                    + static_cast<std::size_t>(m)
                          * static_cast<std::size_t>(npts);
                for (int p = 0; p < n_mo; ++p) {
                    const Cplx* pp =
                        psi.data()
                        + static_cast<std::size_t>(p)
                              * static_cast<std::size_t>(npts);
                    Cplx acc{0.0, 0.0};
                    for (long g = 0; g < npts; ++g)
                        acc += Ws[static_cast<std::size_t>(g)]
                               * (std::conj(pm[g]) * pp[g]);
                    acc *= w_grid;
                    tdm_abs2[s0
                             + static_cast<std::size_t>(m)
                                   * static_cast<std::size_t>(n_mo)
                             + static_cast<std::size_t>(p)] =
                        std::norm(acc); // |t|²
                }
            }
        }

        // MT exchange diagonals Σ_x,p = −Σ_i Re(pi|pi)
        auto& sx = sigma_x_k[static_cast<std::size_t>(ik)];
        std::vector<Cplx> Vmt;
        for (int p = 0; p < n_mo; ++p) {
            const Cplx* pp =
                psi.data()
                + static_cast<std::size_t>(p) * static_cast<std::size_t>(npts);
            double sum = 0.0;
            for (int i = 0; i < n_occ; ++i) {
                const Cplx* pi =
                    psi.data()
                    + static_cast<std::size_t>(i)
                          * static_cast<std::size_t>(npts);
                for (long g = 0; g < npts; ++g)
                    rho[static_cast<std::size_t>(g)] =
                        std::conj(pp[g]) * pi[g];
                poisson_mt_complex(grid, rho.data(), Vmt);
                Cplx s{0.0, 0.0};
                for (long g = 0; g < npts; ++g)
                    s += Vmt[static_cast<std::size_t>(g)]
                         * rho[static_cast<std::size_t>(g)];
                sum -= w_grid * s.real();
            }
            sx[static_cast<std::size_t>(p)] = sum;
        }

        // Local v_xc diagonals
        auto& vxc_d = vxc_diag_k[static_cast<std::size_t>(ik)];
        for (int p = 0; p < n_mo; ++p) {
            const Cplx* pp =
                psi.data()
                + static_cast<std::size_t>(p) * static_cast<std::size_t>(npts);
            double acc = 0.0;
            for (long g = 0; g < npts; ++g)
                acc += std::norm(pp[g]) * vxc_r[static_cast<std::size_t>(g)];
            vxc_d[static_cast<std::size_t>(p)] = w_grid * acc;
        }

        // Newton QP per orbital
        auto& qp = qp_energy_k[static_cast<std::size_t>(ik)];
        for (int p = 0; p < n_mo; ++p) {
            const double vx_m_vxc = sx[static_cast<std::size_t>(p)]
                                   - vxc_d[static_cast<std::size_t>(p)];
            double w = eps[static_cast<std::size_t>(p)];
            const bool ok =
                solve_qp_abs2(eps[static_cast<std::size_t>(p)], vx_m_vxc, p,
                              n_occ, n_mo, eps, Om, tdm_abs2, w);
            if (!ok)
                all_converged = false;
            qp[static_cast<std::size_t>(p)] = w;
        }
    }
    (void)k_weights; // k-diagonal Σ_c; weights validated above
}

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
