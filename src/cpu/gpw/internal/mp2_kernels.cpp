/// @file src/cpu/gpw/internal/mp2_kernels.cpp
/// @brief GPW MP2 kernels — MO-on-grid, (ia|jb) via per-pair Poisson or
///        RI-GPW G-space Coulomb intermediates, and the closed-shell
///        doubles sum. See mp2_kernels.h.

#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"

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
void poisson_complex(const UniformGrid& grid, const Cplx* rho,
                     std::vector<Cplx>& V) {
    const long npts = grid.npts();
    std::vector<double> re(static_cast<std::size_t>(npts));
    std::vector<double> im(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g) {
        re[static_cast<std::size_t>(g)] = rho[g].real();
        im[static_cast<std::size_t>(g)] = rho[g].imag();
    }
    const PoissonResult vr = solve_poisson(grid, re);
    const PoissonResult vi = solve_poisson(grid, im);
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

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
