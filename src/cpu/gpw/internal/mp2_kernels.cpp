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

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/gpw/internal/mt_poisson.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"

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

void direct_rpa_rhf_energy(const std::vector<double>& v_ovov,
                           const std::vector<double>& eps_occ,
                           const std::vector<double>& eps_virt, double& e_corr,
                           double& tr_A, std::vector<double>& omega) {
    using linalg::RealMatrix;
    using linalg::sym_power;

    const int no = static_cast<int>(eps_occ.size());
    const int nv = static_cast<int>(eps_virt.size());
    if (no < 1 || nv < 1)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf_energy: need n_occ >= 1 and n_virt >= 1");
    const std::size_t nov =
        static_cast<std::size_t>(no) * static_cast<std::size_t>(nv);
    if (v_ovov.size() != nov * nov)
        throw std::invalid_argument(
            "gpw direct_rpa_rhf_energy: tensor size does not match "
            "eps_occ/eps_virt");

    // Singlet dRPA (α_hf = 0): A = Δε δ + 2 V, B = 2 V, index p = i*nv + a.
    RealMatrix A = xt::zeros<double>({static_cast<int>(nov), static_cast<int>(nov)});
    RealMatrix B = xt::zeros<double>({static_cast<int>(nov), static_cast<int>(nov)});
    tr_A = 0.0;
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a) {
            const int p = i * nv + a;
            const double de =
                eps_virt[static_cast<std::size_t>(a)]
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
            tr_A += A(p, p);
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
                "gpw direct_rpa_rhf_energy: (A - B) is not positive definite "
                "(lowest eigenvalue = "
                + std::to_string(s(0))
                + " Ha) — dRPA reduction undefined (no fallback)");
    }
    const RealMatrix T =
        sym_power(Mm, 0.5, kAmBFloor, "gpw direct_rpa_rhf_energy");
    RealMatrix M = xt::linalg::dot(T, RealMatrix(xt::linalg::dot(Mp, T)));
    // Symmetrize residual floating noise before the Hermitian eigendecomp.
    for (std::size_t r = 0; r < nov; ++r)
        for (std::size_t c = r + 1; c < nov; ++c) {
            const double v = 0.5 * (M(r, c) + M(c, r));
            M(r, c) = v;
            M(c, r) = v;
        }
    auto [w2, Z] = xt::linalg::eigh(M);
    (void)Z;

    omega.assign(nov, 0.0);
    double sum_om = 0.0;
    for (std::size_t s = 0; s < nov; ++s) {
        if (w2(s) < -kOmega2Tol)
            throw std::runtime_error(
                "gpw direct_rpa_rhf_energy: dRPA omega^2 = "
                + std::to_string(w2(s)) + " Ha^2 (root " + std::to_string(s)
                + ") is negative — imaginary excitation (no fallback)");
        if (w2(s) <= kOmega2Tol)
            throw std::runtime_error(
                "gpw direct_rpa_rhf_energy: dRPA omega^2 = "
                + std::to_string(w2(s)) + " Ha^2 (root " + std::to_string(s)
                + ") is numerically zero (no fallback)");
        const double w = std::sqrt(w2(s));
        omega[s] = w;
        sum_om += w;
    }
    e_corr = 0.5 * (sum_om - tr_A);
}

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
