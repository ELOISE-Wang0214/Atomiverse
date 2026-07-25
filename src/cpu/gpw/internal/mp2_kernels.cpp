/// @file src/cpu/gpw/internal/mp2_kernels.cpp
/// @brief GPW MP2 kernels — MO-on-grid, (ia|jb) via per-pair Poisson, and
///        the canonical closed-shell doubles sum. See mp2_kernels.h.

#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"

#include <cstddef>
#include <stdexcept>
#include <string>

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
            "gpw mo_on_grid: MO range [" + std::to_string(p0) + ", " +
            std::to_string(p0 + np) + ") outside 0.." + std::to_string(n_mo));
    const long npts = cbg.npts;
    std::vector<double> psi(static_cast<std::size_t>(np) *
                                static_cast<std::size_t>(npts),
                            0.0);
    const std::size_t na = cbg.gidx.size();
    for (std::size_t k = 0; k < na; ++k) {
        const long g = cbg.gidx[k];
        for (std::size_t nz = cbg.off[k]; nz < cbg.off[k + 1]; ++nz) {
            const int u = cbg.fu[nz];
            const double val = cbg.fv[nz];
            for (int p = 0; p < np; ++p)
                psi[static_cast<std::size_t>(p) *
                        static_cast<std::size_t>(npts) +
                    static_cast<std::size_t>(g)] += C(u, p0 + p) * val;
        }
    }
    return psi;
}

std::vector<double> mo_eri_ovov(const UniformGrid& grid,
                                const std::vector<double>& psi_occ, int n_occ,
                                const std::vector<double>& psi_virt,
                                int n_virt, bool nonperiodic) {
    const long npts = grid.npts();
    const double w = grid.weight();
    if (psi_occ.size() != static_cast<std::size_t>(n_occ) *
                              static_cast<std::size_t>(npts) ||
        psi_virt.size() != static_cast<std::size_t>(n_virt) *
                               static_cast<std::size_t>(npts))
        throw std::invalid_argument(
            "gpw mo_eri_ovov: MO grid arrays do not match n_occ/n_virt/npts");
    const int nov = n_occ * n_virt;
    std::vector<double> v(static_cast<std::size_t>(nov) *
                          static_cast<std::size_t>(nov));
    std::vector<double> rho(static_cast<std::size_t>(npts));
    for (int ia = 0; ia < nov; ++ia) {
        const double* pi =
            psi_occ.data() + static_cast<std::size_t>(ia / n_virt) *
                                 static_cast<std::size_t>(npts);
        const double* pa =
            psi_virt.data() + static_cast<std::size_t>(ia % n_virt) *
                                  static_cast<std::size_t>(npts);
        for (long g = 0; g < npts; ++g)
            rho[static_cast<std::size_t>(g)] = pi[g] * pa[g];
        const PoissonResult pr =
            nonperiodic ? solve_poisson_mt(grid, rho) : solve_poisson(grid, rho);
        // (ia|jb) = (jb|ia): fill the upper triangle jb >= ia and mirror.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int jb = ia; jb < nov; ++jb) {
            const double* pj =
                psi_occ.data() + static_cast<std::size_t>(jb / n_virt) *
                                     static_cast<std::size_t>(npts);
            const double* pb =
                psi_virt.data() + static_cast<std::size_t>(jb % n_virt) *
                                      static_cast<std::size_t>(npts);
            double s = 0.0;
            for (long g = 0; g < npts; ++g)
                s += pr.v_r[static_cast<std::size_t>(g)] * pj[g] * pb[g];
            const double val = w * s;
            v[static_cast<std::size_t>(ia) * static_cast<std::size_t>(nov) +
              static_cast<std::size_t>(jb)] = val;
            v[static_cast<std::size_t>(jb) * static_cast<std::size_t>(nov) +
              static_cast<std::size_t>(ia)] = val;
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
    const std::size_t nov = static_cast<std::size_t>(no) *
                            static_cast<std::size_t>(nv);
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
                        v_ovov[(static_cast<std::size_t>(i) * nv + a) * nov +
                               static_cast<std::size_t>(j) * nv + b];
                    const double v_ibja =
                        v_ovov[(static_cast<std::size_t>(i) * nv + b) * nov +
                               static_cast<std::size_t>(j) * nv + a];
                    const double denom =
                        eps_occ[static_cast<std::size_t>(i)] +
                        eps_occ[static_cast<std::size_t>(j)] -
                        eps_virt[static_cast<std::size_t>(a)] -
                        eps_virt[static_cast<std::size_t>(b)];
                    e_os += v_iajb * v_iajb / denom;
                    e_ss += v_iajb * (v_iajb - v_ibja) / denom;
                }
}

} // namespace internal
} // namespace gpw
} // namespace cpu
} // namespace ATV
