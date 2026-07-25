/// @file tests/cpu/gpw/test_cpu_gpw_mp2.cpp
/// @brief GPW MP2 kernel units: closed-shell doubles sum vs an explicit
///        spin-orbital MP2 sum, MO-on-grid CSR expansion vs the dense AO
///        evaluation, (ia|jb) vs the Poisson Hartree energy, and the
///        GpwDft::mp2 guard throws. No SCF is run here (unit tier).

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"
#include "atomiverse/cpu/gpw/internal/poisson_fft.h"
#include "atomiverse/cpu/linalg/matrix.h"

using ATV::cpu::gpw::build_compact_basis_grid;
using ATV::cpu::gpw::CompactBasisGrid;
using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::make_uniform_grid;
using ATV::cpu::gpw::solve_poisson;
using ATV::cpu::gpw::UniformGrid;
using ATV::cpu::gpw::internal::mo_eri_ovov;
using ATV::cpu::gpw::internal::mo_eri_ovov_ri;
using ATV::cpu::gpw::internal::mo_on_grid;
using ATV::cpu::gpw::internal::mp2_rhf_energy;
using ATV::cpu::linalg::RealMatrix;

namespace {

// One normalized-ish s Gaussian centered at c.
GridBasisFunction s_gauss(double alpha, double cx, double cy, double cz,
                          double coeff = 1.0) {
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = cx;
    p.center[1] = cy;
    p.center[2] = cz;
    p.coeff = coeff;
    return GridBasisFunction{p};
}

} // namespace

// Closed-shell doubles sum == explicit spin-orbital MP2
// E = 1/4 Σ_{ijab,spin} |<ij||ab>|² / Δ on the same synthetic integrals.
TEST(GpwMp2, ClosedShellSumMatchesSpinOrbitalSum) {
    const int no = 2, nv = 2, nov = no * nv;
    std::vector<double> eps_occ = {-0.9, -0.5};
    std::vector<double> eps_virt = {0.2, 0.7};
    // Synthetic (ia|jb) with the physical symmetry (ia|jb) = (jb|ia).
    std::vector<double> v(static_cast<std::size_t>(nov) * nov, 0.0);
    auto idx = [&](int i, int a, int j, int b) {
        return static_cast<std::size_t>((i * nv + a) * nov + j * nv + b);
    };
    unsigned s = 12345u;
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a)
            for (int j = 0; j < no; ++j)
                for (int b = 0; b < nv; ++b) {
                    s = s * 1664525u + 1013904223u;
                    const double x = 0.05
                                     * (static_cast<double>(s >> 8)
                                            / static_cast<double>(1u << 24)
                                        - 0.5);
                    v[idx(i, a, j, b)] = x;
                    v[idx(j, b, i, a)] = x; // enforce (ia|jb) = (jb|ia)
                }

    double e_os = 0.0, e_ss = 0.0;
    mp2_rhf_energy(v, eps_occ, eps_virt, e_os, e_ss);
    const double e_closed = e_os + e_ss;

    // Spin-orbital sum: spatial p → (p,alpha),(p,beta); <ij|ab> =
    // (ia|jb) δ_{σi σa} δ_{σj σb}.
    const int nso_o = 2 * no, nso_v = 2 * nv;
    double e_so = 0.0;
    for (int I = 0; I < nso_o; ++I)
        for (int J = 0; J < nso_o; ++J)
            for (int A = 0; A < nso_v; ++A)
                for (int B = 0; B < nso_v; ++B) {
                    const int i = I / 2, si = I % 2;
                    const int j = J / 2, sj = J % 2;
                    const int a = A / 2, sa = A % 2;
                    const int b = B / 2, sb = B % 2;
                    const double dir =
                        (si == sa && sj == sb) ? v[idx(i, a, j, b)] : 0.0;
                    const double exc =
                        (si == sb && sj == sa) ? v[idx(i, b, j, a)] : 0.0;
                    const double anti = dir - exc;
                    const double denom =
                        eps_occ[static_cast<std::size_t>(i)]
                        + eps_occ[static_cast<std::size_t>(j)]
                        - eps_virt[static_cast<std::size_t>(a)]
                        - eps_virt[static_cast<std::size_t>(b)];
                    e_so += 0.25 * anti * anti / denom;
                }
    EXPECT_NEAR(e_closed, e_so, 1e-12);
    // Unscaled decomposition is exact: e_corr = e_os + e_ss by definition.
    EXPECT_LT(e_closed, 0.0); // gapped spectrum → strictly negative
}

// mo_on_grid over the compact representation == an explicit Γ-point
// image-summed AO evaluation (independent in-test reference,
// χ_u(r) = Σ_T φ_u(r − T)) contracted with the same coefficients.
TEST(GpwMp2, MoOnGridMatchesImageSummedContraction) {
    const double cell[3] = {6.0, 6.0, 6.0};
    const UniformGrid grid = make_uniform_grid(cell, 30.0);
    std::vector<GridBasisFunction> basis;
    basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    basis.push_back(s_gauss(0.5, 3.8, 3.0, 3.0));
    const CompactBasisGrid cbg =
        build_compact_basis_grid(grid, basis, /*periodic=*/true, 1e-14);
    RealMatrix C = xt::zeros<double>({2, 2});
    C(0, 0) = 0.7;
    C(1, 0) = -0.3;
    C(0, 1) = 0.2;
    C(1, 1) = 0.9;
    const std::vector<double> psi = mo_on_grid(cbg, C, 0, 2);

    // Explicit lattice-image sum over T ∈ [-2,2]³ (Gaussian radius at
    // eps=1e-14 is well inside 2 cells for these exponents).
    auto chi = [&](const GridBasisFunction& bf, const double r[3]) {
        double v = 0.0;
        for (int tx = -2; tx <= 2; ++tx)
            for (int ty = -2; ty <= 2; ++ty)
                for (int tz = -2; tz <= 2; ++tz)
                    for (const auto& p : bf) {
                        const double dx = r[0] - p.center[0] - tx * cell[0];
                        const double dy = r[1] - p.center[1] - ty * cell[1];
                        const double dz = r[2] - p.center[2] - tz * cell[2];
                        v += p.coeff
                             * std::exp(-p.alpha
                                        * (dx * dx + dy * dy + dz * dz));
                    }
        return v;
    };
    const long npts = grid.npts();
    double max_err = 0.0;
    for (long g = 0; g < npts; ++g) {
        double r[3];
        grid.point(g, r);
        const double c0 = chi(basis[0], r), c1 = chi(basis[1], r);
        for (int p = 0; p < 2; ++p) {
            const double ref = c0 * C(0, p) + c1 * C(1, p);
            const double err =
                std::abs(psi[static_cast<std::size_t>(p) * npts + g] - ref);
            max_err = std::max(max_err, err);
        }
    }
    EXPECT_LT(max_err, 1e-9);
}

// (pp|pp) from the pair-integral route == 2 × Poisson Hartree energy of the
// same pair density (E_H = 0.5 Σ w ρ V ⇒ (pp|pp) = Σ w ρ V = 2 E_H).
TEST(GpwMp2, MoEriMatchesPoissonHartreeEnergy) {
    const double cell[3] = {6.0, 6.0, 6.0};
    const UniformGrid grid = make_uniform_grid(cell, 30.0);
    const long npts = grid.npts();
    // Hand-made "MOs": two localized smooth functions.
    std::vector<double> psi_o(static_cast<std::size_t>(npts));
    std::vector<double> psi_v(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g) {
        double r[3];
        grid.point(g, r);
        const double dx = r[0] - 3.0, dy = r[1] - 3.0, dz = r[2] - 3.0;
        const double r2 = dx * dx + dy * dy + dz * dz;
        psi_o[static_cast<std::size_t>(g)] = std::exp(-0.7 * r2);
        psi_v[static_cast<std::size_t>(g)] = dx * std::exp(-0.5 * r2);
    }
    const std::vector<double> v =
        mo_eri_ovov(grid, psi_o, 1, psi_v, 1, /*nonperiodic=*/false);
    ASSERT_EQ(v.size(), 1u);

    std::vector<double> rho(static_cast<std::size_t>(npts));
    for (long g = 0; g < npts; ++g)
        rho[static_cast<std::size_t>(g)] = psi_o[static_cast<std::size_t>(g)]
                                           * psi_v[static_cast<std::size_t>(g)];
    const auto pr = solve_poisson(grid, rho);
    EXPECT_NEAR(v[0], 2.0 * pr.e_hartree, 1e-10);
    EXPECT_GT(v[0], 0.0); // a self-repulsion integral is positive
}

// RI-GPW G-space (ia|jb) ≡ Poisson real-space tensor on the same grid
// (Parseval identity; product path must match the foothold to ~1e-10).
TEST(GpwMp2, MoEriRiMatchesPoisson) {
    const double cell[3] = {6.0, 6.0, 6.0};
    const UniformGrid grid = make_uniform_grid(cell, 30.0);
    const long npts = grid.npts();
    const int no = 2, nv = 2;
    std::vector<double> psi_o(static_cast<std::size_t>(no) * npts);
    std::vector<double> psi_v(static_cast<std::size_t>(nv) * npts);
    for (long g = 0; g < npts; ++g) {
        double r[3];
        grid.point(g, r);
        const double dx = r[0] - 3.0, dy = r[1] - 3.0, dz = r[2] - 3.0;
        const double r2 = dx * dx + dy * dy + dz * dz;
        psi_o[static_cast<std::size_t>(0) * npts + g] = std::exp(-0.7 * r2);
        psi_o[static_cast<std::size_t>(1) * npts + g] =
            dy * std::exp(-0.55 * r2);
        psi_v[static_cast<std::size_t>(0) * npts + g] = dx * std::exp(-0.5 * r2);
        psi_v[static_cast<std::size_t>(1) * npts + g] =
            dz * std::exp(-0.45 * r2);
    }
    const std::vector<double> v_pois =
        mo_eri_ovov(grid, psi_o, no, psi_v, nv, /*nonperiodic=*/false);
    const std::vector<double> v_ri = mo_eri_ovov_ri(grid, psi_o, no, psi_v, nv);
    ASSERT_EQ(v_pois.size(), v_ri.size());
    double max_err = 0.0;
    for (std::size_t i = 0; i < v_pois.size(); ++i)
        max_err = std::max(max_err, std::abs(v_pois[i] - v_ri[i]));
    EXPECT_LT(max_err, 1e-10) << "max |RI − Poisson| = " << max_err;
}

// Guard: mp2() before build() throws logic_error.
TEST(GpwMp2, ThrowsBeforeBuild) {
    GpwDft m;
    EXPECT_THROW(m.mp2({}), std::logic_error);
}

// Guard: UKS configuration is rejected (closed-shell foothold).
TEST(GpwMp2, ThrowsForUks) {
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 6.0;
    m.n_alpha = 1;
    m.n_beta = 1;
    GpwAtom a;
    a.z_val = 2.0;
    a.basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    m.build({a});
    EXPECT_THROW(m.mp2({a}), std::invalid_argument);
}

// Guard: converged run() required before correlation.
TEST(GpwMp2, ThrowsWithoutConvergedRun) {
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 6.0;
    m.n_occ = 1;
    GpwAtom a;
    a.z_val = 2.0;
    a.basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    m.build({a});
    EXPECT_THROW(m.mp2({a}), std::logic_error);
}

// Complex OS/SS reduction: real integrals → same e_os/e_ss as mp2_rhf_energy.
TEST(GpwMp2, ComplexOsSsReducesToReal) {
    const int no = 2, nv = 2, nov = no * nv;
    std::vector<double> eps_occ = {-1.0, -0.4};
    std::vector<double> eps_virt = {0.3, 0.8};
    std::vector<double> v(static_cast<std::size_t>(nov) * nov, 0.0);
    auto idx = [&](int i, int a, int j, int b) {
        return static_cast<std::size_t>((i * nv + a) * nov + j * nv + b);
    };
    v[idx(0, 0, 0, 0)] = 0.11;
    v[idx(0, 0, 0, 1)] = 0.03;
    v[idx(0, 0, 1, 0)] = 0.03;
    v[idx(0, 0, 1, 1)] = 0.07;
    v[idx(0, 1, 0, 0)] = 0.03;
    v[idx(0, 1, 0, 1)] = 0.09;
    v[idx(0, 1, 1, 0)] = 0.02;
    v[idx(0, 1, 1, 1)] = 0.05;
    v[idx(1, 0, 0, 0)] = 0.03;
    v[idx(1, 0, 0, 1)] = 0.02;
    v[idx(1, 0, 1, 0)] = 0.08;
    v[idx(1, 0, 1, 1)] = 0.04;
    v[idx(1, 1, 0, 0)] = 0.07;
    v[idx(1, 1, 0, 1)] = 0.05;
    v[idx(1, 1, 1, 0)] = 0.04;
    v[idx(1, 1, 1, 1)] = 0.06;
    // Symmetrize (ia|jb)=(jb|ia)
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a)
            for (int j = 0; j < no; ++j)
                for (int b = 0; b < nv; ++b)
                    v[idx(j, b, i, a)] = v[idx(i, a, j, b)];

    double e_os = 0, e_ss = 0;
    mp2_rhf_energy(v, eps_occ, eps_virt, e_os, e_ss);

    double e_os_c = 0, e_ss_c = 0;
    for (int i = 0; i < no; ++i)
        for (int j = 0; j < no; ++j)
            for (int a = 0; a < nv; ++a)
                for (int b = 0; b < nv; ++b) {
                    const std::complex<double> v_iajb = v[idx(i, a, j, b)];
                    const std::complex<double> v_ibja = v[idx(i, b, j, a)];
                    const double denom = eps_occ[i] + eps_occ[j] - eps_virt[a]
                                         - eps_virt[b];
                    e_os_c += std::norm(v_iajb) / denom;
                    e_ss_c += (v_iajb * std::conj(v_iajb - v_ibja)).real() / denom;
                }
    EXPECT_NEAR(e_os_c, e_os, 1e-14);
    EXPECT_NEAR(e_ss_c, e_ss, 1e-14);
}

// mo_on_grid_k at Γ (real phases) matches real mo_on_grid on the same AOs.
TEST(GpwMp2, MoOnGridKAtGammaMatchesReal) {
    using ATV::cpu::gpw::internal::build_kpoint_basis_grid;
    using ATV::cpu::gpw::internal::KPoint;
    using ATV::cpu::gpw::internal::mo_on_grid_k;
    using ATV::cpu::linalg::CplxMatrix;

    const double cell[3] = {6.0, 6.0, 6.0};
    const UniformGrid grid = make_uniform_grid(cell, 30.0);
    std::vector<GridBasisFunction> basis;
    basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    basis.push_back(s_gauss(0.5, 3.8, 3.0, 3.0));
    const CompactBasisGrid cbg =
        build_compact_basis_grid(grid, basis, /*periodic=*/true, 1e-14);

    RealMatrix C = xt::zeros<double>({2, 2});
    C(0, 0) = 0.7;
    C(1, 0) = -0.3;
    C(0, 1) = 0.2;
    C(1, 1) = 0.9;
    const std::vector<double> psi_r = mo_on_grid(cbg, C, 0, 2);

    std::vector<KPoint> kpts(1);
    kpts[0].frac[0] = kpts[0].frac[1] = kpts[0].frac[2] = 0.0;
    kpts[0].weight = 1.0;
    const auto kbg = build_kpoint_basis_grid(grid, basis, kpts, 1e-14);
    CplxMatrix Ck = xt::zeros<std::complex<double>>({2, 2});
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            Ck(i, j) = C(i, j);
    const std::vector<std::complex<double>> psi_k = mo_on_grid_k(kbg, 0, Ck, 0, 2);

    const long npts = grid.npts();
    double max_err = 0.0;
    for (std::size_t i = 0; i < psi_r.size(); ++i) {
        max_err = std::max(max_err, std::abs(psi_k[i].real() - psi_r[i]));
        max_err = std::max(max_err, std::abs(psi_k[i].imag()));
    }
    EXPECT_LT(max_err, 1e-10) << "max |complex Γ − real| = " << max_err;
    (void)npts;
}
