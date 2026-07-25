/// @file tests/cpu/gpw/test_cpu_gpw_rpa.cpp
/// @brief GPW dRPA kernel units: Furche plasmon formula on a synthetic
///        (ia|jb) tensor vs an independent dense A/B path, and
///        GpwDft::direct_rpa guard throws. No SCF (unit tier).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::internal::direct_rpa_rhf_energy;
using ATV::cpu::linalg::RealMatrix;
using ATV::cpu::linalg::sym_power;

namespace {

GridBasisFunction s_gauss(double alpha, double cx, double cy, double cz) {
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = cx;
    p.center[1] = cy;
    p.center[2] = cz;
    p.coeff = 1.0;
    return GridBasisFunction{p};
}

} // namespace

namespace {

// Independent dense A/B plasmon path (mirrors molecular OrbitalHessian
// conventions with GPW composite index p = i*nv + a).
double independent_plasmon(const std::vector<double>& v,
                           const std::vector<double>& eps_occ,
                           const std::vector<double>& eps_virt) {
    const int no = static_cast<int>(eps_occ.size());
    const int nv = static_cast<int>(eps_virt.size());
    const int nov = no * nv;
    RealMatrix A = xt::zeros<double>({nov, nov});
    RealMatrix B = xt::zeros<double>({nov, nov});
    double tr_A = 0.0;
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a) {
            const int p = i * nv + a;
            const double de = eps_virt[static_cast<std::size_t>(a)]
                              - eps_occ[static_cast<std::size_t>(i)];
            for (int j = 0; j < no; ++j)
                for (int b = 0; b < nv; ++b) {
                    const int q = j * nv + b;
                    const double two_v =
                        2.0 * v[static_cast<std::size_t>(p) * nov
                                + static_cast<std::size_t>(q)];
                    double Aval = two_v;
                    if (p == q)
                        Aval += de;
                    A(p, q) = Aval;
                    B(p, q) = two_v;
                }
            tr_A += A(p, p);
        }
    const RealMatrix Mp = A + B;
    const RealMatrix Mm = A - B;
    const RealMatrix T = sym_power(Mm, 0.5, 1e-10, "test");
    RealMatrix M = xt::linalg::dot(T, RealMatrix(xt::linalg::dot(Mp, T)));
    for (int r = 0; r < nov; ++r)
        for (int c = r + 1; c < nov; ++c) {
            const double vv = 0.5 * (M(r, c) + M(c, r));
            M(r, c) = vv;
            M(c, r) = vv;
        }
    auto [w2, Z] = xt::linalg::eigh(M);
    (void)Z;
    double sum_om = 0.0;
    for (int s = 0; s < nov; ++s)
        sum_om += std::sqrt(w2(s));
    return 0.5 * (sum_om - tr_A);
}

} // namespace

TEST(GpwRpa, PlasmonMatchesIndependentDensePath) {
    const int no = 2, nv = 3, nov = no * nv;
    std::vector<double> eps_occ = {-1.1, -0.4};
    std::vector<double> eps_virt = {0.15, 0.55, 1.2};
    std::vector<double> v(static_cast<std::size_t>(nov) * nov, 0.0);
    auto idx = [&](int i, int a, int j, int b) {
        return static_cast<std::size_t>((i * nv + a) * nov + j * nv + b);
    };
    unsigned s = 99191u;
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a)
            for (int j = 0; j < no; ++j)
                for (int b = 0; b < nv; ++b) {
                    s = s * 1664525u + 1013904223u;
                    const double x =
                        0.03
                        * (static_cast<double>(s >> 8)
                               / static_cast<double>(1u << 24)
                           - 0.5);
                    v[idx(i, a, j, b)] = x;
                    v[idx(j, b, i, a)] = x;
                }

    double e_corr = 0.0, tr_A = 0.0;
    std::vector<double> omega;
    direct_rpa_rhf_energy(v, eps_occ, eps_virt, e_corr, tr_A, omega);

    EXPECT_EQ(omega.size(), static_cast<std::size_t>(nov));
    EXPECT_NEAR(e_corr, independent_plasmon(v, eps_occ, eps_virt), 1e-12);
    // Gapped spectrum with Coulomb → correlation is finite and typically
    // negative for physical (ia|jb); synthetic noise may not guarantee
    // sign, but Tr A and sum Ω are positive and finite.
    EXPECT_TRUE(std::isfinite(e_corr));
    EXPECT_GT(tr_A, 0.0);
    double sum_om = 0.0;
    for (double w : omega) {
        EXPECT_GT(w, 0.0);
        sum_om += w;
    }
    EXPECT_NEAR(e_corr, 0.5 * (sum_om - tr_A), 1e-12);
}

// Zero integrals ⇒ Ec = 0 (Ω = diagonal of A = Δε).
TEST(GpwRpa, ZeroIntegralsGiveZeroCorrelation) {
    const int no = 2, nv = 2, nov = no * nv;
    std::vector<double> eps_occ = {-0.8, -0.3};
    std::vector<double> eps_virt = {0.2, 0.9};
    std::vector<double> v(static_cast<std::size_t>(nov) * nov, 0.0);
    double e_corr = 0.0, tr_A = 0.0;
    std::vector<double> omega;
    direct_rpa_rhf_energy(v, eps_occ, eps_virt, e_corr, tr_A, omega);
    EXPECT_NEAR(e_corr, 0.0, 1e-12);
    // Ω_s = ε_a − ε_i for each pair (ascending after eigh).
    std::vector<double> expected;
    for (int i = 0; i < no; ++i)
        for (int a = 0; a < nv; ++a)
            expected.push_back(eps_virt[static_cast<std::size_t>(a)]
                               - eps_occ[static_cast<std::size_t>(i)]);
    std::sort(expected.begin(), expected.end());
    ASSERT_EQ(omega.size(), expected.size());
    for (std::size_t k = 0; k < omega.size(); ++k)
        EXPECT_NEAR(omega[k], expected[k], 1e-12);
}

TEST(GpwRpa, DirectRpaBeforeBuildThrows) {
    GpwDft m;
    EXPECT_THROW(m.direct_rpa({}), std::logic_error);
}

TEST(GpwRpa, DirectRpaUksThrows) {
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 6.0;
    m.n_alpha = 1;
    m.n_beta = 1;
    GpwAtom a;
    a.z_val = 2.0;
    a.basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    m.build({a});
    EXPECT_THROW(m.direct_rpa({a}), std::invalid_argument);
}

TEST(GpwRpa, DirectRpaUnconvergedThrows) {
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 6.0;
    m.n_occ = 1;
    GpwAtom a;
    a.z_val = 2.0;
    a.basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    m.build({a});
    EXPECT_THROW(m.direct_rpa({a}), std::logic_error);
}

// Multi-k kernel at N_k=1 with real-valued complex MOs recovers the real
// Γ plasmon formula on the same pair densities (no SCF).
TEST(GpwRpa, KpointsNk1RecoversRealPlasmon) {
    using ATV::cpu::gpw::make_uniform_grid;
    using ATV::cpu::gpw::UniformGrid;
    using ATV::cpu::gpw::internal::direct_rpa_rhf_energy;
    using ATV::cpu::gpw::internal::direct_rpa_rhf_energy_kpoints;
    using ATV::cpu::gpw::internal::mo_eri_ovov;
    using ATV::cpu::gpw::internal::Cplx;

    const double cell[3] = {6.0, 6.0, 6.0};
    const UniformGrid grid = make_uniform_grid(cell, 30.0);
    const long npts = grid.npts();
    const int no = 1, nv = 2;
    std::vector<double> psi_o(static_cast<std::size_t>(no) * npts);
    std::vector<double> psi_v(static_cast<std::size_t>(nv) * npts);
    for (long g = 0; g < npts; ++g) {
        double r[3];
        grid.point(g, r);
        const double dx = r[0] - 3.0, dy = r[1] - 3.0, dz = r[2] - 3.0;
        const double r2 = dx * dx + dy * dy + dz * dz;
        psi_o[static_cast<std::size_t>(g)] = std::exp(-0.7 * r2);
        psi_v[static_cast<std::size_t>(0) * npts + g] = dx * std::exp(-0.5 * r2);
        psi_v[static_cast<std::size_t>(1) * npts + g] = dy * std::exp(-0.45 * r2);
    }
    // Same Poisson route as multi-k (not G-space RI) so the integral tensors
    // are bit-identical for real MOs.
    const std::vector<double> v =
        mo_eri_ovov(grid, psi_o, no, psi_v, nv, /*nonperiodic=*/false);
    std::vector<double> eps_o = {-0.85};
    std::vector<double> eps_v = {0.25, 0.70};
    double e_real = 0.0, tr_real = 0.0;
    std::vector<double> om_real;
    direct_rpa_rhf_energy(v, eps_o, eps_v, e_real, tr_real, om_real);

    std::vector<std::vector<Cplx>> psi_ok(1), psi_vk(1);
    psi_ok[0].resize(static_cast<std::size_t>(no) * npts);
    psi_vk[0].resize(static_cast<std::size_t>(nv) * npts);
    for (std::size_t i = 0; i < psi_o.size(); ++i)
        psi_ok[0][i] = Cplx{psi_o[i], 0.0};
    for (std::size_t i = 0; i < psi_v.size(); ++i)
        psi_vk[0][i] = Cplx{psi_v[i], 0.0};
    std::vector<std::vector<double>> eps_ok = {eps_o};
    std::vector<std::vector<double>> eps_vk = {eps_v};
    std::vector<std::array<double, 3>> k_frac = {{{0.0, 0.0, 0.0}}};
    std::vector<double> k_w = {1.0};
    double e_k = 0.0, tr_k = 0.0;
    std::vector<double> om_k;
    direct_rpa_rhf_energy_kpoints(grid, psi_ok, psi_vk, eps_ok, eps_vk, k_frac,
                                  k_w, e_k, tr_k, om_k);
    EXPECT_NEAR(e_k, e_real, 1e-10);
    EXPECT_NEAR(tr_k, tr_real, 1e-10);
    ASSERT_EQ(om_k.size(), om_real.size());
    for (std::size_t s = 0; s < om_k.size(); ++s)
        EXPECT_NEAR(om_k[s], om_real[s], 1e-10);
}
