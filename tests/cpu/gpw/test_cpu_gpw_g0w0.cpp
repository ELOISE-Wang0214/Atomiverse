/// @file tests/cpu/gpw/test_cpu_gpw_g0w0.cpp
/// @brief GPW G0W0 kernel units: dRPA amplitudes match energy path +
///        X/Y norm; GpwDft::g0w0 guard throws. No SCF (unit tier).

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/mp2_kernels.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::internal::DirectRpaAmplitudes;
using ATV::cpu::gpw::internal::direct_rpa_rhf;
using ATV::cpu::gpw::internal::direct_rpa_rhf_energy;

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

TEST(GpwG0W0, AmplitudesMatchEnergyPath) {
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

    DirectRpaAmplitudes amp;
    direct_rpa_rhf(v, eps_occ, eps_virt, amp);

    EXPECT_NEAR(amp.e_corr, e_corr, 1e-12);
    EXPECT_NEAR(amp.tr_A, tr_A, 1e-12);
    ASSERT_EQ(amp.omega.size(), omega.size());
    for (std::size_t k = 0; k < omega.size(); ++k)
        EXPECT_NEAR(amp.omega[k], omega[k], 1e-12);

    // XᵀX − YᵀY ≈ I for each state (RPA metric norm).
    for (int st = 0; st < nov; ++st) {
        double nrm = 0.0;
        for (int ia = 0; ia < nov; ++ia) {
            const double Xs =
                amp.X[static_cast<std::size_t>(st) * nov
                      + static_cast<std::size_t>(ia)];
            const double Ys =
                amp.Y[static_cast<std::size_t>(st) * nov
                      + static_cast<std::size_t>(ia)];
            nrm += Xs * Xs - Ys * Ys;
        }
        EXPECT_NEAR(nrm, 1.0, 1e-8);
    }
}

TEST(GpwG0W0, G0W0BeforeBuildThrows) {
    GpwDft m;
    EXPECT_THROW(m.g0w0({}), std::logic_error);
}

TEST(GpwG0W0, G0W0UksThrows) {
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 6.0;
    m.n_alpha = 1;
    m.n_beta = 1;
    GpwAtom a;
    a.z_val = 2.0;
    a.basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    m.build({a});
    EXPECT_THROW(m.g0w0({a}), std::invalid_argument);
}

TEST(GpwG0W0, G0W0UnconvergedThrows) {
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 6.0;
    m.n_occ = 1;
    GpwAtom a;
    a.z_val = 2.0;
    a.basis.push_back(s_gauss(0.8, 3.0, 3.0, 3.0));
    m.build({a});
    EXPECT_THROW(m.g0w0({a}), std::logic_error);
}

// Multi-k kernel at N_k=1 with real-valued complex MOs recovers the real
// Γ spectral HOMO/LUMO QP shifts on the same orbitals (no SCF).
TEST(GpwG0W0, KpointsNk1RecoversRealSpectral) {
    using ATV::cpu::gpw::make_uniform_grid;
    using ATV::cpu::gpw::UniformGrid;
    using ATV::cpu::gpw::internal::Cplx;
    using ATV::cpu::gpw::internal::direct_rpa_rhf;
    using ATV::cpu::gpw::internal::DirectRpaAmplitudes;
    using ATV::cpu::gpw::internal::exchange_diag_mt;
    using ATV::cpu::gpw::internal::mo_eri_ovov;
    using ATV::cpu::gpw::internal::mo_pair_potentials_ov;
    using ATV::cpu::gpw::internal::spectral_g0w0_rhf_kpoints;

    const double cell[3] = {6.0, 6.0, 6.0};
    const UniformGrid grid = make_uniform_grid(cell, 25.0);
    const long npts = grid.npts();
    const double gw = grid.weight();
    const int n_occ = 1, n_virt = 1, n_mo = 2;

    std::vector<double> psi_all(static_cast<std::size_t>(n_mo) * npts);
    for (long g = 0; g < npts; ++g) {
        double r[3];
        grid.point(g, r);
        const double dx = r[0] - 3.0, dy = r[1] - 3.0, dz = r[2] - 3.0;
        const double r2 = dx * dx + dy * dy + dz * dz;
        psi_all[static_cast<std::size_t>(g)] = std::exp(-0.7 * r2);
        psi_all[static_cast<std::size_t>(npts) + static_cast<std::size_t>(g)] =
            dx * std::exp(-0.5 * r2);
    }
    // Normalize loosely not required — spectral formula is scale-sensitive
    // but both paths use the same raw values.
    std::vector<double> psi_o(static_cast<std::size_t>(n_occ) * npts);
    std::vector<double> psi_v(static_cast<std::size_t>(n_virt) * npts);
    for (long g = 0; g < npts; ++g) {
        psi_o[static_cast<std::size_t>(g)] = psi_all[static_cast<std::size_t>(g)];
        psi_v[static_cast<std::size_t>(g)] =
            psi_all[static_cast<std::size_t>(npts) + static_cast<std::size_t>(g)];
    }
    std::vector<double> eps = {-0.85, 0.30};
    std::vector<double> eps_o = {eps[0]};
    std::vector<double> eps_v = {eps[1]};

    const std::vector<double> v =
        mo_eri_ovov(grid, psi_o, n_occ, psi_v, n_virt, /*nonperiodic=*/false);
    DirectRpaAmplitudes rpa;
    direct_rpa_rhf(v, eps_o, eps_v, rpa);
    const int nov = n_occ * n_virt;
    const std::vector<double> V_ov =
        mo_pair_potentials_ov(grid, psi_o, n_occ, psi_v, n_virt);
    const double sqrt2 = std::sqrt(2.0);
    std::vector<double> tdm(static_cast<std::size_t>(nov)
                                * static_cast<std::size_t>(n_mo)
                                * static_cast<std::size_t>(n_mo),
                            0.0);
    std::vector<double> Ws(static_cast<std::size_t>(npts));
    for (int s = 0; s < nov; ++s) {
        std::fill(Ws.begin(), Ws.end(), 0.0);
        for (int ia = 0; ia < nov; ++ia) {
            const double z =
                sqrt2
                * (rpa.X[static_cast<std::size_t>(s) * nov
                         + static_cast<std::size_t>(ia)]
                   + rpa.Y[static_cast<std::size_t>(s) * nov
                           + static_cast<std::size_t>(ia)]);
            const double* Via =
                V_ov.data()
                + static_cast<std::size_t>(ia) * static_cast<std::size_t>(npts);
            for (long g = 0; g < npts; ++g)
                Ws[static_cast<std::size_t>(g)] += z * Via[g];
        }
        const std::size_t s0 = static_cast<std::size_t>(s)
                               * static_cast<std::size_t>(n_mo)
                               * static_cast<std::size_t>(n_mo);
        for (int m = 0; m < n_mo; ++m) {
            const double* pm =
                psi_all.data()
                + static_cast<std::size_t>(m) * static_cast<std::size_t>(npts);
            for (int p = 0; p < n_mo; ++p) {
                const double* pp =
                    psi_all.data()
                    + static_cast<std::size_t>(p)
                          * static_cast<std::size_t>(npts);
                double acc = 0.0;
                for (long g = 0; g < npts; ++g)
                    acc += Ws[static_cast<std::size_t>(g)] * pm[g] * pp[g];
                tdm[s0 + static_cast<std::size_t>(m)
                             * static_cast<std::size_t>(n_mo)
                         + static_cast<std::size_t>(p)] = gw * acc;
            }
        }
    }
    const std::vector<double> sx = exchange_diag_mt(grid, psi_all, n_mo, n_occ);
    // Flat zero vxc for the unit recovery (identical on both paths).
    std::vector<double> vxc_r(static_cast<std::size_t>(npts), 0.0);
    std::vector<double> vxc_diag(static_cast<std::size_t>(n_mo), 0.0);

    // Real Newton (inline: residual zero at ω = ε + Σx when Σc=0; not
    // zero in general — use a few Newton steps matching g0w0.cpp).
    auto sigma_c = [&](double omega, int p) {
        double re = 0.0;
        for (int s = 0; s < nov; ++s) {
            const double Os = rpa.omega[static_cast<std::size_t>(s)];
            const std::size_t s0 =
                static_cast<std::size_t>(s) * n_mo * n_mo
                + static_cast<std::size_t>(p);
            for (int m = 0; m < n_occ; ++m) {
                const double t =
                    tdm[s0 + static_cast<std::size_t>(m) * n_mo];
                const double den_r = omega + Os - eps[static_cast<std::size_t>(m)];
                re += (t * t) * den_r / (den_r * den_r + 1e-16);
            }
            for (int m = n_occ; m < n_mo; ++m) {
                const double t =
                    tdm[s0 + static_cast<std::size_t>(m) * n_mo];
                const double den_r = omega - Os - eps[static_cast<std::size_t>(m)];
                re += (t * t) * den_r / (den_r * den_r + 1e-16);
            }
        }
        return re;
    };
    std::vector<double> qp_real(static_cast<std::size_t>(n_mo));
    for (int p = 0; p < n_mo; ++p) {
        double w = eps[static_cast<std::size_t>(p)];
        const double vx = sx[static_cast<std::size_t>(p)]
                          - vxc_diag[static_cast<std::size_t>(p)];
        for (int it = 0; it < 50; ++it) {
            const double sig = sigma_c(w, p);
            const double f = w - eps[static_cast<std::size_t>(p)] - (sig + vx);
            if (std::abs(f) < 1e-10)
                break;
            const double sig_p = sigma_c(w + 1e-6, p);
            w -= f / (1.0 - (sig_p - sig) / 1e-6);
        }
        qp_real[static_cast<std::size_t>(p)] = w;
    }

    std::vector<std::vector<Cplx>> psi_k(1);
    psi_k[0].resize(static_cast<std::size_t>(n_mo) * npts);
    for (std::size_t i = 0; i < psi_all.size(); ++i)
        psi_k[0][i] = Cplx{psi_all[i], 0.0};
    std::vector<std::vector<double>> eps_k = {eps};
    std::vector<double> k_w = {1.0};
    std::vector<std::vector<double>> qp_k, sx_k, vxc_k;
    bool ok = false;
    spectral_g0w0_rhf_kpoints(grid, psi_k, eps_k, k_w, n_occ, n_virt, vxc_r,
                              qp_k, sx_k, vxc_k, ok);
    ASSERT_EQ(qp_k.size(), 1u);
    ASSERT_EQ(qp_k[0].size(), static_cast<std::size_t>(n_mo));
    EXPECT_TRUE(ok);
    for (int p = 0; p < n_mo; ++p) {
        EXPECT_NEAR(qp_k[0][static_cast<std::size_t>(p)],
                    qp_real[static_cast<std::size_t>(p)], 1e-8);
        EXPECT_NEAR(sx_k[0][static_cast<std::size_t>(p)],
                    sx[static_cast<std::size_t>(p)], 1e-8);
    }
}
