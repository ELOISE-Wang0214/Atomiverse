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
