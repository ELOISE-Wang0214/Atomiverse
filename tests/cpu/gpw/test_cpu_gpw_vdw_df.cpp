/// @file tests/cpu/gpw/test_cpu_gpw_vdw_df.cpp
/// @brief Unit tests for nonlocal vdW-DF (DRSLL) / vdW-DF2 (LMKLL) vs libvdwxc.

#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/vdw_df.h"

using ATV::cpu::gpw::kVdwDf2Zab;
using ATV::cpu::gpw::UniformGrid;
using ATV::cpu::gpw::vdw_df_drsll;
using ATV::cpu::gpw::vdw_df_lmkll;
using ATV::cpu::gpw::vdw_df_q0;
using ATV::cpu::gpw::VdwDfResult;

namespace {

/// @brief Synthetic Gaussian blob density used for the offline libvdwxc pin.
void gaussian_blob(const UniformGrid& grid, double a, double A,
                   std::vector<double>& rho, std::vector<double>& sigma) {
    const long ntot = grid.npts();
    rho.assign(static_cast<std::size_t>(ntot), 0.0);
    sigma.assign(static_cast<std::size_t>(ntot), 0.0);
    const double cx = 0.5 * grid.L[0];
    const double cy = 0.5 * grid.L[1];
    const double cz = 0.5 * grid.L[2];
    for (long g = 0; g < ntot; ++g) {
        double r[3];
        grid.point(g, r);
        double dx = r[0] - cx;
        double dy = r[1] - cy;
        double dz = r[2] - cz;
        auto mi = [](double d, double L) {
            d -= L * std::nearbyint(d / L);
            return d;
        };
        dx = mi(dx, grid.L[0]);
        dy = mi(dy, grid.L[1]);
        dz = mi(dz, grid.L[2]);
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double e = std::exp(-a * r2);
        rho[static_cast<std::size_t>(g)] = A * e;
        const double fac = -2.0 * a * A * e;
        const double gx = fac * dx;
        const double gy = fac * dy;
        const double gz = fac * dz;
        sigma[static_cast<std::size_t>(g)] = gx * gx + gy * gy + gz * gz;
    }
}

} // namespace

// Offline pin: libvdwxc 0.5.0 VDWXC_DF1 on 16³ cubic L=8 Bohr Gaussian blob
// (a=0.5, A=1.0, analytic σ). See tools/cp2k_oracle/inputs/gpw-vdw-df-ice/.
TEST(GpwVdwDf, GaussianBlobMatchesLibvdwxc) {
    UniformGrid grid;
    grid.L[0] = grid.L[1] = grid.L[2] = 8.0;
    grid.n[0] = grid.n[1] = grid.n[2] = 16;
    grid.origin[0] = grid.origin[1] = grid.origin[2] = 0.0;

    std::vector<double> rho, sigma;
    gaussian_blob(grid, /*a=*/0.5, /*A=*/1.0, rho, sigma);

    const VdwDfResult r = vdw_df_drsll(grid, rho, sigma);
    const double libvdwxc = 1.2367739353586624e-01; // Ha
    std::printf("[vdW-DF] E_nl=%.16e  libvdwxc=%.16e  diff=%.3e Ha\n", r.energy,
                libvdwxc, r.energy - libvdwxc);
    // float32 kernel tables → expect ~1e-7 absolute; pin at 1e-6 Ha.
    EXPECT_NEAR(r.energy, libvdwxc, 1e-6);
}

TEST(GpwVdwDf, Q0LowDensityIsQCut) {
    std::vector<double> rho = {1e-20};
    std::vector<double> sigma = {0.0};
    std::vector<double> q0;
    vdw_df_q0(rho, sigma, q0);
    EXPECT_NEAR(q0[0], ATV::cpu::gpw::kVdwDfQCut, 1e-15);
}

TEST(GpwVdwDf, Q0FiniteDensityInRange) {
    // Uniform density ρ = 0.1 → q0 should be between 0 and q_cut.
    std::vector<double> rho = {0.1};
    std::vector<double> sigma = {0.0};
    std::vector<double> q0;
    vdw_df_q0(rho, sigma, q0);
    EXPECT_GT(q0[0], 0.0);
    EXPECT_LE(q0[0], ATV::cpu::gpw::kVdwDfQCut + 1e-12);
}

// Offline pin: libvdwxc 0.5.0 VDWXC_DF2 on the same 16³ Gaussian blob.
// See tools/cp2k_oracle/inputs/gpw-vdw-df2-ice/.
TEST(GpwVdwDf, GaussianBlobDf2MatchesLibvdwxc) {
    UniformGrid grid;
    grid.L[0] = grid.L[1] = grid.L[2] = 8.0;
    grid.n[0] = grid.n[1] = grid.n[2] = 16;
    grid.origin[0] = grid.origin[1] = grid.origin[2] = 0.0;

    std::vector<double> rho, sigma;
    gaussian_blob(grid, /*a=*/0.5, /*A=*/1.0, rho, sigma);

    const VdwDfResult r = vdw_df_lmkll(grid, rho, sigma);
    const double libvdwxc = 1.1428161716669384e-01; // Ha
    std::printf("[vdW-DF2] E_nl=%.16e  libvdwxc=%.16e  diff=%.3e Ha\n", r.energy,
                libvdwxc, r.energy - libvdwxc);
    EXPECT_NEAR(r.energy, libvdwxc, 1e-6);
}

TEST(GpwVdwDf, Df2Q0DiffersFromDf1WhenGradient) {
    // With nonzero σ, DF2 Z_ab = −1.887 yields a different q0 than DF1.
    std::vector<double> rho = {0.1};
    std::vector<double> sigma = {0.01};
    std::vector<double> q0_df1, q0_df2;
    vdw_df_q0(rho, sigma, q0_df1);
    vdw_df_q0(rho, sigma, q0_df2, kVdwDf2Zab);
    EXPECT_NE(q0_df1[0], q0_df2[0]);
    EXPECT_GT(q0_df1[0], 0.0);
    EXPECT_GT(q0_df2[0], 0.0);
}
