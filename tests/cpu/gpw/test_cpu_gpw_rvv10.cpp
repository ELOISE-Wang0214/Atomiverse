/// @file tests/cpu/gpw/test_cpu_gpw_rvv10.cpp
/// @brief Unit tests for nonlocal rVV10 RPS kernel (offline CP2K-aligned tables).

#include <cmath>
#include <cstdio>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/rvv10.h"

using ATV::cpu::gpw::kRvv10B;
using ATV::cpu::gpw::kRvv10C;
using ATV::cpu::gpw::kRvv10QCut;
using ATV::cpu::gpw::kRvv10QMin;
using ATV::cpu::gpw::rvv10_energy;
using ATV::cpu::gpw::rvv10_q0;
using ATV::cpu::gpw::Rvv10Result;
using ATV::cpu::gpw::UniformGrid;

namespace {

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

TEST(GpwRvv10, Q0LowDensityIsQCut) {
    std::vector<double> rho = {1e-20};
    std::vector<double> sigma = {0.0};
    std::vector<double> q0;
    rvv10_q0(rho, sigma, q0);
    EXPECT_NEAR(q0[0], kRvv10QCut, 1e-15);
}

TEST(GpwRvv10, Q0FiniteDensityInMesh) {
    std::vector<double> rho = {0.1};
    std::vector<double> sigma = {0.0};
    std::vector<double> q0;
    rvv10_q0(rho, sigma, q0);
    EXPECT_GE(q0[0], kRvv10QMin - 1e-15);
    EXPECT_LE(q0[0], kRvv10QCut + 1e-12);
}

TEST(GpwRvv10, BetaTermPositiveOnUniformBlob) {
    // β∫ρ must be positive for any non-zero density (β > 0 for b = 6.3).
    UniformGrid grid;
    grid.L[0] = grid.L[1] = grid.L[2] = 8.0;
    grid.n[0] = grid.n[1] = grid.n[2] = 16;
    grid.origin[0] = grid.origin[1] = grid.origin[2] = 0.0;
    std::vector<double> rho, sigma;
    gaussian_blob(grid, /*a=*/0.5, /*A=*/1.0, rho, sigma);
    const Rvv10Result r = rvv10_energy(grid, rho, sigma);
    std::printf("[rVV10] E=%.16e  e_kernel=%.16e  e_beta=%.16e Ha\n", r.energy,
                r.e_kernel, r.e_beta);
    EXPECT_TRUE(std::isfinite(r.energy));
    EXPECT_GT(r.e_beta, 0.0);
    // Dual-frozen synthetic pin (16³, L=8 Bohr, a=0.5, A=1.0, b=6.3, C=0.0093).
    // Printed on first green run; float32 kernel → 1e-5 Ha gate.
    constexpr double kPin = 5.4275454996569092e-02; // Ha (dual-frozen)
    EXPECT_NEAR(r.energy, kPin, 1e-5);
}

TEST(GpwRvv10, DefaultParamsMatchSabatini) {
    EXPECT_NEAR(kRvv10B, 6.3, 1e-15);
    EXPECT_NEAR(kRvv10C, 0.0093, 1e-15);
}
