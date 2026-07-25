/// @file tests/cpu/gpw/test_cpu_gpw_uks.cpp
/// @brief Unit tests for Γ-point UKS GPW (PADE LSD + PBE GGA, dual H).
///
/// Fast gate only: tiny H₂ triplet-like (n_alpha=2, n_beta=0) and closed-shell
/// still works. Full O₂ CP2K parity lives in regressions/gpw/gpw-uks-o2-gamma
/// (PADE) and regressions/gpw/gpw-uks-pbe (PBE).

#include <cmath>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/internal/xc_grid.h"
#include "atomiverse/cpu/gpw/internal/xc_grid_gga.h"

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::UniformGrid;
using ATV::cpu::gpw::xc_grid_pade;
using ATV::cpu::gpw::xc_grid_pade_lsd;
using ATV::cpu::gpw::xc_grid_pbe;
using ATV::cpu::gpw::xc_grid_pbe_spin;
using ATV::cpu::gpw::XcKind;

namespace {

GridBasisFunction s_basis(double alpha, const double c[3]) {
    GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = c[0];
    p.center[1] = c[1];
    p.center[2] = c[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = std::pow(2.0 * alpha / M_PI, 0.75);
    return {p};
}

GpwAtom hydrogen(const double pos[3]) {
    GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis(0.5, pos)};
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

std::vector<GpwAtom> small_h2(double box = 8.0) {
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    return {hydrogen(rA), hydrogen(rB)};
}

} // namespace

TEST(GpwUks, PadeLsdMatchesUnpolarizedCollapse) {
    // ρ_α = ρ_β = ρ/2 must recover the spin-unpolarized PADE energy/potential.
    UniformGrid grid;
    grid.n[0] = grid.n[1] = grid.n[2] = 4;
    grid.L[0] = grid.L[1] = grid.L[2] = 2.0;
    const long ng = grid.npts();
    std::vector<double> rho(static_cast<std::size_t>(ng), 0.1);
    std::vector<double> ra(static_cast<std::size_t>(ng), 0.05);
    std::vector<double> rb(static_cast<std::size_t>(ng), 0.05);
    const auto uni = xc_grid_pade(grid, rho);
    const auto lsd = xc_grid_pade_lsd(grid, ra, rb);
    EXPECT_NEAR(uni.e_xc, lsd.e_xc, 1e-12);
    for (long g = 0; g < ng; ++g) {
        EXPECT_NEAR(uni.vxc_r[static_cast<std::size_t>(g)],
                    lsd.vxc_a[static_cast<std::size_t>(g)], 1e-12);
        EXPECT_NEAR(lsd.vxc_a[static_cast<std::size_t>(g)],
                    lsd.vxc_b[static_cast<std::size_t>(g)], 1e-12);
    }
}

TEST(GpwUks, ClosedShellStillConverges) {
    auto atoms = small_h2();
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
    m.cutoff_ry = 120.0;
    m.n_occ = 1;
    m.max_iter = 80;
    m.tol = 1e-7;
    m.mix = 0.5;
    m.build(atoms);
    const GpwResult r = m.run(atoms);
    EXPECT_TRUE(r.converged);
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_NEAR(r.s2, 0.0, 1e-12);
}

TEST(GpwUks, TripletLikeH2Converges) {
    // Two α electrons, zero β — high-spin limit on the H₂ toy (not physical
    // chemistry, but exercises dual H / S² / occupations).
    auto atoms = small_h2();
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
    m.cutoff_ry = 120.0;
    m.n_alpha = 2;
    m.n_beta = 0;
    m.n_occ = 0;
    m.max_iter = 100;
    m.tol = 1e-6;
    m.mix = 0.4;
    m.xc = XcKind::PADE;
    m.build(atoms);
    ASSERT_TRUE(m.is_uks());
    const GpwResult r = m.run(atoms);
    EXPECT_TRUE(r.converged) << "iters=" << r.iters << " E=" << r.e_total;
    EXPECT_TRUE(std::isfinite(r.e_total));
    // Ideal S=1 → ⟨S²⟩ = 2 for pure high-spin (Nβ=0).
    EXPECT_NEAR(r.s2, 2.0, 1e-8);
    EXPECT_GT(static_cast<int>(r.C.shape(0)), 0);
}

TEST(GpwUks, PbeSpinMatchesUnpolarizedCollapse) {
    // ρ_α = ρ_β = ρ/2 must recover closed-shell PBE energy/potential.
    UniformGrid grid;
    grid.n[0] = grid.n[1] = grid.n[2] = 6;
    grid.L[0] = grid.L[1] = grid.L[2] = 4.0;
    const long ng = grid.npts();
    std::vector<double> rho(static_cast<std::size_t>(ng), 0.12);
    std::vector<double> ra(static_cast<std::size_t>(ng), 0.06);
    std::vector<double> rb(static_cast<std::size_t>(ng), 0.06);
    const auto uni = xc_grid_pbe(grid, rho);
    const auto spin = xc_grid_pbe_spin(grid, ra, rb);
    EXPECT_NEAR(uni.e_xc, spin.e_xc, 1e-10 * std::max(1.0, std::abs(uni.e_xc)));
    for (long g = 0; g < ng; ++g) {
        EXPECT_NEAR(uni.vxc_r[static_cast<std::size_t>(g)],
                    spin.vxc_a[static_cast<std::size_t>(g)], 1e-9);
        EXPECT_NEAR(spin.vxc_a[static_cast<std::size_t>(g)],
                    spin.vxc_b[static_cast<std::size_t>(g)], 1e-12);
    }
}

TEST(GpwUks, PbeUksAcceptedAtBuild) {
    auto atoms = small_h2();
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
    m.n_alpha = 1;
    m.n_beta = 1;
    m.xc = XcKind::PBE;
    EXPECT_NO_THROW(m.build(atoms));
    ASSERT_TRUE(m.is_uks());
}

TEST(GpwUks, TripletLikeH2ConvergesPbe) {
    auto atoms = small_h2();
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
    m.cutoff_ry = 120.0;
    m.n_alpha = 2;
    m.n_beta = 0;
    m.n_occ = 0;
    m.max_iter = 100;
    m.tol = 1e-6;
    m.mix = 0.4;
    m.xc = XcKind::PBE;
    m.build(atoms);
    const GpwResult r = m.run(atoms);
    EXPECT_TRUE(r.converged) << "iters=" << r.iters << " E=" << r.e_total;
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_NEAR(r.s2, 2.0, 1e-8);
}
