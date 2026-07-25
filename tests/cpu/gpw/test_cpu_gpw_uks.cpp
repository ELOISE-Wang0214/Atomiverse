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

// ── Multi-k collinear UKS (gpw-kpoints-uks) ────────────────────────────────

TEST(GpwUks, KpointsUksAcceptedAtBuild) {
    auto atoms = small_h2();
    GpwDft m;
    m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
    m.n_alpha = 2;
    m.n_beta = 0;
    m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
    m.xc = XcKind::PADE;
    EXPECT_NO_THROW(m.build(atoms));
}

TEST(GpwUks, KpointsUksGammaLimitMatchesGammaPath) {
    // Explicit k = Γ through the complex k-UKS machinery must recover the
    // legacy real Γ-UKS energy on the same toy H₂ high-spin setup.
    auto atoms = small_h2();
    GpwDft gamma;
    gamma.cell[0] = gamma.cell[1] = gamma.cell[2] = 8.0;
    gamma.cutoff_ry = 100.0;
    gamma.n_alpha = 2;
    gamma.n_beta = 0;
    gamma.max_iter = 80;
    gamma.tol = 1e-7;
    gamma.mix = 0.5;
    gamma.xc = XcKind::PADE;
    gamma.build(atoms);
    const GpwResult r_g = gamma.run(atoms);
    ASSERT_TRUE(r_g.converged) << r_g.e_total;

    GpwDft kpath;
    kpath.cell[0] = kpath.cell[1] = kpath.cell[2] = 8.0;
    kpath.cutoff_ry = 100.0;
    kpath.n_alpha = 2;
    kpath.n_beta = 0;
    kpath.max_iter = 80;
    kpath.tol = 1e-7;
    kpath.mix = 0.5;
    kpath.xc = XcKind::PADE;
    kpath.k_points = {{0.0, 0.0, 0.0}};
    kpath.k_point_weights = {1.0};
    kpath.build(atoms);
    const GpwResult r_k = kpath.run(atoms);
    ASSERT_TRUE(r_k.converged) << r_k.e_total;
    EXPECT_NEAR(r_k.e_total, r_g.e_total, 1e-6);
    EXPECT_EQ(r_k.n_kpts, 1);
    ASSERT_EQ(r_k.eps_k.size(), 1u);
}

TEST(GpwUks, KpointsUksClosedShellMatchesRks) {
    // n_alpha = n_beta = 1 and closed-shell n_occ = 1 on the same MP mesh
    // must agree (spin-unpolarized density is the RKS solution).
    auto atoms = small_h2();
    GpwDft rks;
    rks.cell[0] = rks.cell[1] = rks.cell[2] = 8.0;
    rks.cutoff_ry = 100.0;
    rks.n_occ = 1;
    rks.max_iter = 80;
    rks.tol = 1e-7;
    rks.xc = XcKind::PADE;
    rks.k_mesh[0] = rks.k_mesh[1] = 1;
    rks.k_mesh[2] = 2;
    rks.k_time_reversal = true;
    rks.build(atoms);
    const GpwResult rr = rks.run(atoms);
    ASSERT_TRUE(rr.converged) << rr.e_total;

    GpwDft uks;
    uks.cell[0] = uks.cell[1] = uks.cell[2] = 8.0;
    uks.cutoff_ry = 100.0;
    uks.n_alpha = 1;
    uks.n_beta = 1;
    uks.max_iter = 80;
    uks.tol = 1e-7;
    uks.xc = XcKind::PADE;
    uks.k_mesh[0] = uks.k_mesh[1] = 1;
    uks.k_mesh[2] = 2;
    uks.k_time_reversal = true;
    uks.build(atoms);
    const GpwResult ru = uks.run(atoms);
    ASSERT_TRUE(ru.converged) << ru.e_total;
    EXPECT_NEAR(ru.e_total, rr.e_total, 5e-6);
}

TEST(GpwUks, KpointsUksTimeReversalIdentity) {
    auto atoms = small_h2();
    auto run_mesh = [&](bool tr) {
        GpwDft m;
        m.cell[0] = m.cell[1] = m.cell[2] = 8.0;
        m.cutoff_ry = 100.0;
        m.n_alpha = 2;
        m.n_beta = 0;
        m.max_iter = 80;
        m.tol = 1e-7;
        m.xc = XcKind::PADE;
        m.k_mesh[0] = m.k_mesh[1] = m.k_mesh[2] = 2;
        m.k_time_reversal = tr;
        m.build(atoms);
        return m.run(atoms);
    };
    const GpwResult full = run_mesh(false);
    const GpwResult fold = run_mesh(true);
    ASSERT_TRUE(full.converged);
    ASSERT_TRUE(fold.converged);
    EXPECT_EQ(full.n_kpts, 8);
    EXPECT_EQ(fold.n_kpts, 4);
    EXPECT_NEAR(full.e_total, fold.e_total, 5e-6);
}
