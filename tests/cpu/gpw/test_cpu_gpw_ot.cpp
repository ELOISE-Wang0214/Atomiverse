/// @file tests/cpu/gpw/test_cpu_gpw_ot.cpp
/// @brief Unit tests for GPW Orbital Transformation (OT) minimizer.
///
/// Single-function units only (no full CP2K e2e — that lives in
/// regressions/gpw/ot-minimizer-water). Covers:
///   1. fixed-H OT band energy vs generalized diagonalization
///   2. mini-H2 SCF: OT vs DIIS total energy parity
///   3. GPWEvaluator::use_ot product switch (no ATV_GPW_OT env required)

#include <cmath>
#include <cstdlib>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gtest/gtest.h>
#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/gpw/internal/kpoints.h"
#include "atomiverse/cpu/gpw/internal/ot.h"
#include "atomiverse/cpu/linalg/eigensolve.h"
#include "atomiverse/cpu/linalg/matrix.h"

#include <complex>
#include <stdexcept>
#include <string>

using ATV::cpu::gpw::GpwAtom;
using ATV::cpu::gpw::GpwDft;
using ATV::cpu::gpw::GpwResult;
using ATV::cpu::gpw::GridBasisFunction;
using ATV::cpu::gpw::GridPrimitive;
using ATV::cpu::gpw::internal::ot_minimize_fixed;
using ATV::cpu::linalg::diagonalize_fock_generalized;
using ATV::cpu::linalg::RealMatrix;

namespace {

void pin_single_thread() {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
}

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

GpwDft small_cfg(double box = 8.0) {
    GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = box;
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 1;
    cfg.max_iter = 100;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;
    return cfg;
}

} // namespace

// Fixed-H OT must recover the occupied band energy of the generalized
// eigenproblem (same subspace as diagonalization).
TEST(GpwOt, FixedHMatchesDiagonalizationBandEnergy) {
    pin_single_thread();
    constexpr int nb = 3;
    constexpr int no = 1;
    RealMatrix H = xt::zeros<double>({nb, nb});
    RealMatrix S = xt::zeros<double>({nb, nb});
    RealMatrix T = xt::zeros<double>({nb, nb});
    // Symmetric H with a clear gap; S = I; T = I (FULL_KINETIC precond).
    H(0, 0) = -1.5;
    H(1, 1) = -0.4;
    H(2, 2) = 0.8;
    H(0, 1) = H(1, 0) = 0.12;
    H(0, 2) = H(2, 0) = 0.05;
    H(1, 2) = H(2, 1) = -0.03;
    for (int i = 0; i < nb; ++i) {
        S(i, i) = 1.0;
        T(i, i) = 1.0;
    }

    auto [eig, C] = diagonalize_fock_generalized(H, S);
    RealMatrix C0 = xt::zeros<double>({nb, no});
    for (int u = 0; u < nb; ++u)
        C0(u, 0) = C(u, 0); // lowest eigenvector as guess (already optimal)

    // Perturb guess so OT has work to do.
    RealMatrix Cguess = C0;
    Cguess(1, 0) += 0.3;
    Cguess(2, 0) -= 0.2;

    const auto ot = ot_minimize_fixed(H, S, T, Cguess, /*tol=*/1e-8,
                                      /*max_iter=*/200);
    ASSERT_TRUE(ot.converged)
        << "OT fixed-H did not converge in " << ot.iters << " iters";

    const double e_diag = 2.0 * eig(0); // closed-shell band energy
    EXPECT_NEAR(ot.band_energy, e_diag, 1e-6)
        << "OT band=" << ot.band_energy << " diag=" << e_diag;
}

// Mini-H2 SCF: OT and DIIS must agree on the same gapped fixed point.
TEST(GpwOt, MiniH2OtMatchesDiisEnergy) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft diis = small_cfg();
    diis.use_ot = false;
    diis.build(atoms);
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged)
        << "DIIS failed in " << r_diis.iters << " iters";

    GpwDft ot = small_cfg();
    ot.use_ot = true;
    ot.max_iter = 200;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged) << "OT failed in " << r_ot.iters << " iters";

    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 1e-5)
        << "OT=" << r_ot.e_total << " DIIS=" << r_diis.e_total
        << " iters ot/diis=" << r_ot.iters << "/" << r_diis.iters;
}

// High-spin H₂ toy UKS: OT and DIIS must agree (occ=1 open-shell gradient).
// Full O₂ CP2K pin lives in regressions/gpw/ot-uks-o2.
TEST(GpwOt, UksOtMatchesDiisEnergy) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft diis = small_cfg();
    diis.n_occ = 0;
    diis.n_alpha = 2;
    diis.n_beta = 0;
    diis.use_ot = false;
    diis.max_iter = 120;
    diis.tol = 1e-6;
    diis.mix = 0.4;
    diis.xc = ATV::cpu::gpw::XcKind::PADE;
    diis.build(atoms);
    ASSERT_TRUE(diis.is_uks());
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged)
        << "UKS DIIS failed in " << r_diis.iters << " iters E=" << r_diis.e_total;

    GpwDft ot = small_cfg();
    ot.n_occ = 0;
    ot.n_alpha = 2;
    ot.n_beta = 0;
    ot.use_ot = true;
    ot.max_iter = 200;
    ot.tol = 1e-6;
    ot.xc = ATV::cpu::gpw::XcKind::PADE;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged)
        << "UKS OT failed in " << r_ot.iters << " iters E=" << r_ot.e_total;

    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 1e-5)
        << "OT=" << r_ot.e_total << " DIIS=" << r_diis.e_total
        << " iters ot/diis=" << r_ot.iters << "/" << r_diis.iters;
    EXPECT_NEAR(r_ot.s2, r_diis.s2, 1e-6);
}

// Γ-UKS OT + Dudarev DFT+U: product path (no throw); both OT and DIIS
// converge with a non-trivial E_U. Fixed-point parity on a small high-spin
// H₂ + U toy can land in different basins; the NiO e2e pin
// (regressions/gpw/ot-dftu) is the energy oracle.
TEST(GpwOt, UksOtPlusUConverges) {
    pin_single_thread();
    using ATV::cpu::gpw::GaussianBasis;
    using ATV::cpu::gpw::GthPotential;
    using ATV::cpu::gpw::PlusU;
    using ATV::cpu::gpw::XcKind;
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double L = 8.0;
    const double r0[3] = {0.5 * L - 0.7, 0.5 * L, 0.5 * L};
    const double r1[3] = {0.5 * L + 0.7, 0.5 * L, 0.5 * L};
    std::vector<GpwAtom> atoms{GpwAtom(r0, basis, pot, 1),
                               GpwAtom(r1, basis, pot, 1)};

    auto base = [&]() {
        GpwDft cfg;
        cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = L;
        cfg.cutoff_ry = 100.0;
        cfg.n_alpha = 2;
        cfg.n_beta = 0;
        cfg.max_iter = 80;
        cfg.tol = 1e-6;
        cfg.mix = 0.4;
        cfg.xc = XcKind::PADE;
        cfg.plus_u.push_back(PlusU{molrs::Element::H, /*l=*/0, /*u_ev=*/4.0});
        return cfg;
    };

    GpwDft diis = base();
    diis.use_ot = false;
    diis.build(atoms);
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged)
        << "UKS DIIS+U failed in " << r_diis.iters
        << " iters E=" << r_diis.e_total;
    EXPECT_GT(std::fabs(r_diis.e_dftu), 1e-8);

    GpwDft ot = base();
    ot.use_ot = true;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged) << "UKS OT+U failed in " << r_ot.iters
                                << " iters E=" << r_ot.e_total;
    EXPECT_TRUE(std::isfinite(r_ot.e_total));
    EXPECT_GT(std::fabs(r_ot.e_dftu), 1e-8);
    // Same energy class (mHa-class not required on this pathological toy).
    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 0.5)
        << "OT+U=" << r_ot.e_total << " DIIS+U=" << r_diis.e_total;
}

// Product switch is a public member (default false); ctor sets it without env.
// Full Frame+MOLOPT OT SCF is the e2e regression (≤60 s), not this unit gate.
TEST(GpwOt, EvaluatorUseOtMemberDefaultAndCtor) {
    unsetenv("ATV_GPW_OT");
    const ATV::cpu::gpw::Kind H{
        molrs::Element::H,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};

    ATV::cpu::gpw::GPWEvaluator def({H});
    EXPECT_FALSE(def.use_ot);

    ATV::cpu::gpw::GPWEvaluator on(
        /*kinds*/ {H}, /*cutoff_ry*/ 120.0, /*max_iter*/ 50, /*tol*/ 1e-6,
        /*mix*/ 0.4, ATV::cpu::gpw::XcKind::PADE, /*nonperiodic*/ false,
        /*compute_forces*/ false, /*fd_step*/ 1e-3, /*d3*/ false,
        /*use_ot*/ true);
    EXPECT_TRUE(on.use_ot);
    on.use_ot = false;
    EXPECT_FALSE(on.use_ot);
}

// Fixed complex Hermitian H: OT recovers the occupied band energy of
// diagonalize_k (same subspace as complex generalized diagonalization).
TEST(GpwOt, FixedHKMatchesDiagonalizeKBandEnergy) {
    pin_single_thread();
    using Cplx = std::complex<double>;
    using ATV::cpu::gpw::internal::CplxMatrix;
    using ATV::cpu::gpw::internal::KOrth;
    using ATV::cpu::gpw::internal::diagonalize_k;
    using ATV::cpu::gpw::internal::ot_minimize_fixed_k;

    constexpr int nb = 3;
    constexpr int no = 1;
    CplxMatrix H = xt::zeros<Cplx>({nb, nb});
    CplxMatrix S = xt::zeros<Cplx>({nb, nb});
    CplxMatrix T = xt::zeros<Cplx>({nb, nb});
    // Hermitian H with a clear gap; S = I; T = I.
    H(0, 0) = -1.5;
    H(1, 1) = -0.4;
    H(2, 2) = 0.8;
    H(0, 1) = Cplx(0.10, 0.05);
    H(1, 0) = std::conj(H(0, 1));
    H(0, 2) = Cplx(0.04, -0.02);
    H(2, 0) = std::conj(H(0, 2));
    H(1, 2) = Cplx(-0.03, 0.01);
    H(2, 1) = std::conj(H(1, 2));
    for (int i = 0; i < nb; ++i) {
        S(i, i) = 1.0;
        T(i, i) = 1.0;
    }

    KOrth orth(S, 1e-8);
    auto [eps, C] = diagonalize_k(H, orth);
    CplxMatrix Cguess = xt::zeros<Cplx>({nb, no});
    for (int u = 0; u < nb; ++u)
        Cguess(u, 0) = C(u, 0);
    Cguess(1, 0) += Cplx(0.25, 0.1);
    Cguess(2, 0) -= Cplx(0.15, -0.05);

    const auto ot = ot_minimize_fixed_k(H, S, T, Cguess, /*tol=*/1e-8,
                                        /*max_iter=*/200);
    ASSERT_TRUE(ot.converged) << "complex OT fixed-H iters=" << ot.iters;
    const double e_diag = 2.0 * eps[0];
    EXPECT_NEAR(ot.band_energy, e_diag, 1e-6)
        << "OT band=" << ot.band_energy << " diag=" << e_diag;
}

// Multi-k RKS OT vs DIIS on compact H₂ (explicit non-Γ k).
TEST(GpwOt, MultiKOtMatchesDiisEnergy) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft diis = small_cfg();
    diis.k_points.push_back({0.25, 0.0, 0.0});
    diis.k_point_weights.push_back(1.0);
    diis.use_ot = false;
    diis.max_iter = 80;
    diis.tol = 1e-7;
    diis.build(atoms);
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged) << "multi-k DIIS iters=" << r_diis.iters;

    GpwDft ot = small_cfg();
    ot.k_points.push_back({0.25, 0.0, 0.0});
    ot.k_point_weights.push_back(1.0);
    ot.use_ot = true;
    ot.max_iter = 200;
    ot.tol = 1e-7;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged) << "multi-k OT iters=" << r_ot.iters
                                << " E=" << r_ot.e_total;
    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 1e-5)
        << "OT=" << r_ot.e_total << " DIIS=" << r_diis.e_total;
    EXPECT_EQ(r_ot.n_kpts, 1);
}

// Explicit Γ through the complex OT path matches real Γ OT.
TEST(GpwOt, ExplicitGammaOtMatchesRealGammaOt) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft real_ot = small_cfg();
    real_ot.use_ot = true;
    real_ot.max_iter = 200;
    real_ot.tol = 1e-7;
    real_ot.build(atoms);
    const GpwResult r_real = real_ot.run(atoms);
    ASSERT_TRUE(r_real.converged);

    GpwDft k_ot = small_cfg();
    k_ot.use_ot = true;
    k_ot.k_points.push_back({0.0, 0.0, 0.0});
    k_ot.k_point_weights.push_back(1.0);
    k_ot.max_iter = 200;
    k_ot.tol = 1e-7;
    k_ot.build(atoms);
    const GpwResult r_k = k_ot.run(atoms);
    ASSERT_TRUE(r_k.converged);
    EXPECT_NEAR(r_k.e_total, r_real.e_total, 1e-6)
        << "k-Γ OT=" << r_k.e_total << " real Γ OT=" << r_real.e_total;
}

// Multi-k collinear UKS OT vs DIIS (high-spin H₂, κ = (1/4,0,0)).
// Full e2e identities live in regressions/gpw/ot-uks-kpoints.
TEST(GpwOt, MultiKUksOtMatchesDiisEnergy) {
    pin_single_thread();
    const auto atoms = small_h2();

    GpwDft diis = small_cfg();
    diis.n_occ = 0;
    diis.n_alpha = 2;
    diis.n_beta = 0;
    diis.k_points.push_back({0.25, 0.0, 0.0});
    diis.k_point_weights.push_back(1.0);
    diis.use_ot = false;
    diis.max_iter = 80;
    diis.tol = 1e-7;
    diis.xc = ATV::cpu::gpw::XcKind::PADE;
    diis.build(atoms);
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged)
        << "multi-k UKS DIIS iters=" << r_diis.iters;

    GpwDft ot = small_cfg();
    ot.n_occ = 0;
    ot.n_alpha = 2;
    ot.n_beta = 0;
    ot.k_points.push_back({0.25, 0.0, 0.0});
    ot.k_point_weights.push_back(1.0);
    ot.use_ot = true;
    ot.max_iter = 200;
    ot.tol = 1e-7;
    ot.xc = ATV::cpu::gpw::XcKind::PADE;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged) << "multi-k UKS OT iters=" << r_ot.iters
                                << " E=" << r_ot.e_total;
    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 1e-5)
        << "OT=" << r_ot.e_total << " DIIS=" << r_diis.e_total
        << " iters ot/diis=" << r_ot.iters << "/" << r_diis.iters;
    EXPECT_EQ(r_ot.n_kpts, 1);
    EXPECT_NEAR(r_ot.s2, r_diis.s2, 1e-6);
}

// Multi-k RKS OT + Dudarev DFT+U: product path; OT ≡ DIIS fixed point with
// non-trivial e_dftu. Full e2e pin: regressions/gpw/ot-dftu-kpoints.
TEST(GpwOt, MultiKOtPlusUMatchesDiisEnergy) {
    pin_single_thread();
    using ATV::cpu::gpw::GaussianBasis;
    using ATV::cpu::gpw::GthPotential;
    using ATV::cpu::gpw::PlusU;
    using ATV::cpu::gpw::XcKind;
    const auto basis = GaussianBasis(molrs::Element::H, "DZV-GTH-PADE");
    const auto pot = GthPotential(molrs::Element::H, "GTH-PADE-q1");
    const double L = 8.0;
    const double r0[3] = {0.5 * L - 0.7, 0.5 * L, 0.5 * L};
    const double r1[3] = {0.5 * L + 0.7, 0.5 * L, 0.5 * L};
    std::vector<GpwAtom> atoms{GpwAtom(r0, basis, pot, 1),
                               GpwAtom(r1, basis, pot, 1)};

    auto base = [&]() {
        GpwDft cfg;
        cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = L;
        cfg.cutoff_ry = 100.0;
        cfg.n_occ = 1;
        cfg.max_iter = 120;
        cfg.tol = 1e-7;
        cfg.mix = 0.35;
        cfg.xc = XcKind::PADE;
        cfg.k_points.push_back({0.25, 0.0, 0.0});
        cfg.k_point_weights.push_back(1.0);
        cfg.plus_u.push_back(PlusU{molrs::Element::H, /*l=*/0, /*u_ev=*/4.0});
        return cfg;
    };

    GpwDft diis = base();
    diis.use_ot = false;
    diis.build(atoms);
    const GpwResult r_diis = diis.run(atoms);
    ASSERT_TRUE(r_diis.converged) << "multi-k DIIS+U iters=" << r_diis.iters
                                  << " E=" << r_diis.e_total;
    EXPECT_GT(std::fabs(r_diis.e_dftu), 1e-8);

    GpwDft ot = base();
    ot.use_ot = true;
    ot.max_iter = 200;
    ot.build(atoms);
    const GpwResult r_ot = ot.run(atoms);
    ASSERT_TRUE(r_ot.converged) << "multi-k OT+U iters=" << r_ot.iters
                                << " E=" << r_ot.e_total;
    EXPECT_GT(std::fabs(r_ot.e_dftu), 1e-8);
    EXPECT_NEAR(r_ot.e_total, r_diis.e_total, 1e-4)
        << "OT+U=" << r_ot.e_total << " DIIS+U=" << r_diis.e_total
        << " iters ot/diis=" << r_ot.iters << "/" << r_diis.iters;
    EXPECT_NEAR(r_ot.e_dftu, r_diis.e_dftu, 1e-3);
}
