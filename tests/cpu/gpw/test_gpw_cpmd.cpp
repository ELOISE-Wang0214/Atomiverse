/// @file tests/cpu/gpw/test_gpw_cpmd.cpp
/// @brief Unit tests: Cpmd electronic propagator (RATTLE) + GPWEvaluator guards.
///
/// No multi-step MD / full water SCF trajectory here (that is the e2e
/// regressions/dynamics/cpmd-water and cpmd-rattle-water pins). Matrix-level
/// and build() validation only.

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/generators/xbuilder.hpp>

#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/cpmd.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gpw_solver.h"
#include "atomiverse/cpu/gpw/gth_potential.h"
#include "atomiverse/cpu/gpw/internal/collocate.h"
#include "atomiverse/cpu/linalg/matrix.h"

using ATV::cpu::gpw::Cpmd;
using ATV::cpu::linalg::RealMatrix;

namespace {

double max_abs_mat(const RealMatrix& a) {
    double m = 0.0;
    for (double v : a) {
        const double av = std::fabs(v);
        if (av > m)
            m = av;
    }
    return m;
}

RealMatrix eye_n(int n) {
    RealMatrix I = xt::zeros<double>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    for (int i = 0; i < n; ++i)
        I(i, i) = 1.0;
    return I;
}

} // namespace

// ── Cpmd propagator ─────────────────────────────────────────────────────

TEST(Cpmd, ConstructsWithPositiveMass) {
    const Cpmd cp(400.0);
    EXPECT_DOUBLE_EQ(cp.mu(), 400.0);
    EXPECT_FALSE(cp.ready());
}

TEST(Cpmd, NonPositiveMassThrows) {
    EXPECT_THROW(Cpmd(0.0), std::invalid_argument);
    EXPECT_THROW(Cpmd(-10.0), std::invalid_argument);
}

TEST(Cpmd, SeedThenReadyAtRest) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({4, 2});
    C(0, 0) = 1.0;
    C(1, 1) = 1.0;
    cp.seed(C);
    ASSERT_TRUE(cp.ready());
    EXPECT_EQ(cp.C().shape(0), 4u);
    EXPECT_EQ(cp.C().shape(1), 2u);
    EXPECT_NEAR(cp.kinetic_ha(), 0.0, 1e-15);
}

TEST(Cpmd, AdvanceBeforeSeedThrows) {
    Cpmd cp(400.0);
    RealMatrix H = xt::zeros<double>({2, 2});
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = 1.0;
    S(1, 1) = 1.0;
    EXPECT_THROW(cp.advance(H, S, 0.1), std::logic_error);
}

TEST(Cpmd, AdvanceWithZeroHKeepsRest) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    RealMatrix H = xt::zeros<double>({2, 2});
    RealMatrix S = xt::zeros<double>({2, 2});
    S(0, 0) = 1.0;
    S(1, 1) = 1.0;
    cp.advance(H, S, 0.1);
    EXPECT_NEAR(cp.C()(0, 0), 1.0, 1e-12);
    EXPECT_NEAR(cp.C()(1, 0), 0.0, 1e-12);
    EXPECT_NEAR(cp.kinetic_ha(), 0.0, 1e-12);
}

TEST(Cpmd, EmptySeedThrows) {
    Cpmd cp(400.0);
    RealMatrix empty;
    EXPECT_THROW(cp.seed(empty), std::invalid_argument);
}

TEST(Cpmd, ResetDropsReady) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    cp.reset();
    EXPECT_FALSE(cp.ready());
    EXPECT_THROW(cp.C(), std::logic_error);
}

// ── Electronic RATTLE residual gates ────────────────────────────────────

TEST(Cpmd, RattlePreservesOrthonormalityAndVelocityConstraint) {
    // Two occupied columns in a 3-AO space; S = I; nontrivial diagonal H so
    // the force is nonzero off the BO surface for a non-eigen C.
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({3, 2});
    C(0, 0) = 1.0;
    C(1, 1) = 1.0;
    cp.seed(C);

    RealMatrix S = eye_n(3);
    RealMatrix H = xt::zeros<double>({3, 3});
    H(0, 0) = -1.0;
    H(1, 1) = -0.5;
    H(2, 2) = 0.2;
    H(0, 2) = H(2, 0) = 0.15; // couple virtual direction into force

    cp.advance(H, S, 0.05);
    ASSERT_TRUE(cp.ready());
    EXPECT_GE(cp.last_position_iters(), 0);
    EXPECT_LE(cp.last_constraint_residual(), 1e-9);

    using xt::linalg::dot;
    const RealMatrix G = dot(xt::transpose(cp.C()), dot(S, cp.C()));
    EXPECT_LE(max_abs_mat(G - eye_n(2)), 1e-9);

    // Velocity constraint: sym(Cᵀ S Ċ) ≈ 0.
    const RealMatrix W = dot(xt::transpose(cp.C()), dot(S, cp.Cdot()));
    const RealMatrix Wsym = 0.5 * (W + xt::transpose(W));
    EXPECT_LE(max_abs_mat(Wsym), 1e-9);

    // Residual accessor matches the recomputed ‖σ‖_∞.
    EXPECT_NEAR(cp.last_constraint_residual(), max_abs_mat(G - eye_n(2)),
                1e-15);
}

TEST(Cpmd, RattleAtRestZeroHKeepsZeroVelocity) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    RealMatrix H = xt::zeros<double>({2, 2});
    RealMatrix S = eye_n(2);
    cp.advance(H, S, 0.1);
    EXPECT_NEAR(cp.kinetic_ha(), 0.0, 1e-12);
    EXPECT_NEAR(cp.Cdot()(0, 0), 0.0, 1e-12);
    EXPECT_NEAR(cp.Cdot()(1, 0), 0.0, 1e-12);
    EXPECT_LE(cp.last_constraint_residual(), 1e-12);
}

// ── GPWEvaluator cpmd guards ────────────────────────────────────────────

namespace {

ATV::cpu::Frame water_frame() {
    ATV::cpu::Frame f(3);
    const int z[3] = {8, 1, 1};
    const double x[3] = {5.0000, 5.7586, 4.2414};
    const double y[3] = {5.0000, 5.0000, 5.0000};
    const double zc[3] = {5.0000, 5.5043, 5.5043};
    for (int i = 0; i < 3; ++i) {
        f.span(ATV::kAtom_type)[i] = z[i];
        f.span(ATV::kAtom_x)[i] = x[i];
        f.span(ATV::kAtom_y)[i] = y[i];
        f.span(ATV::kAtom_z)[i] = zc[i];
    }
    f.simbox = ATV::cpu::Box::cubic(10.0);
    return f;
}

ATV::cpu::gpw::GPWEvaluator water_eval() {
    const ATV::cpu::gpw::Kind O{
        molrs::Element::O,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
    const ATV::cpu::gpw::Kind H{
        molrs::Element::H,
        ATV::cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
        ATV::cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};
    return ATV::cpu::gpw::GPWEvaluator(
        {O, H}, /*cutoff_ry*/ 200.0, /*max_iter*/ 60, /*tol*/ 1e-5,
        /*mix*/ 0.3, ATV::cpu::gpw::XcKind::PBE, /*nonperiodic*/ true,
        /*compute_forces*/ true);
}

} // namespace

TEST(GpwEvaluatorCpmd, BuildAcceptsCpmdWithUks) {
    // Γ-UKS CPMD is product (dual α/β propagators); O2-like occupations.
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorCpmd, BuildRejectsUksCpmdWithZeroAlpha) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 0;
    eval.n_beta = 3;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(Cpmd, AdvanceWithOccOneMatchesUksFactor) {
    // occ=1 uses ot_gradient factor 2 H C (UKS spin channel).
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    RealMatrix H = xt::zeros<double>({2, 2});
    H(0, 0) = -1.0;
    H(1, 1) = 0.5;
    H(0, 1) = H(1, 0) = 0.1;
    RealMatrix S = eye_n(2);
    // At BO-ish eigen direction force is small; just ensure finite advance.
    EXPECT_NO_THROW(cp.advance(H, S, 0.05, /*occ=*/1.0));
    EXPECT_TRUE(std::isfinite(cp.kinetic_ha()));
    EXPECT_GE(cp.kinetic_ha(), 0.0);
}

TEST(Cpmd, AdvanceWithNonPositiveOccThrows) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    RealMatrix H = eye_n(2);
    RealMatrix S = eye_n(2);
    EXPECT_THROW(cp.advance(H, S, 0.1, /*occ=*/0.0), std::invalid_argument);
    EXPECT_THROW(cp.advance(H, S, 0.1, /*occ=*/-1.0), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithOt) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.use_ot = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithXlBomd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.xl_bomd = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsNonPositiveMuOrDt) {
    auto frame = water_frame();
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_mu = 0.0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_dt = -0.1;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
}

TEST(GpwEvaluatorCpmd, BuildRejectsCpmdWithoutForces) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.compute_forces = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// ── Nosé–Hoover electronic thermostat ───────────────────────────────────

TEST(Cpmd, ThermostatInvalidParamsThrowOnAdvance) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    cp.thermostat = true;
    cp.ekin_target_ha = 0.0;
    cp.tau_e_fs = 10.0;
    RealMatrix H = eye_n(2);
    RealMatrix S = eye_n(2);
    EXPECT_THROW(cp.advance(H, S, 0.1), std::invalid_argument);
    cp.ekin_target_ha = 1e-4;
    cp.tau_e_fs = -1.0;
    EXPECT_THROW(cp.advance(H, S, 0.1), std::invalid_argument);
}

TEST(Cpmd, SeedWithCdotSetsKineticEnergy) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    RealMatrix Cdot = xt::zeros<double>({2, 1});
    Cdot(1, 0) = 0.01; // off-manifold component is fine for KE check
    cp.seed(C, Cdot);
    const double expect = 0.5 * 400.0 * (0.01 * 0.01);
    EXPECT_NEAR(cp.kinetic_ha(), expect, 1e-15);
    EXPECT_NEAR(cp.thermostat_energy_ha(), 0.0, 1e-15);
    EXPECT_NEAR(cp.xi(), 0.0, 1e-15);
    EXPECT_NEAR(cp.eta(), 0.0, 1e-15);
}

TEST(Cpmd, ThermostatRegulatesHotElectronicKe) {
    // Hot seed Te ≫ K0; NH should pull Te down over several steps.
    Cpmd cp(400.0);
    cp.thermostat = true;
    cp.ekin_target_ha = 1e-6;
    cp.tau_e_fs = 5.0; // short τ for aggressive regulation on unit scale
    RealMatrix C = xt::zeros<double>({3, 2});
    C(0, 0) = 1.0;
    C(1, 1) = 1.0;
    RealMatrix Cdot = xt::zeros<double>({3, 2});
    Cdot(2, 0) = 0.05;
    Cdot(2, 1) = -0.04;
    cp.seed(C, Cdot);
    const double te0 = cp.kinetic_ha();
    ASSERT_GT(te0, 10.0 * cp.ekin_target_ha);

    RealMatrix S = eye_n(3);
    RealMatrix H = xt::zeros<double>({3, 3});
    H(0, 0) = -1.0;
    H(1, 1) = -0.5;
    H(2, 2) = 0.1;

    for (int i = 0; i < 40; ++i)
        cp.advance(H, S, 0.05);

    const double tef = cp.kinetic_ha();
    EXPECT_TRUE(std::isfinite(tef));
    EXPECT_GE(tef, 0.0);
    // Pull toward target: final below initial and not runaway.
    EXPECT_LT(tef, te0);
    EXPECT_LE(tef, 50.0 * cp.ekin_target_ha + 0.1 * te0);
    EXPECT_TRUE(std::isfinite(cp.thermostat_energy_ha()));
    EXPECT_LE(cp.last_constraint_residual(), 1e-9);

    using xt::linalg::dot;
    const RealMatrix G = dot(xt::transpose(cp.C()), dot(S, cp.C()));
    EXPECT_LE(max_abs_mat(G - eye_n(2)), 1e-9);
    const RealMatrix W = dot(xt::transpose(cp.C()), dot(S, cp.Cdot()));
    const RealMatrix Wsym = 0.5 * (W + xt::transpose(W));
    EXPECT_LE(max_abs_mat(Wsym), 1e-9);
}

TEST(Cpmd, ThermostatOffKeepsZeroThermostatEnergy) {
    Cpmd cp(400.0);
    RealMatrix C = xt::zeros<double>({2, 1});
    C(0, 0) = 1.0;
    cp.seed(C);
    EXPECT_FALSE(cp.thermostat);
    EXPECT_NEAR(cp.thermostat_energy_ha(), 0.0, 1e-15);
    RealMatrix H = eye_n(2);
    H(0, 0) = -0.5;
    RealMatrix S = eye_n(2);
    cp.advance(H, S, 0.1);
    EXPECT_NEAR(cp.thermostat_energy_ha(), 0.0, 1e-15);
    EXPECT_NEAR(cp.xi(), 0.0, 1e-15);
    EXPECT_NEAR(cp.eta(), 0.0, 1e-15);
}

TEST(GpwEvaluatorCpmd, BuildRejectsThermostatWithoutCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = false;
    eval.cpmd_thermostat = true;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsBadThermostatParams) {
    auto frame = water_frame();
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_thermostat = true;
        eval.cpmd_ekin_target_ha = 0.0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
    {
        auto eval = water_eval();
        eval.cpmd = true;
        eval.cpmd_thermostat = true;
        eval.cpmd_tau_e = -5.0;
        ATV::cpu::Context rt;
        ATV::cpu::BuildContext bc{rt, frame};
        EXPECT_THROW(eval.build(bc), std::invalid_argument);
    }
}

TEST(GpwEvaluatorCpmd, BuildAcceptsCpmdWithThermostat) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.cpmd_thermostat = true;
    eval.cpmd_ekin_target_ha = 1e-5;
    eval.cpmd_tau_e = 50.0;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

// ── Multi-k Bloch CPMD (CpmdCplx + evaluator guards) ────────────────────

using ATV::cpu::gpw::CpmdCplx;
using ATV::cpu::linalg::CplxMatrix;
using Cplx = std::complex<double>;

namespace {

double max_abs_cmat(const CplxMatrix& a) {
    double m = 0.0;
    for (const Cplx& v : a) {
        const double av = std::abs(v);
        if (av > m)
            m = av;
    }
    return m;
}

CplxMatrix eye_c(int n) {
    CplxMatrix I = xt::zeros<Cplx>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(n)});
    for (int i = 0; i < n; ++i)
        I(i, i) = Cplx(1.0, 0.0);
    return I;
}

// A† B (small helper for residual checks).
CplxMatrix c_dagger_mul_test(const CplxMatrix& A, const CplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    CplxMatrix C = xt::zeros<Cplx>(
        {static_cast<std::size_t>(m), static_cast<std::size_t>(p)});
    for (int i = 0; i < m; ++i)
        for (int k = 0; k < n; ++k) {
            const Cplx a = std::conj(A(k, i));
            for (int j = 0; j < p; ++j)
                C(i, j) += a * B(k, j);
        }
    return C;
}

CplxMatrix c_mul_test(const CplxMatrix& A, const CplxMatrix& B) {
    const int n = static_cast<int>(A.shape(0));
    const int m = static_cast<int>(A.shape(1));
    const int p = static_cast<int>(B.shape(1));
    CplxMatrix C = xt::zeros<Cplx>(
        {static_cast<std::size_t>(n), static_cast<std::size_t>(p)});
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < m; ++k) {
            const Cplx a = A(i, k);
            for (int j = 0; j < p; ++j)
                C(i, j) += a * B(k, j);
        }
    return C;
}

} // namespace

TEST(CpmdCplx, ConstructsWithPositiveMass) {
    const CpmdCplx cp(400.0);
    EXPECT_DOUBLE_EQ(cp.mu(), 400.0);
    EXPECT_FALSE(cp.ready());
}

TEST(CpmdCplx, NonPositiveMassThrows) {
    EXPECT_THROW(CpmdCplx(0.0), std::invalid_argument);
    EXPECT_THROW(CpmdCplx(-1.0), std::invalid_argument);
}

TEST(CpmdCplx, EmptySeedThrows) {
    CpmdCplx cp(400.0);
    CplxMatrix empty;
    EXPECT_THROW(cp.seed(empty), std::invalid_argument);
}

TEST(CpmdCplx, RattlePreservesComplexOrthonormality) {
    CpmdCplx cp(400.0);
    CplxMatrix C = xt::zeros<Cplx>({3, 2});
    C(0, 0) = Cplx(1.0, 0.0);
    C(1, 1) = Cplx(0.0, 1.0); // complex phase on second column
    cp.seed(C);

    CplxMatrix S = eye_c(3);
    CplxMatrix H = xt::zeros<Cplx>({3, 3});
    H(0, 0) = Cplx(-1.0, 0.0);
    H(1, 1) = Cplx(-0.5, 0.0);
    H(2, 2) = Cplx(0.2, 0.0);
    H(0, 2) = Cplx(0.1, 0.05);
    H(2, 0) = std::conj(H(0, 2));

    cp.advance(H, S, 0.05);
    ASSERT_TRUE(cp.ready());
    EXPECT_LE(cp.last_constraint_residual(), 1e-9);

    const CplxMatrix G = c_dagger_mul_test(cp.C(), c_mul_test(S, cp.C()));
    CplxMatrix residual = G - eye_c(2);
    EXPECT_LE(max_abs_cmat(residual), 1e-9);

    // Velocity constraint: Herm(C† S Ċ) ≈ 0.
    const CplxMatrix W = c_dagger_mul_test(cp.C(), c_mul_test(S, cp.Cdot()));
    CplxMatrix Wherm = xt::zeros<Cplx>(W.shape());
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            Wherm(i, j) = 0.5 * (W(i, j) + std::conj(W(j, i)));
    EXPECT_LE(max_abs_cmat(Wherm), 1e-9);
    EXPECT_TRUE(std::isfinite(cp.kinetic_ha()));
    EXPECT_GE(cp.kinetic_ha(), 0.0);
}

TEST(CpmdCplx, AdvanceBeforeSeedThrows) {
    CpmdCplx cp(400.0);
    CplxMatrix H = eye_c(2);
    CplxMatrix S = eye_c(2);
    EXPECT_THROW(cp.advance(H, S, 0.1), std::logic_error);
}

TEST(GpwEvaluatorCpmd, BuildAcceptsMultiKRksCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.k_mesh[0] = 1;
    eval.k_mesh[1] = 1;
    eval.k_mesh[2] = 2;
    eval.k_time_reversal = false;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_NO_THROW(eval.build(bc));
}

TEST(GpwEvaluatorCpmd, BuildRejectsMultiKUksCpmd) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.n_alpha = 5;
    eval.n_beta = 3;
    eval.k_mesh[2] = 2;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

TEST(GpwEvaluatorCpmd, BuildRejectsMultiKCpmdThermostat) {
    auto frame = water_frame();
    auto eval = water_eval();
    eval.cpmd = true;
    eval.cpmd_thermostat = true;
    eval.k_mesh[2] = 2;
    eval.nonperiodic = false;
    ATV::cpu::Context rt;
    ATV::cpu::BuildContext bc{rt, frame};
    EXPECT_THROW(eval.build(bc), std::invalid_argument);
}

// Multi-k fixed-density payload (tiny H2 GpwDft).
namespace {

ATV::cpu::gpw::GridBasisFunction s_basis_cp(double alpha, const double c[3]) {
    ATV::cpu::gpw::GridPrimitive p;
    p.alpha = alpha;
    p.center[0] = c[0];
    p.center[1] = c[1];
    p.center[2] = c[2];
    p.l[0] = p.l[1] = p.l[2] = 0;
    p.coeff = std::pow(2.0 * alpha / M_PI, 0.75);
    return ATV::cpu::gpw::GridBasisFunction{p};
}

ATV::cpu::gpw::GpwAtom hydrogen_cp(const double pos[3]) {
    ATV::cpu::gpw::GpwAtom a;
    a.z_val = 1.0;
    a.z = 1;
    a.pos[0] = pos[0];
    a.pos[1] = pos[1];
    a.pos[2] = pos[2];
    a.basis = {s_basis_cp(0.5, pos)};
    a.rloc = 0.2;
    a.cloc = {-4.17890044, 0.72446331};
    return a;
}

std::vector<ATV::cpu::gpw::GpwAtom> small_h2_cp(double box = 8.0) {
    const double rA[3] = {box / 2 - 0.7, box / 2, box / 2};
    const double rB[3] = {box / 2 + 0.7, box / 2, box / 2};
    return {hydrogen_cp(rA), hydrogen_cp(rB)};
}

ATV::cpu::gpw::GpwDft small_cfg_cp(double box = 8.0) {
    ATV::cpu::gpw::GpwDft cfg;
    cfg.cell[0] = cfg.cell[1] = cfg.cell[2] = box;
    cfg.cutoff_ry = 120.0;
    cfg.n_occ = 1;
    cfg.max_iter = 80;
    cfg.tol = 1e-7;
    cfg.mix = 0.5;
    return cfg;
}

} // namespace

TEST(GpwSolverFixedDensityK, MultiKFixedDensityExportsHkSk) {
    const auto atoms = small_h2_cp(/*box=*/8.0);
    ATV::cpu::gpw::GpwDft cold = small_cfg_cp(/*box=*/8.0);
    cold.k_mesh[0] = 1;
    cold.k_mesh[1] = 1;
    cold.k_mesh[2] = 2;
    cold.k_time_reversal = false;
    cold.nonperiodic = false;
    cold.build(atoms);
    const auto ref = cold.run(atoms);
    ASSERT_TRUE(ref.converged);
    ASSERT_GE(ref.n_kpts, 2);
    ASSERT_EQ(static_cast<int>(ref.P_k.size()), ref.n_kpts);

    ATV::cpu::gpw::GpwDft fix = small_cfg_cp(/*box=*/8.0);
    fix.k_mesh[0] = 1;
    fix.k_mesh[1] = 1;
    fix.k_mesh[2] = 2;
    fix.k_time_reversal = false;
    fix.nonperiodic = false;
    fix.max_iter = 0;
    fix.build(atoms);
    const auto& r = fix.run(atoms, nullptr, nullptr, nullptr, &ref.P_k);
    EXPECT_TRUE(r.converged);
    EXPECT_EQ(r.iters, 0);
    EXPECT_TRUE(std::isfinite(r.e_total));
    EXPECT_EQ(static_cast<int>(r.H_ks_k.size()), r.n_kpts);
    EXPECT_EQ(static_cast<int>(r.S_k.size()), r.n_kpts);
    EXPECT_EQ(static_cast<int>(r.P_k.size()), r.n_kpts);
    // Energy near the cold SCF (same density).
    EXPECT_NEAR(r.e_total, ref.e_total, 1e-6);
}

TEST(GpwSolverFixedDensityK, MultiKFixedDensityWithoutGuessThrows) {
    const auto atoms = small_h2_cp();
    ATV::cpu::gpw::GpwDft fix = small_cfg_cp();
    fix.k_mesh[2] = 2;
    fix.k_time_reversal = false;
    fix.nonperiodic = false;
    fix.max_iter = 0;
    fix.build(atoms);
    EXPECT_THROW(fix.run(atoms), std::invalid_argument);
}
