/// @file test_gradient.cpp
/// @brief Tests for Boys F_n, integral derivatives, and RHF/UHF analytic gradients.

#include <gtest/gtest.h>
#include <cmath>
#include <xtensor/xio.hpp>
#include <xtensor/xmath.hpp>

#include "atomiverse/core/gto.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/constants.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/integral_derivatives.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/gradient.h"

using namespace ATV;

// ── Helper: build H2 frame ───────────────────────────────────────────

static cpu::Frame make_h2(Real bond_length) {
    cpu::Frame f(2);
    f.x = {0.0, bond_length};
    f.y = {0.0, 0.0};
    f.z = {0.0, 0.0};
    f.mass = {1.008, 1.008};
    f.type_id = {1, 1};
    f.box = Mat3::cubic(100.0);
    return f;
}

static cpu::Frame make_heh_plus() {
    // HeH+ at 0.77 Angstrom
    cpu::Frame f(2);
    f.x = {0.0, 0.77};
    f.y = {0.0, 0.0};
    f.z = {0.0, 0.0};
    f.mass = {4.003, 1.008};
    f.type_id = {2, 1};  // He, H
    f.box = Mat3::cubic(100.0);
    return f;
}

// ── Boys function tests ──────────────────────────────────────────────

TEST(BoysFunction, F0KnownValues) {
    EXPECT_NEAR(boys_f0(0.0), 1.0, 1e-12);
    EXPECT_NEAR(boys_f0(1.0), 0.7468241328, 1e-8);
    // F_0(10) = 0.5*sqrt(pi/10)*erf(sqrt(10)) ≈ 0.28025
    EXPECT_NEAR(boys_f0(10.0), 0.28024739051, 1e-8);
}

TEST(BoysFunction, FnLimits) {
    // F_n(0) = 1/(2n+1)
    EXPECT_NEAR(boys_fn(0, 0.0), 1.0, 1e-12);
    EXPECT_NEAR(boys_fn(1, 0.0), 1.0/3.0, 1e-12);
    EXPECT_NEAR(boys_fn(2, 0.0), 1.0/5.0, 1e-12);
    EXPECT_NEAR(boys_fn(3, 0.0), 1.0/7.0, 1e-12);
}

TEST(BoysFunction, F1Numerical) {
    // F_1(x) = [F_0(x) - exp(-x)] / (2x)
    // Verify at several points via direct formula
    for (double x : {0.1, 0.5, 1.0, 5.0, 10.0}) {
        double F0 = boys_f0(x);
        double F1_direct = (F0 - std::exp(-x)) / (2.0 * x);
        EXPECT_NEAR(boys_fn(1, x), F1_direct, 1e-10)
            << "F_1 mismatch at x=" << x;
    }
}

TEST(BoysFunction, RecurrenceConsistency) {
    // Verify: F_n(x) = [(2n-1)*F_{n-1}(x) - exp(-x)] / (2x)
    for (double x : {0.01, 0.5, 2.0, 10.0}) {
        double ex = std::exp(-x);
        for (int n = 1; n <= 4; ++n) {
            double fn = boys_fn(n, x);
            double fn_from_recurrence = ((2*n - 1) * boys_fn(n-1, x) - ex) / (2.0 * x);
            EXPECT_NEAR(fn, fn_from_recurrence, 1e-10)
                << "Recurrence failed at n=" << n << ", x=" << x;
        }
    }
}

// ── Integral derivative tests ────────────────────────────────────────

TEST(IntegralDerivatives, TranslationalInvariance) {
    // For H2: dS/dR_1 + dS/dR_2 = 0 (translational invariance)
    auto frame = make_h2(0.74);
    auto dpack = cpu::scf::build_integral_derivatives_sto3g(frame);

    for (int d = 0; d < 3; ++d)
        for (Int mu = 0; mu < dpack.n_basis; ++mu)
            for (Int nu = 0; nu < dpack.n_basis; ++nu) {
                Real sum_dS = 0, sum_dT = 0, sum_dV = 0;
                for (int A = 0; A < dpack.n_atoms; ++A) {
                    sum_dS += dpack.dS(A, d, mu, nu);
                    sum_dT += dpack.dT(A, d, mu, nu);
                    sum_dV += dpack.dV(A, d, mu, nu);
                }
                EXPECT_NEAR(static_cast<double>(sum_dS), 0.0, 1e-10)
                    << "dS translational invariance violated at d=" << d
                    << " mu=" << mu << " nu=" << nu;
                EXPECT_NEAR(static_cast<double>(sum_dT), 0.0, 1e-10)
                    << "dT translational invariance violated";
                EXPECT_NEAR(static_cast<double>(sum_dV), 0.0, 1e-8)
                    << "dV translational invariance violated";
            }
}

TEST(IntegralDerivatives, ERITranslationalInvariance) {
    auto frame = make_h2(0.74);
    auto dpack = cpu::scf::build_integral_derivatives_sto3g(frame);
    Int n = dpack.n_basis;

    // Sum over all atoms of dERI should be zero
    for (int d = 0; d < 3; ++d) {
        Real max_err = 0;
        for (Int i = 0; i < n; ++i)
          for (Int j = 0; j < n; ++j)
            for (Int k = 0; k < n; ++k)
              for (Int l = 0; l < n; ++l) {
                Real s = 0;
                for (int A = 0; A < dpack.n_atoms; ++A)
                    s += dpack.dERI(A, d, i, j, k, l);
                max_err = std::max(max_err, std::abs(s));
              }
        EXPECT_LT(static_cast<double>(max_err), 1e-6)
            << "ERI translational invariance violated at d=" << d;
    }
}

// ── Integral derivatives vs numerical differences ────────────────────

TEST(IntegralDerivatives, OverlapDerivVsNumerical) {
    // Compare analytic dS/dR against numerical (S(R+h)-S(R-h))/(2h).
    // Note: Real may be float, so FD precision is limited.
    // Use h=1e-3 A for float-safe FD (~1e-3 Bohr perturbation).
    double h = 1e-3;  // Angstrom — use double for FD computation
    auto frame = make_h2(0.74);
    auto dpack = cpu::scf::build_integral_derivatives_sto3g(frame);

    auto frame_fwd = frame; frame_fwd.x(0) += Real(h);
    auto frame_bwd = frame; frame_bwd.x(0) -= Real(h);
    auto pack_fwd = cpu::scf::build_integrals_sto3g(frame_fwd);
    auto pack_bwd = cpu::scf::build_integrals_sto3g(frame_bwd);

    Int n = dpack.n_basis;
    for (Int mu = 0; mu < n; ++mu)
        for (Int nu = 0; nu < n; ++nu) {
            double dS_num = (static_cast<double>(pack_fwd.S(mu, nu))
                           - static_cast<double>(pack_bwd.S(mu, nu)))
                          / (2.0 * h * kAngstromToBohr);
            double dS_ana = static_cast<double>(dpack.dS(0, 0, mu, nu));
            // For float: FD truncation O(h^2) ~ 1e-6, float roundoff ~ 1e-4
            EXPECT_NEAR(dS_ana, dS_num, 1e-3)
                << "dS mismatch at mu=" << mu << " nu=" << nu;
        }
}

TEST(IntegralDerivatives, HcoreDerivVsNumerical) {
    double h = 1e-3;
    auto frame = make_h2(0.74);
    auto dpack = cpu::scf::build_integral_derivatives_sto3g(frame);

    auto frame_fwd = frame; frame_fwd.x(0) += Real(h);
    auto frame_bwd = frame; frame_bwd.x(0) -= Real(h);
    auto pack_fwd = cpu::scf::build_integrals_sto3g(frame_fwd);
    auto pack_bwd = cpu::scf::build_integrals_sto3g(frame_bwd);

    Int n = dpack.n_basis;
    for (Int mu = 0; mu < n; ++mu)
        for (Int nu = 0; nu < n; ++nu) {
            double dH_num = (static_cast<double>(pack_fwd.Hcore(mu, nu))
                           - static_cast<double>(pack_bwd.Hcore(mu, nu)))
                          / (2.0 * h * kAngstromToBohr);
            double dH_ana = static_cast<double>(dpack.dT(0, 0, mu, nu))
                          + static_cast<double>(dpack.dV(0, 0, mu, nu));
            EXPECT_NEAR(dH_ana, dH_num, 1e-2)
                << "dHcore mismatch at mu=" << mu << " nu=" << nu;
        }
}

TEST(IntegralDerivatives, ERIDerivVsNumerical) {
    Real h = 1e-5;
    auto frame = make_h2(0.74);
    auto dpack = cpu::scf::build_integral_derivatives_sto3g(frame);

    auto frame_fwd = frame; frame_fwd.x(0) += h;
    auto frame_bwd = frame; frame_bwd.x(0) -= h;
    auto pack_fwd = cpu::scf::build_integrals_sto3g(frame_fwd);
    auto pack_bwd = cpu::scf::build_integrals_sto3g(frame_bwd);

    Int n = dpack.n_basis;
    double max_err = 0;
    for (Int i = 0; i < n; ++i)
      for (Int j = 0; j < n; ++j)
        for (Int k = 0; k < n; ++k)
          for (Int l = 0; l < n; ++l) {
            Real dERI_num = (pack_fwd.ERI(i,j,k,l) - pack_bwd.ERI(i,j,k,l)) / (2.0 * h * kAngstromToBohr);
            Real dERI_ana = dpack.dERI(0, 0, i, j, k, l);
            double err = std::abs(static_cast<double>(dERI_ana) - static_cast<double>(dERI_num));
            max_err = std::max(max_err, err);
            EXPECT_NEAR(static_cast<double>(dERI_ana), static_cast<double>(dERI_num), 1e-3)
                << "dERI mismatch at (" << i << "," << j << "|" << k << "," << l << ")";
          }
}

TEST(IntegralDerivatives, VnucDerivVsNumerical) {
    double h = 1e-3;
    auto frame = make_h2(0.74);
    auto dpack = cpu::scf::build_integral_derivatives_sto3g(frame);

    auto frame_fwd = frame; frame_fwd.x(0) += Real(h);
    auto frame_bwd = frame; frame_bwd.x(0) -= Real(h);
    auto pack_fwd = cpu::scf::build_integrals_sto3g(frame_fwd);
    auto pack_bwd = cpu::scf::build_integrals_sto3g(frame_bwd);

    double dVnuc_num = (static_cast<double>(pack_fwd.V_nuc)
                      - static_cast<double>(pack_bwd.V_nuc))
                     / (2.0 * h * kAngstromToBohr);
    double dVnuc_ana = static_cast<double>(dpack.dV_nuc(0, 0));
    EXPECT_NEAR(dVnuc_ana, dVnuc_num, 1e-3)
        << "dV_nuc mismatch";
}

// ── Analytic vs finite-difference gradient ───────────────────────────

TEST(RHFGradient, H2AnalyticVsFiniteDifference) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::rhf_solve(integrals, Real(1e-10), 200);
    ASSERT_TRUE(result.converged);

    auto analytic = cpu::scf::rhf_gradient(frame, integrals, result);
    auto fd = cpu::scf::rhf_finite_difference_gradient(frame, 0, Real(1e-10), 200, Real(5e-5));

    for (int A = 0; A < 2; ++A)
        for (int d = 0; d < 3; ++d)
            EXPECT_NEAR(
                static_cast<double>(analytic.forces(A, d)),
                static_cast<double>(fd.forces(A, d)),
                0.5)  // float precision limits FD accuracy
                << "H2 gradient mismatch at atom " << A << " dim " << d;
}

TEST(RHFGradient, HeHPlusAnalyticVsFiniteDifference) {
    auto frame = make_heh_plus();
    auto integrals = cpu::scf::build_integrals_sto3g(frame, +1);
    auto result = cpu::scf::rhf_solve(integrals, Real(1e-10), 200);
    ASSERT_TRUE(result.converged);

    auto analytic = cpu::scf::rhf_gradient(frame, integrals, result);
    auto fd = cpu::scf::rhf_finite_difference_gradient(frame, +1, Real(1e-10), 200, Real(5e-5));

    for (int A = 0; A < 2; ++A)
        for (int d = 0; d < 3; ++d)
            EXPECT_NEAR(
                static_cast<double>(analytic.forces(A, d)),
                static_cast<double>(fd.forces(A, d)),
                0.5)  // float precision limits FD accuracy
                << "HeH+ gradient mismatch at atom " << A << " dim " << d;
}

// ── Zero net force ───────────────────────────────────────────────────

TEST(RHFGradient, ZeroNetForce) {
    // Newton's third law: sum of all forces must be zero
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::rhf_solve(integrals, Real(1e-10), 200);
    auto grad = cpu::scf::rhf_gradient(frame, integrals, result);

    Real fx_sum = 0, fy_sum = 0, fz_sum = 0;
    for (int A = 0; A < 2; ++A) {
        fx_sum += grad.forces(A, 0);
        fy_sum += grad.forces(A, 1);
        fz_sum += grad.forces(A, 2);
    }
    EXPECT_NEAR(static_cast<double>(fx_sum), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(fy_sum), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(fz_sum), 0.0, 1e-6);
}

TEST(RHFGradient, HeHPlusZeroNetForce) {
    auto frame = make_heh_plus();
    auto integrals = cpu::scf::build_integrals_sto3g(frame, +1);
    auto result = cpu::scf::rhf_solve(integrals, Real(1e-10), 200);
    auto grad = cpu::scf::rhf_gradient(frame, integrals, result);

    Real fx_sum = 0, fy_sum = 0, fz_sum = 0;
    for (int A = 0; A < 2; ++A) {
        fx_sum += grad.forces(A, 0);
        fy_sum += grad.forces(A, 1);
        fz_sum += grad.forces(A, 2);
    }
    EXPECT_NEAR(static_cast<double>(fx_sum), 0.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(fy_sum), 0.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(fz_sum), 0.0, 1e-4);
}

// ── Force direction ──────────────────────────────────────────────────

TEST(RHFGradient, H2ForceDirection) {
    // At equilibrium (~0.74 A), forces should be small.
    // At a stretched bond (1.5 A), atoms should attract (force points inward).
    auto frame = make_h2(1.5);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::rhf_solve(integrals, Real(1e-10), 200);
    auto grad = cpu::scf::rhf_gradient(frame, integrals, result);

    // Atom 0 at x=0, atom 1 at x=1.5. Attractive force on atom 0 points +x.
    EXPECT_GT(static_cast<double>(grad.forces(0, 0)), 0.0)
        << "H2 stretched: atom 0 should feel attractive force in +x";
    EXPECT_LT(static_cast<double>(grad.forces(1, 0)), 0.0)
        << "H2 stretched: atom 1 should feel attractive force in -x";
}

// ── Energy-force consistency ─────────────────────────────────────────

TEST(RHFGradient, EnergyForceConsistency) {
    // E(R + delta) - E(R) ~ -F . delta  for small delta
    Real R0 = 0.74;
    Real delta = 0.001;  // Angstrom

    auto frame0 = make_h2(R0);
    auto int0 = cpu::scf::build_integrals_sto3g(frame0);
    auto res0 = cpu::scf::rhf_solve(int0, Real(1e-10), 200);
    auto grad0 = cpu::scf::rhf_gradient(frame0, int0, res0);

    auto frame1 = make_h2(R0 + delta);
    auto int1 = cpu::scf::build_integrals_sto3g(frame1);
    auto res1 = cpu::scf::rhf_solve(int1, Real(1e-10), 200);

    Real dE_actual = res1.energy * kHartree2eV - res0.energy * kHartree2eV;
    // Force on atom 1 in x direction, displacement is +delta in x
    Real dE_predicted = -grad0.forces(1, 0) * delta;

    EXPECT_NEAR(static_cast<double>(dE_actual), static_cast<double>(dE_predicted), 0.001)
        << "Energy-force consistency violated";
}

// ── UHF gradient ─────────────────────────────────────────────────────

TEST(UHFGradient, H2ZeroNetForce) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(integrals, 1, 1, Real(1e-10), 200);
    ASSERT_TRUE(result.converged);
    auto grad = cpu::scf::uhf_gradient(frame, integrals, result);

    Real fx_sum = grad.forces(0, 0) + grad.forces(1, 0);
    Real fy_sum = grad.forces(0, 1) + grad.forces(1, 1);
    Real fz_sum = grad.forces(0, 2) + grad.forces(1, 2);
    EXPECT_NEAR(static_cast<double>(fx_sum), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(fy_sum), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(fz_sum), 0.0, 1e-6);
}

TEST(UHFGradient, H2AnalyticVsFiniteDifference) {
    auto frame = make_h2(0.74);
    auto analytic_integrals = cpu::scf::build_integrals_sto3g(frame);
    auto analytic_result = cpu::scf::uhf_solve(analytic_integrals, 1, 1, Real(1e-10), 200);
    ASSERT_TRUE(analytic_result.converged);
    auto analytic = cpu::scf::uhf_gradient(frame, analytic_integrals, analytic_result);
    auto fd = cpu::scf::uhf_finite_difference_gradient(frame, 1, 1, 0, Real(1e-10), 200, Real(5e-5));

    for (int A = 0; A < 2; ++A)
        for (int d = 0; d < 3; ++d)
            EXPECT_NEAR(
                static_cast<double>(analytic.forces(A, d)),
                static_cast<double>(fd.forces(A, d)),
                0.5)  // float precision limits FD accuracy
                << "UHF H2 gradient mismatch at atom " << A << " dim " << d;
}

TEST(UHFGradient, H2ForceDirection) {
    // Stretched H2: atoms should attract
    auto frame = make_h2(1.5);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(integrals, 1, 1, Real(1e-10), 200);
    auto grad = cpu::scf::uhf_gradient(frame, integrals, result);

    EXPECT_GT(static_cast<double>(grad.forces(0, 0)), 0.0)
        << "UHF H2 stretched: atom 0 should feel attractive force in +x";
    EXPECT_LT(static_cast<double>(grad.forces(1, 0)), 0.0)
        << "UHF H2 stretched: atom 1 should feel attractive force in -x";
}
