/// @file tests/cpu/scf/test_uhf_solve.cpp
/// @brief Tests for UHF solver: open-shell atoms, RHF/UHF comparison, spin diagnostics.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include <xtensor/xnorm.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/rhf_solve.h"
#include "atomiverse/cpu/scf/uhf_solve.h"

using namespace ATV;

class UHFSolveTest : public ::testing::Test {
protected:
    static constexpr Real kTol = (sizeof(Real) == sizeof(float)) ? Real(1e-5) : Real(1e-9);
    static constexpr double kEnergyTol = 1e-3;

    static cpu::Frame make_h_atom() {
        cpu::Frame f(1);
        f.atoms.x[0] = Real(0); f.atoms.y[0] = Real(0); f.atoms.z[0] = Real(0);
        f.atoms.type[0] = 1;
        f.box = Mat3::cubic(100.0);
        return f;
    }

    static cpu::Frame make_he_atom() {
        cpu::Frame f(1);
        f.atoms.x[0] = Real(0); f.atoms.y[0] = Real(0); f.atoms.z[0] = Real(0);
        f.atoms.type[0] = 2;
        f.box = Mat3::cubic(100.0);
        return f;
    }

    static cpu::Frame make_li_atom() {
        cpu::Frame f(1);
        f.atoms.x[0] = Real(0); f.atoms.y[0] = Real(0); f.atoms.z[0] = Real(0);
        f.atoms.type[0] = 3;
        f.box = Mat3::cubic(100.0);
        return f;
    }

    static cpu::Frame make_h2(double bond_length_angstrom) {
        cpu::Frame f(2);
        f.atoms.x[0] = Real(0);                    f.atoms.y[0] = Real(0); f.atoms.z[0] = Real(0);
        f.atoms.x[1] = Real(bond_length_angstrom); f.atoms.y[1] = Real(0); f.atoms.z[1] = Real(0);
        f.atoms.type[0] = 1; f.atoms.type[1] = 1;
        f.box = Mat3::cubic(100.0);
        return f;
    }
};

// ── H atom ──────────────────────────────────────────────────────────

TEST_F(UHFSolveTest, HAtomConverges) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 1);
    EXPECT_EQ(result.n_beta, 0);
}

TEST_F(UHFSolveTest, HAtomEnergy) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(static_cast<double>(result.energy), -0.4666, 5e-3);
}

TEST_F(UHFSolveTest, HAtomElectronCount) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);
    Real tr_alpha = xt::sum(result.P_alpha * ints.S)();
    Real tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-5);
    EXPECT_NEAR(static_cast<double>(tr_beta), 0.0, 1e-5);
}

TEST_F(UHFSolveTest, HAtomSpinContamination) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 0, kTol, 50);
    EXPECT_NEAR(static_cast<double>(result.spin_contamination), 0.0, 1e-6);
}

// ── He atom ─────────────────────────────────────────────────────────

TEST_F(UHFSolveTest, HeAtomEqualsRHF) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto rhf = cpu::scf::rhf_solve(ints, kTol, 50);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);
    EXPECT_TRUE(uhf.converged);
    EXPECT_NEAR(static_cast<double>(uhf.energy),
                static_cast<double>(rhf.energy), kEnergyTol);
}

TEST_F(UHFSolveTest, HeAtomSpinContamination) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);
    EXPECT_NEAR(static_cast<double>(uhf.spin_contamination), 0.0, 1e-4);
}

// ── Li atom ─────────────────────────────────────────────────────────

TEST_F(UHFSolveTest, LiAtomConverges) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 2);
    EXPECT_EQ(result.n_beta, 1);
}

TEST_F(UHFSolveTest, LiAtomElectronCount) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);
    Real tr_alpha = xt::sum(result.P_alpha * ints.S)();
    Real tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 2.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 1.0, 1e-4);
}

TEST_F(UHFSolveTest, LiAtomEnergyRange) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);
    EXPECT_TRUE(result.converged);
    EXPECT_LT(static_cast<double>(result.energy), -7.0);
    EXPECT_GT(static_cast<double>(result.energy), -8.0);
}

// ── H2 equilibrium ──────────────────────────────────────────────────

TEST_F(UHFSolveTest, H2EquilibriumAgreesWithRHF) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto rhf = cpu::scf::rhf_solve(ints, kTol, 200);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 200);
    EXPECT_TRUE(uhf.converged);
    EXPECT_NEAR(static_cast<double>(uhf.energy),
                static_cast<double>(rhf.energy), kEnergyTol);
}

TEST_F(UHFSolveTest, H2EquilibriumSpinContamination) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 200);
    EXPECT_NEAR(static_cast<double>(uhf.spin_contamination), 0.0, 1e-3);
}

// ── H2 stretched ────────────────────────────────────────────────────

TEST_F(UHFSolveTest, H2StretchedUHFLowerThanRHF) {
    auto frame = make_h2(2.0);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto rhf = cpu::scf::rhf_solve(ints, kTol, 50);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 100);
    EXPECT_TRUE(uhf.converged);
    EXPECT_LE(static_cast<double>(uhf.energy),
              static_cast<double>(rhf.energy) + 1e-6);
}

TEST_F(UHFSolveTest, H2StretchedSpinContaminationNonzero) {
    auto frame = make_h2(2.0);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 100);
    EXPECT_TRUE(uhf.converged);
    EXPECT_GT(static_cast<double>(uhf.spin_contamination), 0.0);
}

// ── PES scan ────────────────────────────────────────────────────────

TEST_F(UHFSolveTest, H2PESScanUHFVariational) {
    for (double r : {0.6, 0.74, 1.0, 1.5, 2.0, 2.5, 3.0}) {
        auto frame = make_h2(r);
        auto ints = cpu::scf::build_integrals_sto3g(frame);
        auto rhf = cpu::scf::rhf_solve(ints, kTol, 200);
        auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 400);
        EXPECT_TRUE(uhf.converged) << "UHF failed to converge at R=" << r;
        EXPECT_LE(static_cast<double>(uhf.energy),
                  static_cast<double>(rhf.energy) + 1e-5)
            << "UHF > RHF at R=" << r;
    }
}

// ── Result structure ────────────────────────────────────────────────

TEST_F(UHFSolveTest, ResultContainsAlphaBetaEigenvalues) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);
    EXPECT_EQ(result.eigenvalues_alpha.size(), (size_t)ints.n_basis);
    EXPECT_EQ(result.eigenvalues_beta.size(), (size_t)ints.n_basis);
}

TEST_F(UHFSolveTest, ResultHasConvergenceTrace) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);
    EXPECT_GT(result.trace.size(), 0u);
    EXPECT_EQ((int)result.trace.back().iter, result.iterations);
}

TEST_F(UHFSolveTest, DensityMatricesHaveCorrectShape) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals_sto3g(frame);
    auto result = cpu::scf::uhf_solve(ints, 2, 1, kTol, 100);
    Int n = ints.n_basis;
    EXPECT_EQ(result.P_alpha.shape(0), (size_t)n);
    EXPECT_EQ(result.P_alpha.shape(1), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(0), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(1), (size_t)n);
}

#else

#include <gtest/gtest.h>
TEST(UHFSolve, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — UHF tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
