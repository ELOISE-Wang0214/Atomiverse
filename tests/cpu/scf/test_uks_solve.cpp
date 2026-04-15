/// @file tests/cpu/scf/test_uks_solve.cpp
/// @brief Tests for UKS solver: LSDA/PBE on H, He, Li atoms and H2 molecule.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include <xtensor/reducers/xnorm.hpp>
#include <xtensor-blas/xlinalg.hpp>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/uhf_solve.h"
#include "atomiverse/cpu/scf/uks_solve.h"

using namespace ATV;

class UKSSolveTest : public ::testing::Test {
protected:
    static constexpr double kTol = (sizeof(double) == sizeof(float)) ? double(1e-4) : double(1e-7);
    static constexpr double kEnergyTol = 5e-3;

    static cpu::Frame make_h_atom() {
        cpu::Frame f(1);
        f.atoms.x[0] = double(0); f.atoms.y[0] = double(0); f.atoms.z[0] = double(0);
        f.atoms.type[0] = 1;
        f.box = Mat3::cubic(100.0);
        return f;
    }

    static cpu::Frame make_he_atom() {
        cpu::Frame f(1);
        f.atoms.x[0] = double(0); f.atoms.y[0] = double(0); f.atoms.z[0] = double(0);
        f.atoms.type[0] = 2;
        f.box = Mat3::cubic(100.0);
        return f;
    }

    static cpu::Frame make_li_atom() {
        cpu::Frame f(1);
        f.atoms.x[0] = double(0); f.atoms.y[0] = double(0); f.atoms.z[0] = double(0);
        f.atoms.type[0] = 3;
        f.box = Mat3::cubic(100.0);
        return f;
    }

    static cpu::Frame make_h2(double bond_length_angstrom) {
        cpu::Frame f(2);
        f.atoms.x[0] = double(0);                    f.atoms.y[0] = double(0); f.atoms.z[0] = double(0);
        f.atoms.x[1] = double(bond_length_angstrom); f.atoms.y[1] = double(0); f.atoms.z[1] = double(0);
        f.atoms.type[0] = 1; f.atoms.type[1] = 1;
        f.box = Mat3::cubic(100.0);
        return f;
    }
};

// ── Result structure ────────────────────────────────────────────────

TEST_F(UKSSolveTest, ResultContainsXCDecomposition) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    EXPECT_NE(static_cast<double>(result.one_electron_energy), 0.0);
    EXPECT_NE(static_cast<double>(result.coulomb_energy), 0.0);
    EXPECT_NE(static_cast<double>(result.xc_energy), 0.0);
}

TEST_F(UKSSolveTest, ResultHasIntegratedElectronCounts) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    EXPECT_GT(static_cast<double>(result.integrated_n_alpha), 0.0);
    EXPECT_GT(static_cast<double>(result.integrated_n_beta), 0.0);
}

TEST_F(UKSSolveTest, DensityMatricesHaveCorrectShape) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    int n = ints.n_basis;
    EXPECT_EQ(result.P_alpha.shape(0), (size_t)n);
    EXPECT_EQ(result.P_alpha.shape(1), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(0), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(1), (size_t)n);
}

TEST_F(UKSSolveTest, EnergyDecompositionSumsCorrectly) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    double E_sum = static_cast<double>(result.one_electron_energy)
                 + static_cast<double>(result.coulomb_energy)
                 + static_cast<double>(result.xc_energy)
                 + static_cast<double>(ints.V_nuc);
    EXPECT_NEAR(E_sum, static_cast<double>(result.energy), 1e-6);
}

// ── H atom (doublet) ────────────────────────────────────────────────

TEST_F(UKSSolveTest, HAtomLSDAConverges) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 1);
    EXPECT_EQ(result.n_beta, 0);
}

TEST_F(UKSSolveTest, HAtomLSDAEnergyReasonable) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 50, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -0.40);
    EXPECT_GT(E, -0.55);
}

TEST_F(UKSSolveTest, HAtomSpinDensity) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 0, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    double tr_alpha = xt::sum(result.P_alpha * ints.S)();
    double tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 0.0, 1e-10);
}

// ── He atom (singlet) ───────────────────────────────────────────────

TEST_F(UKSSolveTest, HeAtomLSDAConverges) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
}

TEST_F(UKSSolveTest, HeAtomLSDAEnergyReasonable) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 50, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -2.70);
    EXPECT_GT(E, -2.95);
}

TEST_F(UKSSolveTest, HeClosedShellAlphaBetaSymmetry) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    double diff = xt::norm_l2(result.P_alpha - result.P_beta)();
    EXPECT_NEAR(static_cast<double>(diff), 0.0, 1e-4);
}

TEST_F(UKSSolveTest, HeAtomElectronCount) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(result.converged);
    double tr_alpha = xt::sum(result.P_alpha * ints.S)();
    double tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 1.0, 1e-4);
}

// ── Li atom (doublet) ───────────────────────────────────────────────

TEST_F(UKSSolveTest, LiAtomLSDAConverges) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 2, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 100);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 2);
    EXPECT_EQ(result.n_beta, 1);
}

TEST_F(UKSSolveTest, LiAtomLSDAEnergyReasonable) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 2, 1, cpu::scf::XCType::LSDA, 50, 110, kTol, 100);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -7.0);
    EXPECT_GT(E, -7.8);
}

// ── H2 at equilibrium ───────────────────────────────────────────────

TEST_F(UKSSolveTest, H2EquilibriumLSDAConverges) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 100);
    EXPECT_TRUE(result.converged);
}

TEST_F(UKSSolveTest, H2EquilibriumLSDAEnergyReasonable) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 50, 110, kTol, 100);
    EXPECT_TRUE(result.converged);
    double E = static_cast<double>(result.energy);
    EXPECT_LT(E, -1.05);
    EXPECT_GT(E, -1.25);
}

// ── Cross-method comparisons ────────────────────────────────────────

TEST_F(UKSSolveTest, UKSvsUHFDifferentEnergy) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto uhf = cpu::scf::uhf_solve(ints, 1, 1, kTol, 50);
    auto uks = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_TRUE(uhf.converged);
    EXPECT_TRUE(uks.converged);
    EXPECT_GT(std::abs(static_cast<double>(uks.energy) - static_cast<double>(uhf.energy)), 1e-4);
}

TEST_F(UKSSolveTest, HeAtomLSDAvsPBE) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto lsda = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    auto pbe = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::PBE, 35, 110, kTol, 80);
    EXPECT_TRUE(lsda.converged);
    EXPECT_TRUE(pbe.converged);
    EXPECT_GT(std::abs(static_cast<double>(lsda.energy) - static_cast<double>(pbe.energy)), 1e-4);
}

// ── Convergence and eigenvalues ─────────────────────────────────────

TEST_F(UKSSolveTest, ResultHasConvergenceTrace) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_GT(result.trace.size(), 0u);
    EXPECT_EQ(static_cast<int>(result.trace.back().iter), result.iterations);
}

TEST_F(UKSSolveTest, EigenvaluesHaveCorrectShape) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto result = cpu::scf::uks_solve(ints, 1, 1, cpu::scf::XCType::LSDA, 35, 110, kTol, 80);
    EXPECT_EQ(result.eigenvalues_alpha.size(), (size_t)ints.n_basis);
    EXPECT_EQ(result.eigenvalues_beta.size(), (size_t)ints.n_basis);
}

#else

#include <gtest/gtest.h>
TEST(UKSSolve, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — UKS tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
