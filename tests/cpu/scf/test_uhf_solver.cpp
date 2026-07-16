/// @file tests/cpu/scf/test_uhf_solver.cpp
/// @brief Tests for UHF configurations through the unified SCFSolver.

#ifdef ATV_HAS_RHF

#include <cmath>
#include <gtest/gtest.h>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/reducers/xnorm.hpp>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/integrals.h"
#include "atomiverse/cpu/scf/scf_solver.h"

using namespace ATV;

class UHFSolverTest : public ::testing::Test {
  protected:
    static constexpr double kTol =
        (sizeof(double) == sizeof(float)) ? double(1e-5) : double(1e-9);
    static constexpr double kEnergyTol = 1e-3;

    static cpu::Frame make_h_atom() {
        cpu::Frame f(1);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_type)[0] = 1;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }

    static cpu::Frame make_he_atom() {
        cpu::Frame f(1);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_type)[0] = 2;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }

    static cpu::Frame make_li_atom() {
        cpu::Frame f(1);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_type)[0] = 3;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }

    static cpu::Frame make_h2(double bond_length_angstrom) {
        cpu::Frame f(2);
        f.span(ATV::kAtom_x)[0] = double(0);
        f.span(ATV::kAtom_y)[0] = double(0);
        f.span(ATV::kAtom_z)[0] = double(0);
        f.span(ATV::kAtom_x)[1] = double(bond_length_angstrom);
        f.span(ATV::kAtom_y)[1] = double(0);
        f.span(ATV::kAtom_z)[1] = double(0);
        f.span(ATV::kAtom_type)[0] = 1;
        f.span(ATV::kAtom_type)[1] = 1;
        f.simbox = ATV::cpu::Box::cubic(100.0);
        return f;
    }
};

// ── H atom ──────────────────────────────────────────────────────────

TEST_F(UHFSolverTest, HAtomConverges) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 0, kTol);
    config.max_iter = 50;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 1);
    EXPECT_EQ(result.n_beta, 0);
}

TEST_F(UHFSolverTest, HAtomEnergy) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 0, kTol);
    config.max_iter = 50;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_NEAR(static_cast<double>(result.energy), -0.4666, 5e-3);
}

TEST_F(UHFSolverTest, HAtomElectronCount) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 0, kTol);
    config.max_iter = 50;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    double tr_alpha = xt::sum(result.P_alpha * ints.S)();
    double tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 1.0, 1e-5);
    EXPECT_NEAR(static_cast<double>(tr_beta), 0.0, 1e-5);
}

TEST_F(UHFSolverTest, HAtomSpinContamination) {
    auto frame = make_h_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 0, kTol);
    config.max_iter = 50;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_NEAR(static_cast<double>(result.spin_contamination), 0.0, 1e-6);
}

// ── He atom ─────────────────────────────────────────────────────────

TEST_F(UHFSolverTest, HeAtomEqualsRHF) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto rhf_config = cpu::scf::rhf_config(ints, kTol);
    rhf_config.max_iter = 50;
    cpu::scf::SCFSolver rhf_solver{rhf_config};
    rhf_solver.build(ints);
    const auto& rhf = rhf_solver.run(ints);
    auto uhf_config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    uhf_config.max_iter = 50;
    cpu::scf::SCFSolver uhf_solver{uhf_config};
    uhf_solver.build(ints);
    const auto& uhf = uhf_solver.run(ints);
    EXPECT_TRUE(uhf.converged);
    EXPECT_NEAR(static_cast<double>(uhf.energy),
                static_cast<double>(rhf.energy), kEnergyTol);
}

TEST_F(UHFSolverTest, HeAtomSpinContamination) {
    auto frame = make_he_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    config.max_iter = 50;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& uhf = solver.run(ints);
    EXPECT_NEAR(static_cast<double>(uhf.spin_contamination), 0.0, 1e-4);
}

// ── Li atom ─────────────────────────────────────────────────────────

TEST_F(UHFSolverTest, LiAtomConverges) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 2, 1, kTol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n_alpha, 2);
    EXPECT_EQ(result.n_beta, 1);
}

TEST_F(UHFSolverTest, LiAtomElectronCount) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 2, 1, kTol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    double tr_alpha = xt::sum(result.P_alpha * ints.S)();
    double tr_beta = xt::sum(result.P_beta * ints.S)();
    EXPECT_NEAR(static_cast<double>(tr_alpha), 2.0, 1e-4);
    EXPECT_NEAR(static_cast<double>(tr_beta), 1.0, 1e-4);
}

TEST_F(UHFSolverTest, LiAtomEnergyRange) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 2, 1, kTol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_TRUE(result.converged);
    EXPECT_LT(static_cast<double>(result.energy), -7.0);
    EXPECT_GT(static_cast<double>(result.energy), -8.0);
}

// ── H2 equilibrium ──────────────────────────────────────────────────

TEST_F(UHFSolverTest, H2EquilibriumAgreesWithRHF) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto rhf_config = cpu::scf::rhf_config(ints, kTol);
    rhf_config.max_iter = 200;
    cpu::scf::SCFSolver rhf_solver{rhf_config};
    rhf_solver.build(ints);
    const auto& rhf = rhf_solver.run(ints);
    auto uhf_config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    uhf_config.max_iter = 200;
    cpu::scf::SCFSolver uhf_solver{uhf_config};
    uhf_solver.build(ints);
    const auto& uhf = uhf_solver.run(ints);
    EXPECT_TRUE(uhf.converged);
    EXPECT_NEAR(static_cast<double>(uhf.energy),
                static_cast<double>(rhf.energy), kEnergyTol);
}

TEST_F(UHFSolverTest, H2EquilibriumSpinContamination) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    config.max_iter = 200;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& uhf = solver.run(ints);
    EXPECT_NEAR(static_cast<double>(uhf.spin_contamination), 0.0, 1e-3);
}

// ── H2 stretched ────────────────────────────────────────────────────

TEST_F(UHFSolverTest, H2StretchedUHFLowerThanRHF) {
    auto frame = make_h2(2.0);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto rhf_config = cpu::scf::rhf_config(ints, kTol);
    rhf_config.max_iter = 50;
    cpu::scf::SCFSolver rhf_solver{rhf_config};
    rhf_solver.build(ints);
    const auto& rhf = rhf_solver.run(ints);
    auto uhf_config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    uhf_config.max_iter = 100;
    cpu::scf::SCFSolver uhf_solver{uhf_config};
    uhf_solver.build(ints);
    const auto& uhf = uhf_solver.run(ints);
    EXPECT_TRUE(uhf.converged);
    EXPECT_LE(static_cast<double>(uhf.energy),
              static_cast<double>(rhf.energy) + 1e-6);
}

TEST_F(UHFSolverTest, H2StretchedSpinContaminationNonzero) {
    auto frame = make_h2(2.0);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& uhf = solver.run(ints);
    EXPECT_TRUE(uhf.converged);
    EXPECT_GT(static_cast<double>(uhf.spin_contamination), 0.0);
}

// ── PES scan ────────────────────────────────────────────────────────

TEST_F(UHFSolverTest, H2PESScanUHFVariational) {
    for (double r : {0.6, 0.74, 1.0, 1.5, 2.0, 2.5, 3.0}) {
        auto frame = make_h2(r);
        auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
        auto rhf_config = cpu::scf::rhf_config(ints, kTol);
        rhf_config.max_iter = 200;
        cpu::scf::SCFSolver rhf_solver{rhf_config};
        rhf_solver.build(ints);
        const auto& rhf = rhf_solver.run(ints);
        auto uhf_config = cpu::scf::uhf_config(ints, 1, 1, kTol);
        uhf_config.max_iter = 400;
        cpu::scf::SCFSolver uhf_solver{uhf_config};
        uhf_solver.build(ints);
        const auto& uhf = uhf_solver.run(ints);
        EXPECT_TRUE(uhf.converged) << "UHF failed to converge at R=" << r;
        EXPECT_LE(static_cast<double>(uhf.energy),
                  static_cast<double>(rhf.energy) + 1e-5)
            << "UHF > RHF at R=" << r;
    }
}

// ── Result structure ────────────────────────────────────────────────

TEST_F(UHFSolverTest, ResultContainsAlphaBetaEigenvalues) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 2, 1, kTol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_EQ(result.eigenvalues_alpha.size(), (size_t)ints.n_basis);
    EXPECT_EQ(result.eigenvalues_beta.size(), (size_t)ints.n_basis);
}

TEST_F(UHFSolverTest, ResultHasConvergenceTrace) {
    auto frame = make_h2(0.74);
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 1, 1, kTol);
    config.max_iter = 50;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    EXPECT_GT(result.trace.size(), 0u);
    EXPECT_EQ((int)result.trace.back().iter, result.iterations);
}

TEST_F(UHFSolverTest, DensityMatricesHaveCorrectShape) {
    auto frame = make_li_atom();
    auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(ints, 2, 1, kTol);
    config.max_iter = 100;
    cpu::scf::SCFSolver solver{config};
    solver.build(ints);
    const auto& result = solver.run(ints);
    int n = ints.n_basis;
    EXPECT_EQ(result.P_alpha.shape(0), (size_t)n);
    EXPECT_EQ(result.P_alpha.shape(1), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(0), (size_t)n);
    EXPECT_EQ(result.P_beta.shape(1), (size_t)n);
}

#else

#include <gtest/gtest.h>
TEST(UHFSolver, SkippedNoBuildRHF) {
    GTEST_SKIP() << "ATV_BUILD_RHF=OFF — UHF tests require -DATV_BUILD_RHF=ON";
}

#endif // ATV_HAS_RHF
