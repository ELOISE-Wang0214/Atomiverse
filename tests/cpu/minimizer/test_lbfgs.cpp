/// @file tests/cpu/minimizer/test_lbfgs.cpp
/// @brief Tests for CPU L-BFGS geometry optimizer with SCF forces.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/minimizer/lbfgs.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"

using namespace ATV;

class LBFGSTest : public ::testing::Test {
protected:
    static cpu::Frame make_h2(double bond_length) {
        cpu::Frame f(2);
        f.atoms.x = {0.0, bond_length};
        f.atoms.y = {0.0, 0.0};
        f.atoms.z = {0.0, 0.0};
        f.atoms.mass = {1.008, 1.008};
        f.atoms.type = {1, 1};
        f.box = Mat3::cubic(100.0);
        return f;
    }
};

TEST_F(LBFGSTest, H2RHFGeometryOptimization) {
    auto frame = make_h2(1.2);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, double(0.01), double(0.5), 10});
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.build();
    opt.run(frame, 50);

    double bond = frame.atoms.x(1) - frame.atoms.x(0);
    EXPECT_LT(static_cast<double>(bond), 1.2)
        << "H2 bond should contract from 1.2 Å";
    EXPECT_GT(static_cast<double>(bond), 0.5)
        << "H2 bond should not collapse";
}

TEST_F(LBFGSTest, H2LSDAGeometryOptimization) {
    auto frame = make_h2(1.2);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uks_config(integrals, 1, 1, cpu::scf::XCType::LSDA);

    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, double(0.01), double(0.5), 10});
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.build();
    opt.run(frame, 50);

    double bond = frame.atoms.x(1) - frame.atoms.x(0);
    EXPECT_LT(static_cast<double>(bond), 1.2)
        << "H2 bond should contract from 1.2 Å (LSDA)";
    EXPECT_GT(static_cast<double>(bond), 0.5);
}

TEST_F(LBFGSTest, CompressedH2Expands) {
    auto frame = make_h2(0.5);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::Minimize opt;
    opt.add(cpu::minimizer::LBFGS{5, double(0.01), double(0.5), 10});
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.build();
    opt.run(frame, 30);

    double bond = frame.atoms.x(1) - frame.atoms.x(0);
    EXPECT_GT(static_cast<double>(bond), 0.5)
        << "Compressed H2 should expand";
}

TEST_F(LBFGSTest, H2RHFArmijoLineSearchReachesTightForce) {
    // Armijo-backtracking line search should drive max-force below a tight
    // threshold that the fixed-step default cannot reliably hit. Uses the
    // same H2/STO-3G setup but a much tighter ftol (1 mHa/Bohr ≈ 0.025 eV/Å).
    auto frame = make_h2(1.2);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::LBFGS step{
        /*memory=*/5, /*alpha0=*/0.05,
        /*ftol=*/0.001,       // eV/Å
        /*max_ls=*/20,
        /*rms_ftol=*/0.0005,
        /*max_dtol=*/0.001,
        /*rms_dtol=*/0.0005,
        /*etol=*/1e-7};
    step.line_search();  // Armijo default constants.
    cpu::minimizer::Minimize opt;
    opt.add(std::move(step));
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.build();
    opt.run(frame, 80);

    double bond = frame.atoms.x(1) - frame.atoms.x(0);
    // H2/STO-3G PySCF-minimum bond ≈ 0.7122 Å.
    EXPECT_NEAR(static_cast<double>(bond), 0.7122, 0.01);
}
