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
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::LBFGS opt;
    opt.add(cpu::minimizer::LBFGSStep{5, double(0.01), double(0.5), 10});
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
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::uks_config(integrals, 1, 1, cpu::scf::XCType::LSDA);

    cpu::minimizer::LBFGS opt;
    opt.add(cpu::minimizer::LBFGSStep{5, double(0.01), double(0.5), 10});
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
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::minimizer::LBFGS opt;
    opt.add(cpu::minimizer::LBFGSStep{5, double(0.01), double(0.5), 10});
    opt.on(ATV::minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
    opt.build();
    opt.run(frame, 30);

    double bond = frame.atoms.x(1) - frame.atoms.x(0);
    EXPECT_GT(static_cast<double>(bond), 0.5)
        << "Compressed H2 should expand";
}
