/// @file tests/cpu/scf/test_scf_evaluator.cpp
/// @brief Integration tests for SCFEvaluator component, BOMD, and ASPC predictor.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

using namespace ATV;

// ── Fixture: H2 frame with optional velocity ────────────────────────

class SCFEvaluatorTest : public ::testing::Test {
protected:
    static cpu::Frame make_h2(Real bond_length, Real v_stretch = 0) {
        cpu::Frame f(2);
        f.atoms.x = {0.0, bond_length};
        f.atoms.y = {0.0, 0.0};
        f.atoms.z = {0.0, 0.0};
        f.atoms.vx = {-v_stretch/2, v_stretch/2};
        f.atoms.vy = {0.0, 0.0};
        f.atoms.vz = {0.0, 0.0};
        f.atoms.fx = {0.0, 0.0};
        f.atoms.fy = {0.0, 0.0};
        f.atoms.fz = {0.0, 0.0};
        f.atoms.mass = {1.008, 1.008};
        f.atoms.type = {1, 1};
        f.energy = 0;
        f.kinetic = 0;
        f.box = Mat3::cubic(100.0);
        return f;
    }
};

TEST_F(SCFEvaluatorTest, SinglePointRHF) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 1);

    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);
    EXPECT_NEAR(static_cast<double>(frame.atoms.fx(0)),
                -static_cast<double>(frame.atoms.fx(1)), 0.01);
}

TEST_F(SCFEvaluatorTest, SinglePointUHF) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::uhf_config(integrals, 1, 1);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 1);

    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);
}

TEST_F(SCFEvaluatorTest, SinglePointKSLSDA) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::uks_config(integrals, 1, 1, cpu::scf::XCType::LSDA);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 1);

    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);
    EXPECT_NEAR(static_cast<double>(frame.atoms.fx(0)),
                -static_cast<double>(frame.atoms.fx(1)), 0.02);
}

TEST_F(SCFEvaluatorTest, SinglePointKSPBE) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::uks_config(integrals, 1, 1, cpu::scf::XCType::PBE);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 1);

    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -33.0);
}

// ── BOMD integration tests ──────────────────────────────────────────

class BOMDTest : public SCFEvaluatorTest {};

TEST_F(BOMDTest, RHFEnergyConservation) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_GT(static_cast<double>(frame.energy), -40.0);
}

TEST_F(BOMDTest, RHFBondOscillation) {
    auto frame = make_h2(0.80, 0.0);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    Real bond = frame.atoms.x(1) - frame.atoms.x(0);
    EXPECT_LT(static_cast<double>(bond), 0.80)
        << "H2 bond should contract from stretched state";
    EXPECT_GT(static_cast<double>(bond), 0.50);
}

TEST_F(BOMDTest, MomentumConservation) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    Real px0 = frame.atoms.mass(0) * frame.atoms.vx(0) + frame.atoms.mass(1) * frame.atoms.vx(1);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    Real px1 = frame.atoms.mass(0) * frame.atoms.vx(0) + frame.atoms.mass(1) * frame.atoms.vx(1);
    EXPECT_NEAR(static_cast<double>(px1), static_cast<double>(px0), 1e-4);
}

TEST_F(BOMDTest, RHFWithASPC) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals_sto3g(frame);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFEvaluator evaluator{config};
    evaluator.predictor(3, 2, Real(1e-4), true);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{Real(0.5)});
    sim.on(md::Phase::force).add(std::move(evaluator));
    sim.build();
    sim.run(frame, 5);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_GT(static_cast<double>(frame.energy), -40.0);
}

// ── ASPC predictor unit tests ───────────────────────────────────────

class ASPCPredictorTest : public ::testing::Test {};

TEST_F(ASPCPredictorTest, BasicFunctionality) {
    cpu::scf::ASPCPredictor pred(3);

    xt::xtensor<Real, 2> P1 = {{1.0, 0.0}, {0.0, 1.0}};
    xt::xtensor<Real, 2> P2 = {{1.1, 0.0}, {0.0, 1.1}};
    xt::xtensor<Real, 2> P3 = {{1.2, 0.0}, {0.0, 1.2}};
    xt::xtensor<Real, 2> P4 = {{1.3, 0.0}, {0.0, 1.3}};

    xt::xtensor<Real, 2> P_pred;

    pred.record(P1); EXPECT_FALSE(pred.predict(P_pred));
    pred.record(P2); EXPECT_FALSE(pred.predict(P_pred));
    pred.record(P3); EXPECT_FALSE(pred.predict(P_pred));
    pred.record(P4); EXPECT_TRUE(pred.predict(P_pred));

    EXPECT_NEAR(static_cast<double>(P_pred(0, 0)), 1.4, 0.01);
}

TEST_F(ASPCPredictorTest, CorrectorBlending) {
    cpu::scf::ASPCPredictor pred(3);

    xt::xtensor<Real, 2> P_scf = {{2.0}};
    xt::xtensor<Real, 2> P_pred = {{1.8}};
    xt::xtensor<Real, 2> P_corr;

    pred.correct(P_scf, P_pred, P_corr);

    EXPECT_NEAR(static_cast<double>(P_corr(0, 0)), 5.0/9.0 * 2.0 + 4.0/9.0 * 1.8, 1e-8);
}
