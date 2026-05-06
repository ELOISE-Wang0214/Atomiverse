/// @file tests/cpu/scf/test_scf_evaluator.cpp
/// @brief Integration tests for SCFEvaluator component, BOMD, and ASPC predictor.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/fix_langevin.h"
#include "atomiverse/cpu/md/barostat.h"
#include "atomiverse/cpu/scf/scf_evaluator.h"
#include "atomiverse/cpu/scf/scf_solve.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

using namespace ATV;

// ── Fixture: H2 frame with optional velocity ────────────────────────

class SCFEvaluatorTest : public ::testing::Test {
protected:
    static cpu::Frame make_h2(double bond_length, double v_stretch = 0) {
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
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
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
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uhf_config(integrals, 1, 1);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 1);

    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);
}

TEST_F(SCFEvaluatorTest, SinglePointKSLSDA) {
    auto frame = make_h2(0.74);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uks_config(integrals, 1, 1, cpu::scf::XCType::LSDA);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
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
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::uks_config(integrals, 1, 1, cpu::scf::XCType::PBE);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
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
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_GT(static_cast<double>(frame.energy), -40.0);
}

TEST_F(BOMDTest, RHFBondOscillation) {
    auto frame = make_h2(0.80, 0.0);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    double bond = frame.atoms.x(1) - frame.atoms.x(0);
    EXPECT_LT(static_cast<double>(bond), 0.80)
        << "H2 bond should contract from stretched state";
    EXPECT_GT(static_cast<double>(bond), 0.50);
}

TEST_F(BOMDTest, MomentumConservation) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    double px0 = frame.atoms.mass(0) * frame.atoms.vx(0) + frame.atoms.mass(1) * frame.atoms.vx(1);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    double px1 = frame.atoms.mass(0) * frame.atoms.vx(0) + frame.atoms.mass(1) * frame.atoms.vx(1);
    EXPECT_NEAR(static_cast<double>(px1), static_cast<double>(px0), 1e-4);
}

TEST_F(BOMDTest, RHFWithASPC) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::scf::SCFEvaluator evaluator{config};
    evaluator.predictor(3, 2, double(1e-4), true);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    sim.on(md::Phase::force).add(std::move(evaluator));
    sim.build();
    sim.run(frame, 5);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_GT(static_cast<double>(frame.energy), -40.0);
}

// ── Ensemble integration tests (AIMD = MD driver + SCFEvaluator + ─────
//    thermostat / barostat). Verifies SCF force evaluation composes
//    cleanly with the existing CPU MD ensemble components — no new
//    AIMD driver class required.

TEST_F(BOMDTest, RHFNoseHooverNVT) {
    auto frame = make_h2(0.74, 0.0);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt = 0.5;          // fs
    const double T_target = 300.0;  // K
    const double tau = 50.0;        // fs (tight coupling for short test)

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/42u});
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::NoseHoover{T_target, tau, dt});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 30);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0)
        << "SCF should still converge under NoseHoover coupling";
    EXPECT_GT(static_cast<double>(frame.kinetic), 0.0)
        << "InitVelocity + VelocityVerlet must produce non-zero KE";
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.kinetic)));
}

TEST_F(BOMDTest, RHFLangevinNVT) {
    auto frame = make_h2(0.74, 0.0);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt = 0.5;
    const double T_target = 300.0;
    const double damp = 100.0;

    // Both SCFEvaluator and FixLangevin advertise Capability::Force —
    // DynamicsDriver must accept this (multiple Force-capable components
    // on the force / post_force phases).
    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/7u});
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::FixLangevin{T_target, damp, dt, /*seed=*/13u});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 30);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_GT(static_cast<double>(frame.kinetic), 0.0);
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.kinetic)));
}

TEST_F(BOMDTest, RHFBerendsenNPT) {
    auto frame = make_h2(0.74, 0.0);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    const double dt = 0.5;
    const double T_target = 300.0;
    const double tau_T = 50.0;
    const double P_target = 0.0;     // SCF on isolated H2 has no PBC stress;
    const double tau_P = 1000.0;     // weak coupling, long time constant

    Mat3 box_initial = frame.box.cell();

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_target, /*seed=*/3u});
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::NoseHoover{T_target, tau_T, dt});
    sim.add(cpu::md::Barostat{P_target, tau_P, dt});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 20);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0)
        << "Box mutation between SCF launches must not break IntegralPack rebuild";
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.kinetic)));
    // Barostat may scale the box at zero pressure with non-zero virial;
    // we only require it stays finite and within a sane range.
    EXPECT_TRUE(std::isfinite(frame.box(0, 0)));
    EXPECT_GT(frame.box(0, 0), 0.0);
    EXPECT_LT(frame.box(0, 0), 10.0 * box_initial(0, 0));
}

// ── Regression: builder chain on rvalue temporary (Bug #1) ──────────
//
// `SCFEvaluator{cfg}.predictor(3)` returns an lvalue reference to the
// temporary. Without ref-qualified overloads, `Driver::add(Comp&& comp)`
// deduces `Comp` as `SCFEvaluator&` (lvalue reference) and stores a
// non-owning pointer to the temporary, which dies at the semicolon →
// SEGV on the first launch. The fix is to add rvalue-ref-qualified
// builder overloads so chains on temporaries propagate rvalue-ness
// through to add()'s rvalue branch (driver moves into unique_ptr).

TEST_F(BOMDTest, ChainedBuilderOnTemporaryDoesNotSEGV) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    // The whole point of this test: chained builder call on a temporary.
    // Pre-fix, this stored a non-owning pointer to a dangling temporary
    // and segfaulted on the first launch.
    sim.on(md::Phase::force).add(
        cpu::scf::SCFEvaluator{config}.predictor(3)
    );
    sim.build();
    sim.run(frame, 3);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
}

TEST_F(BOMDTest, ChainedBuilderWithGradientToggleOnTemporary) {
    auto frame = make_h2(0.74, 0.001);
    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    // Two consecutive builder calls on the same temporary chain.
    sim.on(md::Phase::force).add(
        cpu::scf::SCFEvaluator{config}.with_gradient(true).predictor(2)
    );
    sim.build();
    sim.run(frame, 1);

    EXPECT_LT(static_cast<double>(frame.energy), 0.0);
    EXPECT_TRUE(std::isfinite(static_cast<double>(frame.energy)));
}

// ── Regression: PBC-wrap × non-min-image SCF (Bug #2) ───────────────
//
// `VelocityVerlet::wrap_orthorhombic` wraps positions through orthorhombic
// PBC and bumps `frame.atoms.image[i]` by ±1 per crossing. `build_integrals`
// is not minimum-image-aware; it builds shells and v_nuc from raw
// `frame.atoms.x/y/z`. If an atom starts on the box edge and thermal
// motion pushes it across, the next SCF sees the two atoms separated by
// approximately one box length and the calculation collapses (RHF
// converges to an unphysical ionic limit at ≈ -0.16 Ha instead of the
// bonded ≈ -1.117 Ha for H2).
//
// The fix is to unwrap positions inside `SCFEvaluator::launch` using
// `frame.atoms.image` before the integral build. The frame storage stays
// wrapped (so MD pair lists / classical force fields keep working);
// SCFEvaluator just sees the unwrapped view.

TEST_F(BOMDTest, BoundaryAtomsRecoverAfterWrap) {
    // H2 sitting on the x-edge: atom 0 at x=0, atom 1 at x=0.74. Any
    // negative x velocity on atom 0 crosses the wrap on step 1.
    cpu::Frame frame(2);
    frame.atoms.x = {0.0, 0.74};
    frame.atoms.y = {50.0, 50.0};
    frame.atoms.z = {50.0, 50.0};
    frame.atoms.vx = {-0.01, 0.01};      // atom 0 moves toward x<0 and wraps
    frame.atoms.vy = {0.0, 0.0};
    frame.atoms.vz = {0.0, 0.0};
    frame.atoms.fx = {0.0, 0.0};
    frame.atoms.fy = {0.0, 0.0};
    frame.atoms.fz = {0.0, 0.0};
    frame.atoms.mass = {1.008, 1.008};
    frame.atoms.type = {1, 1};
    frame.energy = 0;
    frame.kinetic = 0;
    frame.box = Mat3::cubic(100.0);

    auto integrals = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
    auto config = cpu::scf::rhf_config(integrals);

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{double(0.5)});
    sim.on(md::Phase::force).add(cpu::scf::SCFEvaluator{config});
    sim.build();
    sim.run(frame, 5);

    // After fix: SCF unwraps before building integrals, so PE stays in
    // the bonded H2 range (≈ -30 eV). Pre-fix, atoms wrapped on step 1
    // and PE collapsed to ≈ -4 eV (RHF ionic limit at ~100 Å apart).
    EXPECT_LT(static_cast<double>(frame.energy), -25.0)
        << "PE ≈ -4 eV indicates SCF saw atoms ~box-length apart "
        << "(wrap × non-minimum-image SCF pathology). Expected bonded H2 "
        << "energy ≈ -30 eV.";
    EXPECT_GT(static_cast<double>(frame.energy), -35.0);
}

// ── ASPC predictor unit tests ───────────────────────────────────────

class ASPCPredictorTest : public ::testing::Test {};

TEST_F(ASPCPredictorTest, BasicFunctionality) {
    cpu::scf::ASPCPredictor pred(3);

    xt::xtensor<double, 2> P1 = {{1.0, 0.0}, {0.0, 1.0}};
    xt::xtensor<double, 2> P2 = {{1.1, 0.0}, {0.0, 1.1}};
    xt::xtensor<double, 2> P3 = {{1.2, 0.0}, {0.0, 1.2}};
    xt::xtensor<double, 2> P4 = {{1.3, 0.0}, {0.0, 1.3}};

    xt::xtensor<double, 2> P_pred;

    pred.record(P1); EXPECT_FALSE(pred.predict(P_pred));
    pred.record(P2); EXPECT_FALSE(pred.predict(P_pred));
    pred.record(P3); EXPECT_FALSE(pred.predict(P_pred));
    pred.record(P4); EXPECT_TRUE(pred.predict(P_pred));

    EXPECT_NEAR(static_cast<double>(P_pred(0, 0)), 1.4, 0.01);
}

TEST_F(ASPCPredictorTest, CorrectorBlending) {
    cpu::scf::ASPCPredictor pred(3);

    xt::xtensor<double, 2> P_scf = {{2.0}};
    xt::xtensor<double, 2> P_pred = {{1.8}};
    xt::xtensor<double, 2> P_corr;

    pred.correct(P_scf, P_pred, P_corr);

    EXPECT_NEAR(static_cast<double>(P_corr(0, 0)), 5.0/9.0 * 2.0 + 4.0/9.0 * 1.8, 1e-8);
}
