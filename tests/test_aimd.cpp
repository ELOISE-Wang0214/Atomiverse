/// @file test_aimd.cpp
/// @brief Integration tests for BOMD and second-gen CPMD via RHFForce/UHFForce.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/scf/rhf_force.h"
#include "atomiverse/cpu/scf/uhf_force.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

using namespace ATV;

// ── Helper ───────────────────────────────────────────────────────────

static cpu::Frame make_h2_with_velocity(Real bond_length, Real v_stretch = 0) {
    cpu::Frame f(2);
    f.x = {0.0, bond_length};
    f.y = {0.0, 0.0};
    f.z = {0.0, 0.0};
    f.vx = {-v_stretch/2, v_stretch/2};
    f.vy = {0.0, 0.0};
    f.vz = {0.0, 0.0};
    f.fx = {0.0, 0.0};
    f.fy = {0.0, 0.0};
    f.fz = {0.0, 0.0};
    f.mass = {1.008, 1.008};
    f.type_id = {1, 1};
    f.energy = 0;
    f.kinetic = 0;
    f.box = Mat3::diagonal(100.0, 100.0, 100.0);
    return f;
}

// ── RHFForce single-point ────────────────────────────────────────────

TEST(RHFForce, SinglePointH2) {
    auto frame = make_h2_with_velocity(0.74);
    cpu::scf::RHFForce rhf{Real(1e-8), 100};
    rhf.compute_forces(frame);

    // Energy should be close to -1.117 Hartree * 27.21 = -30.4 eV
    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);

    // Forces should be small at equilibrium
    EXPECT_NEAR(static_cast<double>(frame.fx(0)), -static_cast<double>(frame.fx(1)), 1e-6);
}

// ── BOMD energy conservation ─────────────────────────────────────────

TEST(BOMD, H2EnergyConservation) {
    // Start H2 near equilibrium with small stretching velocity
    auto frame = make_h2_with_velocity(0.74, 0.001);  // 0.001 A/fs stretch

    cpu::scf::RHFForce rhf{Real(1e-8), 100};

    std::vector<double> total_energies;
    rhf.on_step([&](const cpu::scf::AIMDStepRecord& rec) {
        total_energies.push_back(static_cast<double>(rec.total_energy));
    });

    rhf.run(frame, 5, Real(0.5));  // 5 steps, 0.5 fs

    ASSERT_GE(total_energies.size(), 5u);

    // Check energy drift: max - min should be small
    double E_min = *std::min_element(total_energies.begin(), total_energies.end());
    double E_max = *std::max_element(total_energies.begin(), total_energies.end());
    double drift = E_max - E_min;

    EXPECT_LT(drift, 0.1)  // < 0.1 eV drift over 5 steps
        << "Energy drift = " << drift << " eV over 5 BOMD steps";
}

// ── Bond oscillation ─────────────────────────────────────────────────

TEST(BOMD, H2BondOscillation) {
    auto frame = make_h2_with_velocity(0.80, 0.0);  // slightly stretched, zero velocity

    cpu::scf::RHFForce rhf{Real(1e-8), 100};
    rhf.run(frame, 5, Real(0.5));

    // After some steps, bond should have contracted (atoms attracted)
    Real bond = frame.x(1) - frame.x(0);
    EXPECT_LT(static_cast<double>(bond), 0.80)
        << "H2 bond should contract from stretched state";
    EXPECT_GT(static_cast<double>(bond), 0.50)
        << "Bond shouldn't collapse";
}

// ── Momentum conservation ────────────────────────────────────────────

TEST(BOMD, MomentumConservation) {
    auto frame = make_h2_with_velocity(0.74, 0.001);

    Real px0 = frame.mass(0) * frame.vx(0) + frame.mass(1) * frame.vx(1);
    Real py0 = frame.mass(0) * frame.vy(0) + frame.mass(1) * frame.vy(1);
    Real pz0 = frame.mass(0) * frame.vz(0) + frame.mass(1) * frame.vz(1);

    cpu::scf::RHFForce rhf{Real(1e-8), 100};
    rhf.run(frame, 5, Real(0.5));

    Real px1 = frame.mass(0) * frame.vx(0) + frame.mass(1) * frame.vx(1);
    Real py1 = frame.mass(0) * frame.vy(0) + frame.mass(1) * frame.vy(1);
    Real pz1 = frame.mass(0) * frame.vz(0) + frame.mass(1) * frame.vz(1);

    EXPECT_NEAR(static_cast<double>(px1), static_cast<double>(px0), 1e-6);
    EXPECT_NEAR(static_cast<double>(py1), static_cast<double>(py0), 1e-6);
    EXPECT_NEAR(static_cast<double>(pz1), static_cast<double>(pz0), 1e-6);
}

// ── ASPC Predictor ───────────────────────────────────────────────────

TEST(ASPCPredictor, BasicFunctionality) {
    cpu::scf::ASPCPredictor pred(3);

    // Need K+1 = 4 history points before prediction works
    xt::xtensor<Real, 2> P1 = {{1.0, 0.0}, {0.0, 1.0}};
    xt::xtensor<Real, 2> P2 = {{1.1, 0.0}, {0.0, 1.1}};
    xt::xtensor<Real, 2> P3 = {{1.2, 0.0}, {0.0, 1.2}};
    xt::xtensor<Real, 2> P4 = {{1.3, 0.0}, {0.0, 1.3}};

    xt::xtensor<Real, 2> P_pred;

    pred.record(P1);
    EXPECT_FALSE(pred.predict(P_pred));

    pred.record(P2);
    EXPECT_FALSE(pred.predict(P_pred));

    pred.record(P3);
    EXPECT_FALSE(pred.predict(P_pred));

    pred.record(P4);
    EXPECT_TRUE(pred.predict(P_pred));

    // For linear sequence 1.0, 1.1, 1.2, 1.3, extrapolation should give ~1.4
    EXPECT_NEAR(static_cast<double>(P_pred(0, 0)), 1.4, 0.01);
}

TEST(ASPCPredictor, CorrectorBlending) {
    cpu::scf::ASPCPredictor pred(3);

    xt::xtensor<Real, 2> P_scf = {{2.0}};
    xt::xtensor<Real, 2> P_pred = {{1.8}};
    xt::xtensor<Real, 2> P_corr;

    pred.correct(P_scf, P_pred, P_corr);

    // (K+2)/(2K+3) * 2.0 + (K+1)/(2K+3) * 1.8  with K=3
    // = 5/9 * 2.0 + 4/9 * 1.8 = 1.111 + 0.800 = 1.911
    EXPECT_NEAR(static_cast<double>(P_corr(0, 0)), 5.0/9.0 * 2.0 + 4.0/9.0 * 1.8, 1e-8);
}

// ── UHFForce basic test ──────────────────────────────────────────────

TEST(UHFForce, SinglePointH2) {
    auto frame = make_h2_with_velocity(0.74);
    cpu::scf::UHFForce uhf{1, 1, Real(1e-8), 100};
    uhf.compute_forces(frame);

    // Energy should be in the right range
    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);

    // Forces should be equal and opposite
    EXPECT_NEAR(static_cast<double>(frame.fx(0)), -static_cast<double>(frame.fx(1)), 1e-4);
}

TEST(UHFForce, H2BOMD) {
    auto frame = make_h2_with_velocity(0.80, 0.0);
    cpu::scf::UHFForce uhf{1, 1, Real(1e-8), 100};
    uhf.run(frame, 3, Real(0.5));

    Real bond = frame.x(1) - frame.x(0);
    EXPECT_LT(static_cast<double>(bond), 0.80)
        << "UHF H2 bond should contract from stretched state";
}

// ── CPMD2 mode ───────────────────────────────────────────────────────

TEST(CPMD2, RHFForceWithPredictor) {
    auto frame = make_h2_with_velocity(0.74, 0.001);

    cpu::scf::RHFForce rhf{Real(1e-6), 50};
    rhf.predictor(3, 2, Real(1e-4), true);

    // Should run without crashing
    rhf.run(frame, 5, Real(0.5));

    // Energy should be reasonable
    EXPECT_LT(static_cast<double>(frame.energy), -29.0);
    EXPECT_GT(static_cast<double>(frame.energy), -32.0);

    // Trace should exist
    EXPECT_EQ(rhf.trace().size(), 5u);
}
