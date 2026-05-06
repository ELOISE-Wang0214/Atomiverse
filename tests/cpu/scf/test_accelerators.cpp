/// @file tests/cpu/scf/test_accelerators.cpp
/// @brief Isolated unit tests for DIIS extrapolator and ASPC predictor.
///
/// These SCF accelerators are currently only tested transitively through
/// scf_solve convergence; this file adds direct, isolated verification
/// of their internal mechanics.

#ifdef ATV_HAS_RHF

#include <gtest/gtest.h>
#include <cmath>

#include <xtensor/containers/xtensor.hpp>

#include "atomiverse/cpu/scf/internal/diis_mixer.h"
#include "atomiverse/cpu/scf/aspc_predictor.h"

using namespace ATV::cpu::scf;
using ATV::cpu::scf::internal::DIISMixer;

// ── DIIS (internal::DIISMixer) ───────────────────────────────────────

TEST(DIIS, SingleEntry_ReturnsSameMatrix) {
    DIISMixer diis(6);
    RealMatrix F = {{1.0, 0.5}, {0.5, 2.0}};
    RealMatrix e = {{0.01, 0.0}, {0.0, -0.01}};
    auto result = diis.step(F, e);
    EXPECT_NEAR(result(0, 0), 1.0, 1e-12);
    EXPECT_NEAR(result(0, 1), 0.5, 1e-12);
    EXPECT_NEAR(result(1, 1), 2.0, 1e-12);
}

TEST(DIIS, TwoEntries_WeightsSmallError) {
    DIISMixer diis(6);
    RealMatrix F1 = {{1.0, 0.0}, {0.0, 2.0}};
    RealMatrix e1 = {{0.1, 0.0}, {0.0, -0.1}};
    RealMatrix F2 = {{1.1, 0.0}, {0.0, 1.9}};
    RealMatrix e2 = {{0.01, 0.0}, {0.0, -0.01}};
    diis.step(F1, e1);  // push first pair
    auto result = diis.step(F2, e2);  // push second + extrapolate
    EXPECT_NE(result(0, 0), F1(0, 0));
    EXPECT_NE(result(0, 0), F2(0, 0));
}

TEST(DIIS, DualChannel_BothExtracted) {
    DIISMixer diis(6);
    RealMatrix Fa = {{1.0}}, Fb = {{2.0}};
    RealMatrix ea = {{0.01}}, eb = {{0.02}};
    auto [ra, rb] = diis.step(Fa, ea, Fb, eb);
    EXPECT_NEAR(ra(0, 0), 1.0, 1e-12);
    EXPECT_NEAR(rb(0, 0), 2.0, 1e-12);
}

TEST(DIIS, ResetClearsHistory) {
    DIISMixer diis(6);
    RealMatrix F = {{1.0}}, e = {{0.01}};
    diis.step(F, e);
    EXPECT_EQ(diis.size(), 1);
    diis.reset();
    EXPECT_EQ(diis.size(), 0);
}

// ── ASPC ─────────────────────────────────────────────────────────────

TEST(ASPC, PredictReturnsFalseBeforeHistory) {
    ASPCPredictor aspc(3);
    xt::xtensor<double, 2> pred = xt::zeros<double>({2, 2});
    bool ok = aspc.predict(pred);
    EXPECT_FALSE(ok) << "predict() should return false with empty history";
}

TEST(ASPC, RecordAndPredict) {
    ASPCPredictor aspc(2);  // order 2 needs 3 history points
    xt::xtensor<double, 2> P1 = {{1.0, 0.0}, {0.0, 1.0}};
    xt::xtensor<double, 2> P2 = {{1.1, 0.0}, {0.0, 1.1}};
    xt::xtensor<double, 2> P3 = {{1.2, 0.0}, {0.0, 1.2}};
    aspc.record(P1);
    aspc.record(P2);
    aspc.record(P3);
    xt::xtensor<double, 2> pred = xt::zeros<double>({2, 2});
    bool ok = aspc.predict(pred);
    EXPECT_TRUE(ok) << "predict() should succeed with 3 history points at order 2";
    // Linear trend: each step adds 0.1 → prediction should be ≥ 1.2
    EXPECT_GT(pred(0, 0), 1.15)
        << "ASPC should extrapolate past last value for linear trend";
}

TEST(ASPC, CorrectBlends) {
    ASPCPredictor aspc(3);
    // K=3: alpha_scf = (K+2)/(2K+3) = 5/9, alpha_pred = (K+1)/(2K+3) = 4/9
    xt::xtensor<double, 2> P_scf  = {{1.0, 0.0}, {0.0, 1.0}};
    xt::xtensor<double, 2> P_pred = {{2.0, 0.0}, {0.0, 2.0}};
    xt::xtensor<double, 2> P_out  = xt::zeros<double>({2, 2});
    aspc.correct(P_scf, P_pred, P_out);
    double expected = (5.0/9.0) * 1.0 + (4.0/9.0) * 2.0;
    EXPECT_NEAR(P_out(0, 0), expected, 1e-12);
}

TEST(ASPC, ResetClearsHistory) {
    ASPCPredictor aspc(3);
    xt::xtensor<double, 2> P = {{5.0, 0.0}, {0.0, 5.0}};
    aspc.record(P);
    EXPECT_EQ(aspc.history_size(), 1);
    aspc.reset();
    EXPECT_EQ(aspc.history_size(), 0);
}

#endif // ATV_HAS_RHF
