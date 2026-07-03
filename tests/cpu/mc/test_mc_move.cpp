/// @file tests/cpu/mc/test_mc_move.cpp
/// @brief Fast unit tests for the MC move + driver assembly
///        (mc-sampling-02-move): TranslateMove's per-launch behaviour and the
///        MC driver's build-time validation + config-constructor guards.
///
/// These are true unit tests — they do NOT drive `MC::run()` (a full
/// backend/plan pipeline, i.e. an integration/e2e path). The Mover is exercised
/// by calling `launch()` directly against a hand-built context; the driver is
/// exercised only through `build()` (plan compile + validation, no stepping).
/// End-to-end sampler behaviour (bootstrap energy, accept/reject/rollback,
/// determinism, Metropolis distribution) belongs to the validation tier in
/// sub-spec mc-sampling-03-validate, not to this ms-scale suite.

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <random>

#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/mc/mc.h"
#include "atomiverse/cpu/mc/transaction.h"
#include "atomiverse/cpu/mc/translate_move.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"

using namespace ATV;

namespace {
constexpr double kEps = 0.0104, kSigma = 3.4, kCut = 8.0, kSkin = 2.0;
} // namespace

// ── TranslateMove::launch — unit behaviour (no driver) ────────────────────

TEST(TranslateMoveUnit, DisplacesOneAtomAndRecordsChangeset) {
    cpu::Frame frame(3);
    for (int i = 0; i < 3; ++i) {
        frame.atoms.x[i] = static_cast<double>(i);
        frame.atoms.y[i] = 0.0;
        frame.atoms.z[i] = 0.0;
    }
    cpu::Frame before = frame;

    std::mt19937_64 rng(123u);
    cpu::mc::Transaction txn;
    cpu::Context ctx;
    ctx.rng = &rng;
    ctx.proposal = &txn;

    cpu::mc::TranslateMove mv{0.5};
    cpu::LaunchContext lctx{ctx, frame, /*stage=*/0, /*step=*/0};
    mv.launch(lctx);

    // Exactly one atom recorded, index in range.
    ASSERT_EQ(txn.size(), 1u);
    const int i = txn.moved[0];
    ASSERT_GE(i, 0);
    ASSERT_LT(i, 3);

    // old_* is the pre-move position; new_* is the (current) post-move position.
    EXPECT_EQ(txn.old_x[0], before.atoms.x[i]);
    EXPECT_EQ(txn.old_y[0], before.atoms.y[i]);
    EXPECT_EQ(txn.old_z[0], before.atoms.z[i]);
    EXPECT_EQ(txn.new_x[0], frame.atoms.x[i]);
    EXPECT_EQ(txn.new_y[0], frame.atoms.y[i]);
    EXPECT_EQ(txn.new_z[0], frame.atoms.z[i]);

    // Displacement is within the symmetric box [-δ, δ] per axis.
    EXPECT_LE(std::abs(txn.new_x[0] - txn.old_x[0]), 0.5);
    EXPECT_LE(std::abs(txn.new_y[0] - txn.old_y[0]), 0.5);
    EXPECT_LE(std::abs(txn.new_z[0] - txn.old_z[0]), 0.5);

    // Only atom i moved; the others are untouched.
    for (int k = 0; k < 3; ++k) {
        if (k == i)
            continue;
        EXPECT_EQ(frame.atoms.x[k], before.atoms.x[k]);
        EXPECT_EQ(frame.atoms.y[k], before.atoms.y[k]);
        EXPECT_EQ(frame.atoms.z[k], before.atoms.z[k]);
    }
}

TEST(TranslateMoveUnit, BootstrapPassIsNoOp) {
    cpu::Frame frame(2);
    frame.atoms.x[0] = 0.0;
    frame.atoms.x[1] = 4.0;
    cpu::Frame before = frame;

    std::mt19937_64 rng(1u);
    cpu::mc::Transaction txn;
    cpu::Context ctx;
    ctx.rng = &rng;
    ctx.proposal = &txn;

    cpu::mc::TranslateMove mv{0.5};
    cpu::LaunchContext lctx{ctx, frame, /*stage=*/0, /*step=*/-1}; // bootstrap
    mv.launch(lctx);

    EXPECT_EQ(txn.size(), 0u); // proposed nothing
    EXPECT_EQ(frame.atoms.x[0], before.atoms.x[0]);
    EXPECT_EQ(frame.atoms.x[1], before.atoms.x[1]);
}

TEST(TranslateMoveUnit, RejectsNonPositiveDelta) {
    EXPECT_THROW((cpu::mc::TranslateMove{0.0}), std::invalid_argument);
    EXPECT_THROW((cpu::mc::TranslateMove{-0.1}), std::invalid_argument);
}

// ── MC driver — build-time validation + config guards (no run) ────────────

TEST(MCBuildUnit, RejectsNonPositiveTemperature) {
    EXPECT_THROW((cpu::mc::MC{0.0, 1u}), std::invalid_argument);
    EXPECT_THROW((cpu::mc::MC{-300.0, 1u}), std::invalid_argument);
}

TEST(MCBuildUnit, ThrowsWithoutMover) {
    cpu::mc::MC mc{300.0, 1u};
    mc.set_pair(cpu::PairLJCut{kEps, kSigma, kCut});
    mc.set_neighbor(cpu::md::NeighborList<>{kCut, kSkin, 1, 0, true});
    EXPECT_THROW(mc.build(), std::runtime_error);
}

TEST(MCBuildUnit, ThrowsWithoutEnergyEvaluator) {
    cpu::mc::MC mc{300.0, 1u};
    mc.set_mover(cpu::mc::TranslateMove{0.2});
    EXPECT_THROW(mc.build(), std::runtime_error);
}

TEST(MCBuildUnit, ThrowsWithPairButNoNeighborList) {
    cpu::mc::MC mc{300.0, 1u};
    mc.set_pair(cpu::PairLJCut{kEps, kSigma, kCut});
    mc.set_mover(cpu::mc::TranslateMove{0.2});
    EXPECT_THROW(mc.build(), std::runtime_error);
}

TEST(MCBuildUnit, SucceedsWhenFullyAssembled) {
    cpu::mc::MC mc{300.0, 1u};
    mc.set_pair(cpu::PairLJCut{kEps, kSigma, kCut});
    mc.set_neighbor(cpu::md::NeighborList<>{kCut, kSkin, 1, 0, true});
    mc.set_mover(cpu::mc::TranslateMove{0.2});
    EXPECT_NO_THROW(mc.build());
}

TEST(MCBuildUnit, SetMoverIsIdempotentOnSameClass) {
    cpu::mc::MC mc{300.0, 1u};
    mc.set_mover(cpu::mc::TranslateMove{0.2});
    // Singular-by-name: re-setting the same Mover class is a no-op (v1 ships
    // only TranslateMove, so the different-class throw path is not exercisable).
    EXPECT_NO_THROW(mc.set_mover(cpu::mc::TranslateMove{0.5}));
}
