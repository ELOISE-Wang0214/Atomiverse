/// @file tests/cpu/test_analysis_executor.cpp
/// @brief RED tests for `cpu::AnalysisExecutor` lifecycle and backpressure policy.
///
/// Spec: .claude/specs/parallel-md-engine.md § B.2.
/// Verifies:
///   - clean init/shutdown (no thread leak)
///   - round-robin slot IDs + drop_oldest counter
///   - drop_oldest / block / grow backpressure under load
///   - shutdown drains and joins cleanly

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cstdint>
#include <thread>
#include <vector>

#include "atomiverse/core/types.h"
#include "atomiverse/cpu/analysis_executor.h"
#include "atomiverse/cpu/async_component_base.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame_snapshot.h"

using namespace ATV;

namespace {

/// @brief Test-only async component used to measure submission behaviour.
class SlowCapture : public cpu::AsyncComponentBase {
  public:
    explicit SlowCapture(
        std::chrono::milliseconds delay = std::chrono::milliseconds(50))
        : delay_(delay) {}

    std::string name() const override { return "slow_capture"; }
    cpu::CaptureMask capture_mask() const override {
        cpu::CaptureMask m;
        m.positions = true;
        return m;
    }

    void on_snapshot(const cpu::FrameSnapshot& /*snap*/) override {
        std::this_thread::sleep_for(delay_);
        ++n_handled_;
    }

    int handled() const { return n_handled_.load(std::memory_order_relaxed); }

  private:
    std::chrono::milliseconds delay_;
    std::atomic<int> n_handled_{0};
};

} // namespace

// ── clean init/shutdown ──────────────────────────────────────────────

TEST(AnalysisExecutor, InitShutdownClean) {
    cpu::AnalysisExecutor exec;
    cpu::Context ctx;

    cpu::AnalysisExecutor::Config cfg;
    cfg.n_threads = 2;
    cfg.ring_size = 3;
    cfg.policy = ATV::BackpressurePolicy::drop_oldest;

    exec.init(cfg, ctx);
    // No implicit leak: shutdown must join all workers.
    exec.shutdown();
    SUCCEED();
}

// ── round-robin slot IDs + drop counter ──────────────────────────────

TEST(AnalysisExecutor, AcquireSlotRoundRobin) {
    cpu::AnalysisExecutor exec;
    cpu::Context ctx;
    cpu::AnalysisExecutor::Config cfg;
    cfg.n_threads = 0; // no workers — just test the ring mechanics
    cfg.ring_size = 3;
    cfg.policy = ATV::BackpressurePolicy::drop_oldest;
    exec.init(cfg, ctx);

    int s0 = exec.acquire_slot();
    int s1 = exec.acquire_slot();
    int s2 = exec.acquire_slot();
    EXPECT_EQ(s0, 0);
    EXPECT_EQ(s1, 1);
    EXPECT_EQ(s2, 2);

    // 4th acquire with drop_oldest wraps and bumps the drop counter.
    int s3 = exec.acquire_slot();
    EXPECT_EQ(s3, 0);
    EXPECT_EQ(exec.n_dropped(), 1);

    exec.shutdown();
}

// ── drop_oldest backpressure ─────────────────────────────────────────

TEST(AnalysisExecutor, DropOldestBackpressure) {
    cpu::AnalysisExecutor exec;
    cpu::Context ctx;
    cpu::AnalysisExecutor::Config cfg;
    cfg.n_threads = 1;
    cfg.ring_size = 3;
    cfg.policy = ATV::BackpressurePolicy::drop_oldest;
    exec.init(cfg, ctx);

    SlowCapture sub{std::chrono::milliseconds(100)};
    for (int k = 0; k < 5; ++k) {
        int slot = exec.acquire_slot();
        exec.publish(slot, k);
        exec.submit(&sub, slot);
    }

    // With ring_size=3 and 5 submissions, 2 overwrites must be counted.
    EXPECT_EQ(exec.n_dropped(), 2);

    exec.shutdown();
}

// ── block backpressure ───────────────────────────────────────────────

TEST(AnalysisExecutor, BlockBackpressure) {
    cpu::AnalysisExecutor exec;
    cpu::Context ctx;
    cpu::AnalysisExecutor::Config cfg;
    cfg.n_threads = 1;
    cfg.ring_size = 3;
    cfg.policy = ATV::BackpressurePolicy::block;
    exec.init(cfg, ctx);

    SlowCapture sub{std::chrono::milliseconds(80)};

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < 5; ++k) {
        int slot = exec.acquire_slot();
        exec.publish(slot, k);
        exec.submit(&sub, slot);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    auto elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    // With block policy, the 4th/5th acquire must wait for worker to free a slot.
    EXPECT_GE(elapsed_ms, 80)
        << "block policy did not block hot-path submit (elapsed " << elapsed_ms
        << " ms)";
    EXPECT_EQ(exec.n_dropped(), 0);

    exec.shutdown();
}

// ── grow backpressure ────────────────────────────────────────────────

TEST(AnalysisExecutor, GrowBackpressure) {
    cpu::AnalysisExecutor exec;
    cpu::Context ctx;
    cpu::AnalysisExecutor::Config cfg;
    cfg.n_threads = 1;
    cfg.ring_size = 3;
    cfg.policy = ATV::BackpressurePolicy::grow;
    exec.init(cfg, ctx);

    SlowCapture sub{std::chrono::milliseconds(40)};
    for (int k = 0; k < 6; ++k) {
        int slot = exec.acquire_slot();
        exec.publish(slot, k);
        exec.submit(&sub, slot);
    }

    EXPECT_EQ(exec.n_dropped(), 0);
    EXPECT_TRUE(exec.ring_size_grew())
        << "grow policy did not expand ring capacity under load";

    exec.shutdown();
}

// ── shutdown drains in-flight jobs ───────────────────────────────────

TEST(AnalysisExecutor, ShutdownWithInFlightJobs) {
    cpu::AnalysisExecutor exec;
    cpu::Context ctx;
    cpu::AnalysisExecutor::Config cfg;
    cfg.n_threads = 2;
    cfg.ring_size = 3;
    cfg.policy = ATV::BackpressurePolicy::block;
    exec.init(cfg, ctx);

    SlowCapture sub{std::chrono::milliseconds(30)};
    for (int k = 0; k < 3; ++k) {
        int slot = exec.acquire_slot();
        exec.publish(slot, k);
        exec.submit(&sub, slot);
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    exec.shutdown();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto elapsed_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    // Must join within a generous deadline; no hang.
    EXPECT_LT(elapsed_ms, 2000);
}
