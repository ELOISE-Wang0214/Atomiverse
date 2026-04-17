/// @file tests/cpu/test_snapshot_publisher.cpp
/// @brief RED tests for `cpu::SnapshotPublisher` zero-subscriber fast-path
///        and field-copy correctness.
///
/// Spec: .claude/specs/parallel-md-engine.md § B.3.
/// Required:
///   - zero subscribers → mean launch() cost < 100 ns (10 000 iters)
///   - positions, box, step propagate into FrameSnapshot
///   - CaptureMask respected (excluded fields remain empty/unmodified)

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <mutex>
#include <string>

#include "atomiverse/cpu/async_component_base.h"
#include "atomiverse/cpu/backend.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/frame_snapshot.h"
#include "atomiverse/cpu/md/async_observer.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/core/mat3.h"

using namespace ATV;

namespace {

cpu::Frame make_small_frame(int n) {
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = 1.0 + 0.1 * i;
        frame.atoms.y[i] = 2.0 + 0.1 * i;
        frame.atoms.z[i] = 3.0 + 0.1 * i;
        frame.atoms.vx[i] = 0.01 * i;
        frame.atoms.vy[i] = 0.02 * i;
        frame.atoms.vz[i] = 0.03 * i;
        frame.atoms.mass[i] = 1.0;
    }
    frame.box = Mat3::cubic(10.0);
    return frame;
}

/// @brief Capture target — records the last snapshot received.
///
/// Inherits `cpu::md::AsyncObserver` to get the MD-flavoured default
/// binding (`md::Phase::end_of_step @ +20`), so `sim.add(cap)` works
/// without explicit `on()` wiring.
class PositionCapture : public cpu::md::AsyncObserver {
public:
    std::string name() const override { return "position_capture"; }
    cpu::CaptureMask capture_mask() const override {
        cpu::CaptureMask m;
        m.positions = true;
        m.velocities = false;
        return m;
    }
    void on_snapshot(const cpu::FrameSnapshot& snap) override {
        std::lock_guard<std::mutex> g(mtx_);
        last_ = snap;        // copy
        ++n_handled_;
    }
    int handled() const { return n_handled_.load(std::memory_order_relaxed); }
    cpu::FrameSnapshot copy_last() const {
        std::lock_guard<std::mutex> g(mtx_);
        return last_;
    }
private:
    mutable std::mutex mtx_;
    cpu::FrameSnapshot last_;
    std::atomic<int> n_handled_{0};
};

} // namespace

// ── zero-subscriber fast path ────────────────────────────────────────

TEST(SnapshotPublisher, ZeroSubscribersEarlyReturn) {
    auto frame = make_small_frame(100);
    cpu::Context ctx;
    cpu::SnapshotPublisher pub;
    cpu::BuildContext bctx{ctx, frame};
    pub.build(bctx);

    constexpr int iters = 10'000;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (int k = 0; k < iters; ++k) {
        cpu::LaunchContext lctx{ctx, frame, 0, k};
        pub.launch(lctx);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    auto mean_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / iters;
    EXPECT_LT(mean_ns, 100)
        << "zero-subscriber mean launch() cost " << mean_ns
        << " ns exceeds 100 ns budget";
}

// ── positions propagate end-to-end ───────────────────────────────────

TEST(SnapshotPublisher, SnapshotCopiesPositions) {
    auto frame = make_small_frame(64);

    // `cap` MUST outlive `sim`: sim only stores a non-owning pointer (lvalue
    // overload of Driver::add). Inspecting cap after the run requires it
    // to still be alive in the same scope.
    PositionCapture cap;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{1e-6});
    sim.add(cap);
    sim.build();
    sim.run(frame, 5);  // give the publisher a step where async fires (step >= 1)

    const int n = frame.atoms.size();
    ASSERT_GE(cap.handled(), 1) << "no snapshot was delivered";
    auto snap = cap.copy_last();
    ASSERT_EQ(static_cast<int>(snap.x.size()), n);
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(snap.x[i], frame.atoms.x[i], 1e-12);
        EXPECT_NEAR(snap.y[i], frame.atoms.y[i], 1e-12);
        EXPECT_NEAR(snap.z[i], frame.atoms.z[i], 1e-12);
    }
}

// ── box + step propagate ─────────────────────────────────────────────

TEST(SnapshotPublisher, SnapshotCopiesBoxAndStep) {
    auto frame = make_small_frame(16);

    PositionCapture cap;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{1e-6});
    sim.add(cap);
    sim.build();
    sim.run(frame, 5);

    ASSERT_GE(cap.handled(), 1);
    auto snap = cap.copy_last();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(snap.box(i, j), frame.box(i, j), 1e-12);
    EXPECT_GE(snap.step, 1) << "publisher fires after step 0; first capture must be step >= 1";
}

// ── CaptureMask respected ────────────────────────────────────────────

TEST(SnapshotPublisher, CaptureMaskRespected) {
    auto frame = make_small_frame(32);

    PositionCapture cap;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{1e-6});
    sim.add(cap);
    sim.build();
    sim.run(frame, 5);

    ASSERT_GE(cap.handled(), 1);
    auto snap = cap.copy_last();
    EXPECT_EQ(static_cast<int>(snap.x.size()), frame.atoms.size());
    EXPECT_EQ(static_cast<std::size_t>(snap.vx.size()), static_cast<std::size_t>(0))
        << "velocities must remain empty when CaptureMask.velocities=false";
    EXPECT_EQ(static_cast<std::size_t>(snap.vy.size()), static_cast<std::size_t>(0));
    EXPECT_EQ(static_cast<std::size_t>(snap.vz.size()), static_cast<std::size_t>(0));
}
