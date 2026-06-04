/// @file tests/cpu/md/test_frame_tag_image.cpp
/// @brief RED tests for Phase 1 `cpu::Frame::Atoms` additions: `tag` and `image`.
///
/// Spec: .claude/specs/parallel-md-engine.md § A.2 (Frame additions).
/// These tests fail to compile until the production code adds `tag` (xtensor<int64_t,1>)
/// and `image` (xtensor<int8_t,2>, shape [n,3]) to `cpu::Frame::Atoms`, with the
/// default-construction invariant `tag[i]=i` and `image=0`.

#include <gtest/gtest.h>

#include <cstdint>
#include <type_traits>

#include "atomiverse/cpu/frame.h"

using namespace ATV;

// ── identity tag defaults ────────────────────────────────────────────

TEST(FrameTagImage, DefaultConstructSetsIdentityTags) {
    constexpr int n = 100;
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(frame.atoms.tag(i), static_cast<std::int64_t>(i))
            << "tag at index " << i << " should equal i by default";
    }
}

// ── image defaults to zero ───────────────────────────────────────────

TEST(FrameTagImage, DefaultConstructZerosImage) {
    constexpr int n = 100;
    cpu::Frame frame(n);
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
            EXPECT_EQ(frame.atoms.image(i, d), static_cast<std::int8_t>(0))
                << "image(" << i << "," << d << ") should default to 0";
        }
    }
}

// ── tag value_type is int64_t ────────────────────────────────────────

TEST(FrameTagImage, TagDtypeIsInt64) {
    cpu::Frame frame(4);
    using tag_t = typename std::decay_t<decltype(frame.atoms.tag)>::value_type;
    static_assert(std::is_same_v<tag_t, std::int64_t>,
                  "Frame::Atoms::tag must be xtensor<int64_t, 1>");
    SUCCEED();
}

// ── image shape [n, 3] ───────────────────────────────────────────────

TEST(FrameTagImage, ImageShapeIsNx3) {
    constexpr int n = 37;
    cpu::Frame frame(n);
    // Per-atom SoA fields are sized to capacity `nmax` (LAMMPS DELTA growth,
    // floored at 64), not exactly `n`. The invariant is that `image` is a
    // `[nmax, 3]` array with the same row count as the other per-atom fields.
    EXPECT_GE(static_cast<int>(frame.atoms.image.shape()[0]), n);
    EXPECT_EQ(frame.atoms.image.shape()[0], frame.atoms.x.shape()[0]);
    EXPECT_EQ(static_cast<int>(frame.atoms.image.shape()[1]), 3);
}
