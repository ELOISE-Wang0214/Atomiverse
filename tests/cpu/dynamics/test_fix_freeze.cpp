/// @file tests/cpu/dynamics/test_fix_freeze.cpp
/// @brief Unit tests for cpu::md::FixFreeze (type-selected immobility).

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/dynamics/fix_freeze.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

cpu::Frame make_toy() {
    // 2 wall (type 2) + 2 fluid (type 1)
    cpu::Frame f(4);
    f.span(kAtom_type)[0] = 2;
    f.span(kAtom_type)[1] = 2;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_type)[3] = 1;
    for (int i = 0; i < 4; ++i) {
        f.span(kAtom_x)[i] = double(i);
        f.span(kAtom_y)[i] = 0.0;
        f.span(kAtom_z)[i] = 0.0;
        f.span(kAtom_vx)[i] = 1.0;
        f.span(kAtom_vy)[i] = 1.0;
        f.span(kAtom_vz)[i] = 1.0;
        f.span(kAtom_fx)[i] = 3.0;
        f.span(kAtom_fy)[i] = 3.0;
        f.span(kAtom_fz)[i] = 3.0;
        f.span(kAtom_mass)[i] = 12.0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.simbox = cpu::Box::cubic(10.0);
    return f;
}

} // namespace

TEST(FixFreeze, ZeroForceOnFrozenType) {
    auto frame = make_toy();
    cpu::Context runtime{};
    cpu::md::FixFreeze freeze(/*atom_type=*/2);

    cpu::BuildContext bctx{runtime, frame};
    freeze.build(bctx);
    EXPECT_EQ(freeze.n_constraints(runtime, frame), 3 * 2);

    cpu::LaunchContext lctx{runtime, frame, ATV::slot_of(md::Phase::post_force),
                            /*step=*/0};
    freeze.launch(lctx);

    EXPECT_DOUBLE_EQ(frame.span(kAtom_fx)[0], 0.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fy)[1], 0.0);
    // Fluid forces untouched.
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fx)[2], 3.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_fz)[3], 3.0);
}

TEST(FixFreeze, RestorePositionAndZeroVelocity) {
    auto frame = make_toy();
    cpu::Context runtime{};
    cpu::md::FixFreeze freeze(2);
    cpu::BuildContext bctx{runtime, frame};
    freeze.build(bctx);

    // Move walls and keep nonzero velocity.
    frame.span(kAtom_x)[0] = 99.0;
    frame.span(kAtom_y)[1] = -7.0;
    frame.span(kAtom_vx)[0] = 5.0;
    frame.span(kAtom_vz)[1] = -2.0;

    cpu::LaunchContext lctx{runtime, frame,
                            ATV::slot_of(md::Phase::initial_integrate), 1};
    freeze.launch(lctx);

    EXPECT_DOUBLE_EQ(frame.span(kAtom_x)[0], 0.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_y)[1], 0.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_vx)[0], 0.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_vz)[1], 0.0);
    // Fluid free.
    EXPECT_DOUBLE_EQ(frame.span(kAtom_x)[2], 2.0);
    EXPECT_DOUBLE_EQ(frame.span(kAtom_vx)[2], 1.0);
}
