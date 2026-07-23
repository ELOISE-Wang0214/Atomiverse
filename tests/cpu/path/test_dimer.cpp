/// @file tests/cpu/path/test_dimer.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/path/dimer.h"

using namespace ATV;

TEST(Dimer, InvertedHarmonicFindsOrigin) {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.3;
    f.span(kAtom_y)[0] = 0.4;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_mass)[0] = 1.0;
    f.span(kAtom_type)[0] = 0;
    f.span(kAtom_tag)[0] = 1;
    f.simbox = cpu::Box::cubic(20.0);

    
    cpu::path::Dimer dimer(0.01, 0.02, 10, 1e-3, 0.05, 7);
    dimer.use_inverted_harmonic(1.0, 1.0);
    const int steps = dimer.run(f, 200);

    EXPECT_GT(steps, 0);
    EXPECT_TRUE(dimer.converged());
    const double r =
        std::hypot(f.span(kAtom_x)[0], f.span(kAtom_y)[0], f.span(kAtom_z)[0]);
    EXPECT_LT(r, 0.15);
}

TEST(Dimer, RejectsBadDelta) {
    EXPECT_THROW((cpu::path::Dimer{0.0, 0.05, 8, 1e-3, 0.05, 1}),
                 std::invalid_argument);
}
