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

    cpu::path::DimerConfig cfg;
    cfg.delta = 0.01;
    cfg.dR = 0.02;
    cfg.max_rotations = 10;
    cfg.f_tol = 0.05;
    cfg.seed = 7;

    cpu::path::Dimer dimer(cfg);
    dimer.use_inverted_harmonic(1.0, 1.0);
    const int steps = dimer.run(f, 200);

    EXPECT_GT(steps, 0);
    EXPECT_TRUE(dimer.converged());
    const double r =
        std::hypot(f.span(kAtom_x)[0], f.span(kAtom_y)[0], f.span(kAtom_z)[0]);
    EXPECT_LT(r, 0.15);
}

TEST(Dimer, RejectsBadDelta) {
    cpu::path::DimerConfig cfg;
    cfg.delta = 0.0;
    EXPECT_THROW(cpu::path::Dimer{cfg}, std::invalid_argument);
}
