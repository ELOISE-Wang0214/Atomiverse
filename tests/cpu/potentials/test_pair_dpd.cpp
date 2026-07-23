/// @file tests/cpu/potentials/test_pair_dpd.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/potentials/pair_dpd.h"

using namespace ATV;

TEST(PairDpd, SoftConservativeEnergy) {
    const double a = 25.0, gamma = 0.0, T = 1.0, rc = 1.0, dt = 0.01, r = 0.5;
    // γ=0 → no D/R; U = ½ a w², w = 1 − r/rc = 0.5 → U = 0.5*25*0.25 = 3.125
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = 1.0;
        f.span(kAtom_type)[i] = 0;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_vx)[i] = f.span(kAtom_vy)[i] = f.span(kAtom_vz)[i] = 0;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.simbox = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;

    cpu::Context ctx;
    cpu::PairDpd dpd(a, gamma, T, rc, dt, 1);
    dpd.set_step(0);
    dpd.build({ctx, f});
    dpd.launch({ctx, f, 0, 0});
    const double w = 1.0 - r / rc;
    EXPECT_NEAR(f.meta.get(kEnergy), 0.5 * a * w * w, 1e-12);
    // Repulsive: F0x < 0
    EXPECT_LT(f.span(kAtom_fx)[0], 0.0);
    EXPECT_NEAR(f.span(kAtom_fx)[0], -f.span(kAtom_fx)[1], 1e-12);
    dpd.destroy({ctx, f});
}
