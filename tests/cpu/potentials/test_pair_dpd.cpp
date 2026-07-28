/// @file tests/cpu/potentials/test_pair_dpd.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/potentials/pair_dpd.h"

using namespace ATV;

TEST(PairDpd, SoftConservativeEnergy) {
    // rc deliberately != 1: the Groot-Warren conservative energy is
    // U = ½ a r_c (1 − r/r_c)², and the shipped kernel dropped the r_c factor
    // for a long time without any gate noticing, because every DPD test and
    // regression in the tree used r_c = 1 Å, where the bug is invisible.
    // a is a force amplitude [eV/Å], so U must carry one power of length.
    const double a = 25.0, gamma = 0.0, T = 1.0, rc = 2.5, dt = 0.01, r = 1.25;
    // γ=0 → no D/R; w = 1 − r/rc = 0.5 → U = 0.5*25*2.5*0.25 = 7.8125 eV
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
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;

    cpu::Context ctx;
    cpu::PairDpd dpd(a, gamma, T, rc, dt, 1);
    dpd.set_step(0);
    dpd.build({ctx, f});
    dpd.launch({ctx, f, 0, 0});
    const double w = 1.0 - r / rc;
    EXPECT_NEAR(f.meta.get(kEnergy), 0.5 * a * rc * w * w, 1e-12);
    // Repulsive: F0x < 0
    EXPECT_LT(f.span(kAtom_fx)[0], 0.0);
    EXPECT_NEAR(f.span(kAtom_fx)[0], -f.span(kAtom_fx)[1], 1e-12);
    // The force the same loop applies must be −dU/dr: F = a w. Pinning both
    // together is what makes the missing r_c impossible to reintroduce —
    // energy and force disagreed by exactly r_c while this test passed.
    EXPECT_NEAR(-f.span(kAtom_fx)[0], a * w, 1e-12);
    dpd.destroy({ctx, f});
}
