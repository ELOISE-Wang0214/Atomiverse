/// @file tests/cpu/semiempirical/test_dftb2.cpp
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

TEST(DFTB2, H2EnergyFinite) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.simbox = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;

    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb.last_energy(), 1e-12);
    // Force along the bond should be non-zero at this separation
    EXPECT_NE(f.span(kAtom_fx)[0], 0.0);
}

TEST(DFTB2, RejectsUnsupportedElement) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26; // Fe — not in static H/C/N/O tables
    f.simbox = cpu::Box::cubic(10.0);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.build({ctx, f});
    EXPECT_THROW(dftb.launch({ctx, f, 0, 0}), std::runtime_error);
}
