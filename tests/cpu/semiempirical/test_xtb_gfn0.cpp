#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/xtb_gfn0.h"
#include <gtest/gtest.h>
using namespace ATV;
TEST(XtbGfn0, H2OEnergyFinite) {
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.96;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = -0.24;
    f.span(kAtom_y)[2] = 0.93;
    f.span(kAtom_z)[2] = 0;
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.simbox = cpu::Box::cubic(20);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::XtbGfn0 xtb;
    xtb.build({ctx, f});
    xtb.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(xtb.last_energy()));
}
