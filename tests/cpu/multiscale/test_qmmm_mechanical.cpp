#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_mechanical.h"
#include <gtest/gtest.h>
using namespace ATV;
TEST(QmmmMechanical, H2PlusAr) {
    cpu::Frame f(3);
    // H2 QM + Ar MM
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = 5.0;
    f.span(kAtom_y)[2] = 0;
    f.span(kAtom_z)[2] = 0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 18;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.span(kAtom_tag)[2] = 3;
    f.span(kAtom_mass)[0] = 1;
    f.span(kAtom_mass)[1] = 1;
    f.span(kAtom_mass)[2] = 40;
    f.box = cpu::Box::cubic(30);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 3; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    cpu::Context ctx;
    cpu::QmmmMechanical qmmm({1, 2}, 0.01, 3.4, 10.0);
    qmmm.build({ctx, f});
    qmmm.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(qmmm.last_qm_energy()));
    EXPECT_TRUE(std::isfinite(qmmm.last_mm_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy),
                qmmm.last_qm_energy() + qmmm.last_mm_energy(), 1e-10);
}
