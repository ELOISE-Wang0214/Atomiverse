#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/multiscale/qmmm_electrostatic.h"
#include <gtest/gtest.h>
using namespace ATV;
TEST(QmmmElectrostatic, ChargeShiftsEnergy) {
    auto make = [](double q) {
        cpu::Frame f(3);
        f.span(kAtom_x)[0] = 0;
        f.span(kAtom_y)[0] = 0;
        f.span(kAtom_z)[0] = 0;
        f.span(kAtom_x)[1] = 0.74;
        f.span(kAtom_y)[1] = 0;
        f.span(kAtom_z)[1] = 0;
        f.span(kAtom_x)[2] = 3.0;
        f.span(kAtom_y)[2] = 0;
        f.span(kAtom_z)[2] = 0;
        f.span(kAtom_type)[0] = 1;
        f.span(kAtom_type)[1] = 1;
        f.span(kAtom_type)[2] = 1;
        f.span(kAtom_tag)[0] = 1;
        f.span(kAtom_tag)[1] = 2;
        f.span(kAtom_tag)[2] = 3;
        f.box = cpu::Box::cubic(20);
        f.meta.get(kEnergy) = 0;
        cpu::Context ctx;
        cpu::QmmmElectrostatic qm({1, 2}, q, 0.01, 3.4, 10.0);
        qm.build({ctx, f});
        qm.launch({ctx, f, 0, 0});
        return f.meta.get(kEnergy);
    };
    const double e0 = make(0.0);
    const double e1 = make(0.5);
    EXPECT_TRUE(std::isfinite(e0) && std::isfinite(e1));
    EXPECT_GT(std::fabs(e1 - e0), 1e-6);
}
