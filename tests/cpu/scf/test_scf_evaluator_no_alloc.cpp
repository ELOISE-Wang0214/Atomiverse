#ifdef ATV_HAS_RHF
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/scf/methods.h"
#include <gtest/gtest.h>
using namespace ATV;
static cpu::Frame h2() {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = 1.008;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_vx)[0] = 0.001;
    f.span(kAtom_vx)[1] = -0.001;
    f.simbox = cpu::Box::cubic(20);
    return f;
}
TEST(SCFEvaluatorNoAlloc, MultiStepRuns) {
    auto frame = h2();
    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(md::Phase::force).add(cpu::scf::RHF{});
    sim.build();
    sim.run(frame, 3);
    EXPECT_LT(frame.meta.get(kEnergy), 0.0);
}
#endif
