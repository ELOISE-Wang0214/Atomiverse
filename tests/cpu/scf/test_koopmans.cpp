#ifdef ATV_HAS_RHF
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"
#include <gtest/gtest.h>
using namespace ATV;
static cpu::Frame water() {
    cpu::Frame f(3);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_x)[1] = 0.7572;
    f.span(kAtom_x)[2] = -0.7572;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_y)[2] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_z)[1] = -0.5860;
    f.span(kAtom_z)[2] = -0.5860;
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_mass)[0] = 15.999;
    f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_mass)[2] = 1.008;
    f.simbox = cpu::Box::cubic(20);
    return f;
}
TEST(Koopmans, HomoNegative) {
    auto f = water();
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    const auto& r = rhf.run(f);
    ASSERT_TRUE(r.converged);
    EXPECT_LT(r.eigenvalues_alpha(r.n_alpha - 1), 0.0);
}
#endif
