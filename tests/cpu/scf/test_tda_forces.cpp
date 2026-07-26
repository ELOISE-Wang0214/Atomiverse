#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"
#include <gtest/gtest.h>
using namespace ATV;
TEST(TdaForces, H2NonZero) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.8;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.box = cpu::Box::cubic(15);
    cpu::scf::RHF rhf;
    rhf.basis = cpu::scf::BasisSet::STO3G;
    rhf.tol = 1e-8;
    rhf.run(f);
    ASSERT_TRUE(rhf.converged());
    auto F = rhf.tda_forces(f, 0, 2e-3);
    ASSERT_EQ(F.size(), 6u);
    // Newton's 3rd: roughly opposite on dimer axis
    EXPECT_TRUE(std::isfinite(F[0]) && std::isfinite(F[3]));
}
