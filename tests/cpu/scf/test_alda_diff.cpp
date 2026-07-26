#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"
#include <cmath>
#include <gtest/gtest.h>
using namespace ATV;
TEST(AldaDiff, LsdaDiffersFromPbePathFlag) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.box = cpu::Box::cubic(20);
    cpu::scf::RKS lsda(cpu::scf::XCType::LSDA);
    lsda.basis = cpu::scf::BasisSet::STO3G;
    lsda.tol = 1e-8;
    lsda.run(f);
    auto e_lsda = lsda.tda(1);
    EXPECT_TRUE(e_lsda.alda);
    cpu::scf::RKS pbe(cpu::scf::XCType::PBE);
    pbe.basis = cpu::scf::BasisSet::STO3G;
    pbe.tol = 1e-8;
    pbe.run(f);
    auto e_pbe = pbe.tda(1);
    EXPECT_FALSE(e_pbe.alda);
    // Full-grid ALDA moves the LSDA S1 relative to pure CIS@PBE.
    EXPECT_GT(std::fabs(e_lsda.energy_ev[0] - e_pbe.energy_ev[0]), 1e-3);
}
