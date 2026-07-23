/// @file tests/cpu/scf/test_ri_mp2.cpp
/// @brief RI-MP2 vs canonical MP2 on H2/STO-3G.
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

TEST(RiMp2, H2Sto3gCloseToCanonical) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[0] = f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[0] = f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(20);

    cpu::scf::RHF can;
    can.basis = cpu::scf::BasisSet::STO3G;
    can.tol = 1e-10;
    auto rc = can.mp2(f);
    EXPECT_FALSE(rc.ri);
    EXPECT_LT(rc.e_corr, 0.0);

    cpu::scf::RHF ri;
    ri.basis = cpu::scf::BasisSet::STO3G;
    ri.tol = 1e-10;
    ri.ri = cpu::scf::RI::JK;
    ri.aux_basis = cpu::scf::AuxBasis::Def2JKFIT;
    auto rr = ri.mp2(f);
    EXPECT_TRUE(rr.ri);
    EXPECT_LT(rr.e_corr, 0.0);
    EXPECT_NEAR(rr.e_corr, rc.e_corr, 1e-4); // Ha; DF error small on H2
}
