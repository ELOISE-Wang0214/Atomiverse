/// @file tests/cpu/scf/test_rks_tda.cpp
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/methods.h"

using namespace ATV;

TEST(RksTda, ThrowsBeforeRun) {
    cpu::scf::RKS rks(cpu::scf::XCType::LSDA);
    EXPECT_THROW(rks.tda(1), std::runtime_error);
}

TEST(RksTda, H2FirstSingletPositive) {
    // Minimal H2/STO-3G LSDA — first singlet TDA energy > 0 eV.
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.simbox = cpu::Box::cubic(20.0);

    cpu::scf::RKS rks(cpu::scf::XCType::LSDA);
    rks.basis = cpu::scf::BasisSet::STO3G;
    rks.tol = 1e-8;
    rks.max_iter = 80;
    rks.run(f);
    ASSERT_TRUE(rks.converged());
    auto es = rks.tda(1, cpu::scf::ExcitationChannel::singlet);
    ASSERT_FALSE(es.energy_ev.empty());
    EXPECT_GT(es.energy_ev[0], 0.5); // H2 first singlet is several eV
    EXPECT_TRUE(es.tda);
}
