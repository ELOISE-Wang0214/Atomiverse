/// @file tests/cpu/semiempirical/test_dftb2.cpp
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double r = 0.74) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.008;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    return f;
}

/// O₂ triplet geometry (experimental-ish bond length, Å).
cpu::Frame make_o2(double r = 1.21) {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 8;
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = r;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 15.999;
    f.span(kAtom_tag)[0] = 1;
    f.span(kAtom_tag)[1] = 2;
    f.box = cpu::Box::cubic(30.0);
    f.meta.get(kEnergy) = 0;
    for (int i = 0; i < 2; ++i)
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    return f;
}

} // namespace

TEST(DFTB2, H2EnergyFinite) {
    auto f = make_h2();
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb.last_energy(), 1e-12);
    // Force along the bond should be non-zero at this separation
    EXPECT_NE(f.span(kAtom_fx)[0], 0.0);
    // Closed shell: zero Mulliken spin
    ASSERT_EQ(dftb.last_mulliken_spin().size(), 2u);
    EXPECT_NEAR(dftb.last_mulliken_spin()[0], 0.0, 1e-10);
    EXPECT_NEAR(dftb.last_mulliken_spin()[1], 0.0, 1e-10);
}

TEST(DFTB2, RejectsUnsupportedElement) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26; // Fe — not in static H/C/N/O tables
    f.box = cpu::Box::cubic(10.0);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.build({ctx, f});
    EXPECT_THROW(dftb.launch({ctx, f, 0, 0}), std::runtime_error);
}

TEST(DFTB2, O2TripletSpinPopulation) {
    auto f = make_o2();
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.multiplicity = 3;
    dftb.with_forces = false;
    dftb.scc_tol = 1e-10;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
    ASSERT_EQ(dftb.last_mulliken_spin().size(), 2u);
    const double m0 = dftb.last_mulliken_spin()[0];
    const double m1 = dftb.last_mulliken_spin()[1];
    const double m_tot = m0 + m1;
    // Triplet: total Mulliken spin population ≈ 2 (two unpaired e⁻)
    EXPECT_NEAR(m_tot, 2.0, 0.25);
    EXPECT_GT(std::fabs(m0), 0.3);
    EXPECT_GT(std::fabs(m1), 0.3);
    // Symmetric O₂: equal spin on each O
    EXPECT_NEAR(m0, m1, 0.15);
}

TEST(DFTB2, RejectsBadMultiplicity) {
    auto f = make_o2();
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.multiplicity = 4; // 12 valence e → odd unpaired not consistent
    dftb.with_forces = false;
    dftb.build({ctx, f});
    EXPECT_THROW(dftb.launch({ctx, f, 0, 0}), std::runtime_error);
}
