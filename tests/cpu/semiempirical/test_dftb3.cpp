/// @file tests/cpu/semiempirical/test_dftb3.cpp
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"
#include "atomiverse/cpu/semiempirical/dftb3.h"

using namespace ATV;

namespace {

cpu::Frame make_water() {
    cpu::Frame f(3);
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_x)[1] = 0.96;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_x)[2] = -0.24;
    f.span(kAtom_y)[2] = 0.93;
    f.span(kAtom_z)[2] = 0.0;
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_mass)[i] = (i == 0) ? 15.999 : 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

} // namespace

TEST(DFTB3, H2OEnergyFinite) {
    auto f = make_water();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.with_forces = false;
    dftb3.scc_tol = 1e-10;
    dftb3.build({ctx, f});
    dftb3.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb3.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb3.last_energy(), 1e-12);
    EXPECT_GT(std::fabs(dftb3.last_energy()), 1e-3);
    EXPECT_TRUE(std::isfinite(dftb3.last_e2()));
    EXPECT_TRUE(std::isfinite(dftb3.last_e3()));
}

TEST(DFTB3, ThirdOrderDiffersFromDftb2) {
    auto f3 = make_water();
    auto f2 = make_water();
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 d3;
    d3.with_forces = false;
    d3.scc_tol = 1e-10;
    d3.build({ctx, f3});
    d3.launch({ctx, f3, 0, 0});
    cpu::semiempirical::DFTB2 d2;
    d2.with_forces = false;
    d2.scc_tol = 1e-10;
    d2.build({ctx, f2});
    d2.launch({ctx, f2, 0, 0});
    // Polar water: nonzero Mulliken Δq → E^(3) and full SCC shift differ.
    EXPECT_GT(std::fabs(d3.last_energy() - d2.last_energy()), 1e-6);
    EXPECT_GT(std::fabs(d3.last_e3()), 0.0);
}

TEST(DFTB3, RejectsUnsupportedElement) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26; // Fe
    f.box = cpu::Box::cubic(10.0);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB3 dftb3;
    dftb3.build({ctx, f});
    EXPECT_THROW(dftb3.launch({ctx, f, 0, 0}), std::runtime_error);
}
