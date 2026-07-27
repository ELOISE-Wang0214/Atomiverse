/// @file tests/cpu/semiempirical/test_dftb0.cpp
/// Unit tests for non-SCC DFTB0 (H2 + table gates).
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb0.h"

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

} // namespace

TEST(DFTB0, H2EnergyFinite) {
    auto f = make_h2();
    cpu::Context ctx;
    cpu::semiempirical::DFTB0 dftb0;
    dftb0.build({ctx, f});
    dftb0.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb0.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb0.last_energy(), 1e-12);
    EXPECT_NE(f.span(kAtom_fx)[0], 0.0);
}

TEST(DFTB0, RejectsUnsupportedElement) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26; // Fe
    f.box = cpu::Box::cubic(10.0);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB0 dftb0;
    dftb0.build({ctx, f});
    EXPECT_THROW(dftb0.launch({ctx, f, 0, 0}), std::runtime_error);
}

TEST(DFTB0, ForceMatchesFiniteDifference) {
    // Slightly stretch H2 so F is nonzero; FD of same non-SCC energy path.
    auto f = make_h2(0.90);
    cpu::Context ctx;
    cpu::semiempirical::DFTB0 dftb0;
    dftb0.with_forces = true;
    dftb0.build({ctx, f});
    dftb0.launch({ctx, f, 0, 0});
    const double fx0 = static_cast<double>(f.span(kAtom_fx)[0]);

    const double h = 1e-4; // Å
    auto e_at = [&](double dx) {
        auto g = make_h2(0.90);
        g.span(kAtom_x)[0] += static_cast<StorageReal>(dx);
        g.meta.get(kEnergy) = 0;
        for (int i = 0; i < 2; ++i)
            g.span(kAtom_fx)[i] = g.span(kAtom_fy)[i] = g.span(kAtom_fz)[i] = 0;
        cpu::semiempirical::DFTB0 d;
        d.with_forces = false;
        d.build({ctx, g});
        d.launch({ctx, g, 0, 0});
        return d.last_energy();
    };
    const double fd = -(e_at(h) - e_at(-h)) / (2.0 * h);
    EXPECT_NEAR(fx0, fd, 0.05) << "force=" << fx0 << " fd=" << fd;
}

TEST(DFTB0, EnergyDependsOnBondLength) {
    cpu::Context ctx;
    auto e_at = [&](double r) {
        auto f = make_h2(r);
        cpu::semiempirical::DFTB0 d;
        d.with_forces = false;
        d.build({ctx, f});
        d.launch({ctx, f, 0, 0});
        return d.last_energy();
    };
    EXPECT_GT(std::fabs(e_at(0.74) - e_at(3.0)), 1e-4);
}
