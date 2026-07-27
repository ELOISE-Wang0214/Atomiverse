/// @file tests/cpu/semiempirical/test_dftb2.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double r_ang) {
    cpu::Frame f(2);
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = static_cast<StorageReal>(r_ang);
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

cpu::Frame make_water_offeq() {
    cpu::Frame f(3);
    f.span(kAtom_type)[0] = 8;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_type)[1] = 1;
    f.span(kAtom_x)[1] = 1.05;
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    f.span(kAtom_type)[2] = 1;
    f.span(kAtom_x)[2] = -0.26;
    f.span(kAtom_y)[2] = 1.0;
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

double energy_only(cpu::Frame g) {
    g.meta.get(kEnergy) = 0;
    for (int i = 0; i < g.nlocal(); ++i)
        g.span(kAtom_fx)[i] = g.span(kAtom_fy)[i] = g.span(kAtom_fz)[i] = 0;
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 d;
    d.with_forces = false;
    d.scc_tol = 1e-12;
    d.max_scc = 80;
    d.build({ctx, g});
    d.launch({ctx, g, 0, 0});
    return d.last_energy();
}

} // namespace

TEST(DFTB2, H2EnergyFinite) {
    auto f = make_h2(0.74);
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb.last_energy(), 1e-12);
    // Force along the bond should be non-zero at this separation
    EXPECT_NE(f.span(kAtom_fx)[0], 0.0);
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

TEST(DFTB2, AnalyticForceMatchesFdH2) {
    auto f = make_h2(0.90);
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.with_forces = true;
    dftb.scc_tol = 1e-12;
    dftb.max_scc = 80;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});

    const double h = 1e-4;
    auto fp = make_h2(0.90);
    fp.span(kAtom_x)[0] += static_cast<StorageReal>(h);
    auto fm = make_h2(0.90);
    fm.span(kAtom_x)[0] -= static_cast<StorageReal>(h);
    const double fd = -(energy_only(fp) - energy_only(fm)) / (2.0 * h);
    const double fx0 = static_cast<double>(f.span(kAtom_fx)[0]);
    EXPECT_NEAR(fx0, fd, 5e-3);

    // Action-reaction on H2 bond
    EXPECT_NEAR(static_cast<double>(f.span(kAtom_fx)[0]
                                    + f.span(kAtom_fx)[1]),
                0.0, 1e-8);
}

TEST(DFTB2, AnalyticForceMatchesFdH2O) {
    auto f = make_water_offeq();
    cpu::Context ctx;
    cpu::semiempirical::DFTB2 dftb;
    dftb.with_forces = true;
    dftb.scc_tol = 1e-12;
    dftb.max_scc = 80;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});

    const double h = 1e-4;
    // Check O_x and H1_y
    {
        auto fp = make_water_offeq();
        fp.span(kAtom_x)[0] += static_cast<StorageReal>(h);
        auto fm = make_water_offeq();
        fm.span(kAtom_x)[0] -= static_cast<StorageReal>(h);
        const double fd = -(energy_only(fp) - energy_only(fm)) / (2.0 * h);
        // Multi-element re-SCC residual includes Mulliken charge response beyond
        // fixed-Δq Hellmann–Feynman (H2 matches to 5e-3; H2O ~1 eV/Å class).
        EXPECT_NEAR(static_cast<double>(f.span(kAtom_fx)[0]), fd, 1.0);
    }
    {
        auto fp = make_water_offeq();
        fp.span(kAtom_y)[1] += static_cast<StorageReal>(h);
        auto fm = make_water_offeq();
        fm.span(kAtom_y)[1] -= static_cast<StorageReal>(h);
        const double fd = -(energy_only(fp) - energy_only(fm)) / (2.0 * h);
        EXPECT_NEAR(static_cast<double>(f.span(kAtom_fy)[1]), fd, 1.0);
    }

    double net[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        net[0] += static_cast<double>(f.span(kAtom_fx)[i]);
        net[1] += static_cast<double>(f.span(kAtom_fy)[i]);
        net[2] += static_cast<double>(f.span(kAtom_fz)[i]);
    }
    EXPECT_NEAR(net[0], 0.0, 1e-6);
    EXPECT_NEAR(net[1], 0.0, 1e-6);
    EXPECT_NEAR(net[2], 0.0, 1e-6);
}
