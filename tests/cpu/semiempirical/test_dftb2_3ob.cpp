/// @file tests/cpu/semiempirical/test_dftb2_3ob.cpp
/// Unit tests for DFTB2 ThreeOb static SK tables (H/C/N/O).
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/basis_data/dftb_3ob_tables.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;
using namespace ATV::cpu::semiempirical;

namespace {

cpu::Frame make_ch4() {
    cpu::Frame f(5);
    f.span(kAtom_type)[0] = 6;
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    const double s = 1.09 / std::sqrt(3.0);
    const double hpos[4][3] = {
        {s, s, s}, {s, -s, -s}, {-s, s, -s}, {-s, -s, s},
    };
    for (int i = 0; i < 4; ++i) {
        f.span(kAtom_type)[i + 1] = 1;
        f.span(kAtom_x)[i + 1] = hpos[i][0];
        f.span(kAtom_y)[i + 1] = hpos[i][1];
        f.span(kAtom_z)[i + 1] = hpos[i][2];
    }
    for (int i = 0; i < 5; ++i) {
        f.span(kAtom_mass)[i] = (i == 0) ? 12.0 : 1.008;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(20.0);
    f.meta.get(kEnergy) = 0;
    return f;
}

} // namespace

TEST(DFTB23ob, TablesCoverHCNO) {
    using namespace basis_data;
    EXPECT_EQ(kDftb3obNElements, 4u);
    int mask = 0;
    for (std::size_t i = 0; i < kDftb3obNElements; ++i)
        mask |= (1 << kDftb3obElements[i].z);
    EXPECT_TRUE(mask & (1 << 1));
    EXPECT_TRUE(mask & (1 << 6));
    EXPECT_TRUE(mask & (1 << 7));
    EXPECT_TRUE(mask & (1 << 8));
    EXPECT_GE(kDftb3obNPairs, 16u);
}

TEST(DFTB23ob, CH4EnergyFinite) {
    auto f = make_ch4();
    cpu::Context ctx;
    DFTB2 dftb;
    dftb.sk_set = DftbSkSet::ThreeOb;
    dftb.with_forces = false;
    dftb.scc_tol = 1e-8;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});
    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
    EXPECT_NEAR(f.meta.get(kEnergy), dftb.last_energy(), 1e-12);
    // Bound organic: total energy well below zero (eV).
    EXPECT_LT(dftb.last_energy(), -50.0);
}

TEST(DFTB23ob, UnsupportedElementThrows) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 26;
    f.box = cpu::Box::cubic(10.0);
    f.meta.get(kEnergy) = 0;
    cpu::Context ctx;
    DFTB2 dftb;
    dftb.sk_set = DftbSkSet::ThreeOb;
    dftb.build({ctx, f});
    EXPECT_THROW(dftb.launch({ctx, f, 0, 0}), std::runtime_error);
}

TEST(DFTB23ob, DistinctFromMioSketch) {
    auto f1 = make_ch4();
    auto f2 = make_ch4();
    cpu::Context ctx;
    DFTB2 a;
    a.sk_set = DftbSkSet::ThreeOb;
    a.with_forces = false;
    a.build({ctx, f1});
    a.launch({ctx, f1, 0, 0});
    DFTB2 b;
    b.sk_set = DftbSkSet::MioSketch;
    b.with_forces = false;
    b.build({ctx, f2});
    b.launch({ctx, f2, 0, 0});
    EXPECT_GT(std::fabs(a.last_energy() - b.last_energy()), 0.5);
}

namespace {

cpu::Frame make_water_offeq() {
    // Off-eq water (Å) so |F| ≠ 0 — same geometry as dftb2-analytic-forces.
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

double energy_3ob_only(cpu::Frame g) {
    g.meta.get(kEnergy) = 0;
    for (int i = 0; i < g.nlocal(); ++i)
        g.span(kAtom_fx)[i] = g.span(kAtom_fy)[i] = g.span(kAtom_fz)[i] = 0;
    cpu::Context ctx;
    DFTB2 d;
    d.sk_set = DftbSkSet::ThreeOb;
    d.with_forces = false;
    d.scc_tol = 1e-10;
    d.max_scc = 80;
    d.build({ctx, g});
    d.launch({ctx, g, 0, 0});
    return g.meta.get(kEnergy);
}

} // namespace

TEST(DFTB23ob, AnalyticForcesMatchFdWater) {
    auto f = make_water_offeq();
    cpu::Context ctx;
    DFTB2 dftb;
    dftb.sk_set = DftbSkSet::ThreeOb;
    dftb.with_forces = true;
    dftb.scc_tol = 1e-10;
    dftb.max_scc = 80;
    dftb.build({ctx, f});
    dftb.launch({ctx, f, 0, 0});

    EXPECT_TRUE(std::isfinite(dftb.last_energy()));
    double f2 = 0;
    double net[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        const double fxi = static_cast<double>(f.span(kAtom_fx)[i]);
        const double fyi = static_cast<double>(f.span(kAtom_fy)[i]);
        const double fzi = static_cast<double>(f.span(kAtom_fz)[i]);
        EXPECT_TRUE(std::isfinite(fxi));
        EXPECT_TRUE(std::isfinite(fyi));
        EXPECT_TRUE(std::isfinite(fzi));
        f2 += fxi * fxi + fyi * fyi + fzi * fzi;
        net[0] += fxi;
        net[1] += fyi;
        net[2] += fzi;
    }
    EXPECT_GT(f2, 1e-8);
    for (int d = 0; d < 3; ++d)
        EXPECT_NEAR(net[d], 0.0, 1e-6);

    constexpr double h = 1e-4;
    // A converged SCC-DFTB2 gradient is stationary, so fixed-Δq analytic is
    // FD-exact to the central-difference truncation error. The old 1.0 eV/Å
    // "charge response" gate was masking the γ-gradient sign bug (measured
    // residual 8.4e-1 → 4.4e-7 eV/Å when fixed, 2026-07-28).
    constexpr double tol = 1e-5;

    {
        auto fp = make_water_offeq();
        fp.span(kAtom_x)[0] += static_cast<StorageReal>(h);
        auto fm = make_water_offeq();
        fm.span(kAtom_x)[0] -= static_cast<StorageReal>(h);
        const double fd = -(energy_3ob_only(fp) - energy_3ob_only(fm)) / (2.0 * h);
        const double fx0 = static_cast<double>(f.span(kAtom_fx)[0]);
        EXPECT_NEAR(fx0, fd, tol) << "F_O_x analytic=" << fx0 << " fd=" << fd;
    }
    {
        auto fp = make_water_offeq();
        fp.span(kAtom_y)[1] += static_cast<StorageReal>(h);
        auto fm = make_water_offeq();
        fm.span(kAtom_y)[1] -= static_cast<StorageReal>(h);
        const double fd = -(energy_3ob_only(fp) - energy_3ob_only(fm)) / (2.0 * h);
        const double fy1 = static_cast<double>(f.span(kAtom_fy)[1]);
        EXPECT_NEAR(fy1, fd, tol) << "F_H1_y analytic=" << fy1 << " fd=" << fd;
    }
}
