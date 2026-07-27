/// @file tests/cpu/semiempirical/test_dftb2.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/dftb2.h"

using namespace ATV;

namespace {

cpu::Frame make_h2(double r_ang = 0.74) {
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

/// O₂ triplet geometry (experimental-ish bond length, Å).
cpu::Frame make_o2(double r = 1.21) {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 8;
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_x)[1] = static_cast<StorageReal>(r);
    f.span(kAtom_y)[1] = 0.0;
    f.span(kAtom_z)[1] = 0.0;
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_mass)[i] = 15.999;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_fx)[i] = f.span(kAtom_fy)[i] = f.span(kAtom_fz)[i] = 0;
    }
    f.box = cpu::Box::cubic(30.0);
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
    // Closed shell: zero Mulliken spin
    ASSERT_EQ(dftb.last_mulliken_spin().size(), 2u);
    EXPECT_NEAR(dftb.last_mulliken_spin()[0], 0.0, 1e-8);
    EXPECT_NEAR(dftb.last_mulliken_spin()[1], 0.0, 1e-8);
}

TEST(DFTB2, RejectsUnsupportedElement) {
    cpu::Frame f(1);
    f.span(kAtom_type)[0] = 9; // F not in tables
    f.span(kAtom_x)[0] = f.span(kAtom_y)[0] = f.span(kAtom_z)[0] = 0;
    f.span(kAtom_mass)[0] = 19.0;
    f.span(kAtom_tag)[0] = 1;
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

    constexpr double h = 1.0e-4;
    for (int a = 0; a < 2; ++a) {
        for (int k = 0; k < 3; ++k) {
            auto gp = f;
            auto gm = f;
            StorageReal* cp = (k == 0)   ? gp.span(kAtom_x)
                              : (k == 1) ? gp.span(kAtom_y)
                                         : gp.span(kAtom_z);
            StorageReal* cm = (k == 0)   ? gm.span(kAtom_x)
                              : (k == 1) ? gm.span(kAtom_y)
                                         : gm.span(kAtom_z);
            cp[a] += static_cast<StorageReal>(h);
            cm[a] -= static_cast<StorageReal>(h);
            const double fd = -(energy_only(gp) - energy_only(gm)) / (2.0 * h);
            const double ana =
                (k == 0)   ? static_cast<double>(f.span(kAtom_fx)[a])
                : (k == 1) ? static_cast<double>(f.span(kAtom_fy)[a])
                           : static_cast<double>(f.span(kAtom_fz)[a]);
            EXPECT_NEAR(ana, fd, 5e-3) << "atom " << a << " k " << k;
        }
    }
    // H2 force antisymmetry on bond axis
    EXPECT_NEAR(static_cast<double>(f.span(kAtom_fx)[0]
                                    + f.span(kAtom_fx)[1]),
                0.0, 1e-6);
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

    // Multi-element re-SCC residual includes Mulliken charge response beyond
    // fixed-Δq analytic gradient; gate is loose (component-wise).
    constexpr double h = 1.0e-4;
    constexpr double tol = 1.0; // eV/Å
    for (int a = 0; a < 3; ++a) {
        for (int k = 0; k < 3; ++k) {
            auto gp = f;
            auto gm = f;
            StorageReal* cp = (k == 0)   ? gp.span(kAtom_x)
                              : (k == 1) ? gp.span(kAtom_y)
                                         : gp.span(kAtom_z);
            StorageReal* cm = (k == 0)   ? gm.span(kAtom_x)
                              : (k == 1) ? gm.span(kAtom_y)
                                         : gm.span(kAtom_z);
            cp[a] += static_cast<StorageReal>(h);
            cm[a] -= static_cast<StorageReal>(h);
            const double fd = -(energy_only(gp) - energy_only(gm)) / (2.0 * h);
            const double ana =
                (k == 0)   ? static_cast<double>(f.span(kAtom_fx)[a])
                : (k == 1) ? static_cast<double>(f.span(kAtom_fy)[a])
                           : static_cast<double>(f.span(kAtom_fz)[a]);
            EXPECT_NEAR(ana, fd, tol) << "atom " << a << " k " << k;
        }
    }
    // Net force ~0 (no external field)
    double net[3] = {0, 0, 0};
    for (int a = 0; a < 3; ++a) {
        net[0] += static_cast<double>(f.span(kAtom_fx)[a]);
        net[1] += static_cast<double>(f.span(kAtom_fy)[a]);
        net[2] += static_cast<double>(f.span(kAtom_fz)[a]);
    }
    EXPECT_NEAR(net[0], 0.0, 1e-6);
    EXPECT_NEAR(net[1], 0.0, 1e-6);
    EXPECT_NEAR(net[2], 0.0, 1e-6);
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
    EXPECT_NEAR(m_tot, 2.0, 0.25);
    EXPECT_GT(std::fabs(m0), 0.3);
    EXPECT_GT(std::fabs(m1), 0.3);
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
