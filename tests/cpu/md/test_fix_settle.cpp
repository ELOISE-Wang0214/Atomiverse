/// @file tests/cpu/md/test_fix_settle.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/dynamics/fix_settle.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

cpu::Frame tip3p_perturbed() {
    // Ideal TIP3P then perturb H atoms.
    const double d_oh = 0.9572;
    const double d_hh = 1.5139;
    const double half = 0.5 * d_hh;
    const double h = std::sqrt(d_oh * d_oh - half * half);
    cpu::Frame f(3);
    // O, H1, H2
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = h;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = -half;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_x)[2] = half;
    f.span(kAtom_y)[2] = 0;
    f.span(kAtom_z)[2] = 0;
    // Perturb
    f.span(kAtom_x)[1] += 0.05;
    f.span(kAtom_y)[2] += 0.03;
    f.span(kAtom_mass)[0] = 16.0;
    f.span(kAtom_mass)[1] = 1.0;
    f.span(kAtom_mass)[2] = 1.0;
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_type)[i] = i == 0 ? 1 : 2;
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_vx)[i] = 0.01 * (i + 1);
        f.span(kAtom_vy)[i] = -0.02 * (i + 1);
        f.span(kAtom_vz)[i] = 0.005;
    }
    f.simbox = cpu::Box::cubic(20.0);
    return f;
}

double dist(const cpu::Frame& f, int i, int j) {
    const double dx = f.span(kAtom_x)[i] - f.span(kAtom_x)[j];
    const double dy = f.span(kAtom_y)[i] - f.span(kAtom_y)[j];
    const double dz = f.span(kAtom_z)[i] - f.span(kAtom_z)[j];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

} // namespace

TEST(FixSettle, RestoresGeometry) {
    auto f = tip3p_perturbed();
    cpu::md::FixSettle settle(2.0, 0.9572, 1.5139);
    settle.add_water(1, 2, 3);
    cpu::Context ctx;
    settle.build({ctx, f});
    settle.launch({ctx, f, ATV::slot_of(md::Phase::initial_integrate), 0});
    settle.launch({ctx, f, ATV::slot_of(md::Phase::post_integrate), 0});
    EXPECT_NEAR(dist(f, 0, 1), 0.9572, 1e-10);
    EXPECT_NEAR(dist(f, 0, 2), 0.9572, 1e-10);
    EXPECT_NEAR(dist(f, 1, 2), 1.5139, 1e-10);
}

TEST(FixSettle, VelocityConstraint) {
    auto f = tip3p_perturbed();
    cpu::md::FixSettle settle(2.0, 0.9572, 1.5139);
    settle.add_water(1, 2, 3);
    cpu::Context ctx;
    settle.build({ctx, f});
    settle.launch({ctx, f, ATV::slot_of(md::Phase::initial_integrate), 0});
    settle.launch({ctx, f, ATV::slot_of(md::Phase::post_integrate), 0});
    settle.launch({ctx, f, ATV::slot_of(md::Phase::final_integrate), 0});
    auto bond_rv = [&](int i, int j) {
        const double dx = f.span(kAtom_x)[i] - f.span(kAtom_x)[j];
        const double dy = f.span(kAtom_y)[i] - f.span(kAtom_y)[j];
        const double dz = f.span(kAtom_z)[i] - f.span(kAtom_z)[j];
        const double dvx = f.span(kAtom_vx)[i] - f.span(kAtom_vx)[j];
        const double dvy = f.span(kAtom_vy)[i] - f.span(kAtom_vy)[j];
        const double dvz = f.span(kAtom_vz)[i] - f.span(kAtom_vz)[j];
        return dx * dvx + dy * dvy + dz * dvz;
    };
    // One pass projection — should be much smaller than unconstrained.
    EXPECT_NEAR(bond_rv(0, 1), 0.0, 1e-6);
}
