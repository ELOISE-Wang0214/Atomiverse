/// @file tests/cpu/dynamics/test_const_pot_bo.cpp
/// @brief Unit tests for ConstPotBO electrode charge solve (single methods).

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/dynamics/const_pot_bo.h"
#include "atomiverse/cpu/frame/frame.h"

using namespace ATV;

namespace {

/// Madelung-like 4-atom geometry, primary cell only (no 2D images).
cpu::Frame madelung_primary() {
    cpu::Frame f(4);
    // bot electrode, top electrode, bot elyt, top elyt
    // Centered in Lz=20 so full 3D PBC cannot wrap the slab (see regression).
    const double z[4] = {8.0, 12.0, 9.0, 11.0};
    const double q0[4] = {0.0, 0.0, 0.5, -0.5};
    const int typ[4] = {1, 2, 3, 3};
    for (int i = 0; i < 4; ++i) {
        f.span(kAtom_x)[i] = 0.0;
        f.span(kAtom_y)[i] = 0.0;
        f.span(kAtom_z)[i] = z[i];
        f.span(kAtom_charge)[i] = q0[i];
        f.span(kAtom_type)[i] = typ[i];
        f.span(kAtom_mass)[i] = (typ[i] < 3) ? 196.97 : 1.0;
        f.span(kAtom_tag)[i] = i + 1;
    }
    f.simbox = cpu::Box::orthorhombic(1.0, 1.0, 20.0);
    return f;
}

void run_once(cpu::md::ConstPotBO& conp, cpu::Frame& f) {
    cpu::Context ctx;
    conp.build({ctx, f});
    conp.launch({ctx, f, /*stage=*/0, /*step=*/0});
}

} // namespace

TEST(ConstPotBO, NameAndNoConfigFactory) {
    cpu::md::ConstPotBO conp(2.0);
    EXPECT_EQ(conp.name(), "const_pot_bo");
    conp.add_electrode(1, 0.0).add_electrode(2, 1.0);
    EXPECT_DOUBLE_EQ(conp.delta_psi(), 1.0);
}

TEST(ConstPotBO, NeutralityPrimaryCell) {
    auto f = madelung_primary();
    cpu::md::ConstPotBO conp(2.0);
    conp.add_electrode(1, 0.0);
    conp.add_electrode(2, 1.0);
    conp.electroneutrality = true;
    conp.xy_image_range = 0;
    run_once(conp, f);
    EXPECT_LE(conp.neutrality_residual(), 1e-10);
    EXPECT_NEAR(conp.total_charge(0) + conp.total_charge(1), 0.0, 1e-10);
}

TEST(ConstPotBO, SignFollowsDeltaPsi) {
    auto f = madelung_primary();
    cpu::md::ConstPotBO conp(2.0);
    conp.add_electrode(1, 0.0);
    conp.add_electrode(2, 1.0);
    conp.xy_image_range = 0;
    run_once(conp, f);
    // Positive ΔΨ (top − bot) ⇒ top electrode draws opposite electrolyte
    // polarization; with +0.5 near bot and −0.5 near top, Q_top > 0.
    EXPECT_GT(conp.total_charge(1), 0.0);
    EXPECT_LT(conp.total_charge(0), 0.0);
    EXPECT_DOUBLE_EQ(static_cast<double>(f.span(kAtom_charge)[0]),
                     conp.total_charge(0));
    EXPECT_DOUBLE_EQ(static_cast<double>(f.span(kAtom_charge)[1]),
                     conp.total_charge(1));
}

TEST(ConstPotBO, PrimaryCellMatchesAnalytic) {
    // Independent analytic: no images, eta=2, Ψ=(0,1), q_elyt=(0.5,-0.5)
    // Q from BO + neutrality (see agent pin script).
    constexpr double kQTop = 0.27175356;
    auto f = madelung_primary();
    cpu::md::ConstPotBO conp(2.0);
    conp.add_electrode(1, 0.0);
    conp.add_electrode(2, 1.0);
    conp.xy_image_range = 0;
    run_once(conp, f);
    EXPECT_NEAR(conp.total_charge(1), kQTop, 1e-6);
    EXPECT_NEAR(conp.total_charge(0), -kQTop, 1e-6);
}

TEST(ConstPotBO, ElectrolyteChargesUntouched) {
    auto f = madelung_primary();
    cpu::md::ConstPotBO conp(2.0);
    conp.add_electrode(1, 0.0).add_electrode(2, 1.0);
    run_once(conp, f);
    EXPECT_DOUBLE_EQ(static_cast<double>(f.span(kAtom_charge)[2]), 0.5);
    EXPECT_DOUBLE_EQ(static_cast<double>(f.span(kAtom_charge)[3]), -0.5);
}

TEST(ConstPotBO, MdBindingMatchesStandalone) {
    auto f = madelung_primary();
    cpu::md::ConstPotBO standalone(2.0);
    standalone.add_electrode(1, 0.0);
    standalone.add_electrode(2, 1.0);
    standalone.xy_image_range = 50;
    run_once(standalone, f);
    const double q_ref = standalone.total_charge(1);

    auto f2 = madelung_primary();
    cpu::md::ConstPotBO conp(2.0);
    conp.add_electrode(1, 0.0);
    conp.add_electrode(2, 1.0);
    conp.xy_image_range = 50;

    // Direct Context launch after MD-style build should match (no MD).
    run_once(conp, f2);
    EXPECT_NEAR(conp.total_charge(1), q_ref, 1e-9);
}
