/// @file tests/cpu/planewave/test_pw_spin.cpp
/// @brief Unit tests for collinear spin-polarized PW members and O2 mag.
#include <cmath>
#include <stdexcept>

#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

using namespace ATV;
using namespace ATV::cpu::pw;

namespace {

cpu::Frame o2_frame(double L_ang) {
    cpu::Frame f(2);
    for (int i = 0; i < 2; ++i) {
        f.span(kAtom_type)[i] = 8;
        f.span(kAtom_mass)[i] = 15.999;
        f.span(kAtom_x)[i] = 0.0;
        f.span(kAtom_y)[i] = 0.0;
    }
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_z)[1] = 1.2075; // Å
    f.simbox = cpu::Box::cubic(L_ang);
    return f;
}

} // namespace

TEST(PWSpin, DefaultRestricted) {
    PW pw;
    EXPECT_FALSE(pw.spin_polarized);
    EXPECT_EQ(pw.multiplicity, 1);
    EXPECT_DOUBLE_EQ(pw.last_result().magnetization, 0.0);
}

TEST(PWSpin, MultiplicityParityThrows) {
    auto f = o2_frame(8.0);
    PW pw;
    pw.spin_polarized = true;
    pw.multiplicity = 2; // odd mag with even n_elec → invalid
    pw.ecut_ry = 8.0;
    pw.k_mesh = 1;
    pw.max_iter = 1;
    // empty upf_path = toy V_loc (valence_from_z)
    pw.build(f);
    EXPECT_THROW(pw.run(f), std::invalid_argument);
}

TEST(PWSpin, O2TripletMagnetizationIdentity) {
    // Fixed occupations: multiplicity 3 → M = 2 on PWResult (occupation
    // identity; does not require a tight SCF). Tiny ecut / short budget so
    // the unit suite stays in the seconds regime.
    auto f = o2_frame(8.0);
    PW pw;
    pw.ecut_ry = 6.0;
    pw.k_mesh = 1;
    pw.tol = 1e-3;
    pw.max_iter = 4;
    pw.mix_beta = 0.5;
    pw.spin_polarized = true;
    pw.multiplicity = 3;

    pw.build(f);
    const auto& r = pw.run(f);
    EXPECT_TRUE(std::isfinite(r.energy));
    EXPECT_TRUE(std::isfinite(r.energy_ha));
    EXPECT_NEAR(r.magnetization, 2.0, 1e-12);
}

TEST(PWSpin, O2TripletMultiKMagnetizationIdentity) {
    // k-point collinear UKS (pw-uks-kmesh-fe): fixed occupations still set
    // M = multiplicity − 1 with Monkhorst–Pack mesh ≥ 2. One SCF iter is
    // enough for the occupation identity; keep wall in the unit-test budget.
    auto f = o2_frame(5.5);
    PW pw;
    pw.ecut_ry = 4.0;
    pw.k_mesh = 2;
    pw.tol = 1e-1;
    pw.max_iter = 1;
    pw.mix_beta = 0.5;
    pw.spin_polarized = true;
    pw.multiplicity = 3;

    pw.build(f);
    const auto& r = pw.run(f);
    EXPECT_TRUE(std::isfinite(r.energy));
    EXPECT_TRUE(std::isfinite(r.energy_ha));
    EXPECT_NEAR(r.magnetization, 2.0, 1e-12);
}
