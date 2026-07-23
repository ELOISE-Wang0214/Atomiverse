/// @file tests/cpu/planewave/test_pw_method.cpp
/// Concrete PW method API: members + build/run(frame), no PWConfig bag.
#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/planewave/pw_solver.h"

using namespace ATV;
using namespace ATV::cpu::pw;

TEST(PWMethod, DefaultMembersConstructible) {
    PW pw;
    pw.ecut_ry = 12.0;
    pw.k_mesh = 2;
    pw.tol = 1e-5;
    pw.max_iter = 30;
    pw.mix_beta = 0.5;
    pw.pulay_depth = 4;
    EXPECT_DOUBLE_EQ(pw.ecut_ry, 12.0);
    EXPECT_EQ(pw.k_mesh, 2);
    EXPECT_TRUE(pw.upf_path.empty());
}

TEST(PWMethod, RunToyUpfFiniteEnergy) {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = 14;
    f.span(kAtom_type)[1] = 14;
    f.span(kAtom_mass)[0] = 28.0855;
    f.span(kAtom_mass)[1] = 28.0855;
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 1.3575;
    f.span(kAtom_y)[1] = 1.3575;
    f.span(kAtom_z)[1] = 1.3575;
    f.simbox = cpu::Box::cubic(5.43);

    PW pw;
    pw.ecut_ry = 10.0;
    pw.k_mesh = 2;
    pw.tol = 1e-6;
    pw.max_iter = 20;
    pw.build(f);
    const auto& r = pw.run(f);
    EXPECT_TRUE(std::isfinite(r.energy));
    EXPECT_TRUE(std::isfinite(r.energy_ha));
    EXPECT_EQ(&r, &pw.last_result());
}
