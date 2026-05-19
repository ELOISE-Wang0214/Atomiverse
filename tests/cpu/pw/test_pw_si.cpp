/// @file test_pw_si.cpp
/// @brief Integration test for pw::pw_scf().
#include <gtest/gtest.h>
#include <cmath>
#include "atomiverse/cpu/pw/pw_scf.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/mat3.h"

using namespace ATV;
using namespace ATV::cpu::pw;

namespace {
cpu::Frame make_si_fcc(double a) {
    cpu::Frame f(2);
    f.atoms.type[0] = 14; f.atoms.type[1] = 14;
    f.atoms.mass[0] = 28.0855; f.atoms.mass[1] = 28.0855;
    f.atoms.x[0] = 0;        f.atoms.y[0] = 0;        f.atoms.z[0] = 0;
    f.atoms.x[1] = a / 4.0;  f.atoms.y[1] = a / 4.0;  f.atoms.z[1] = a / 4.0;
    double h = a / 2.0;
    Mat3 cell;
    cell(0,0) = 0; cell(0,1) = h; cell(0,2) = h;
    cell(1,0) = h; cell(1,1) = 0; cell(1,2) = h;
    cell(2,0) = h; cell(2,1) = h; cell(2,2) = 0;
    f.box = cell;
    return f;
}
} // anonymous

TEST(PWScf, CompilesAndRuns) {
    auto frame = make_si_fcc(5.43);
    auto r = pw_scf(frame, "", 10.0, 2, 1e-6, 30);
    EXPECT_TRUE(std::isfinite(r.energy)) << "Energy should be finite, got: " << r.energy;
}

TEST(PWScf, RunsWithoutCrash) {
    auto frame = make_si_fcc(5.43);
    auto r = pw_scf(frame, "", 10.0, 2, 1e-6, 20);
    EXPECT_TRUE(std::isfinite(r.energy)) << "Energy should be finite, got: " << r.energy;
}
