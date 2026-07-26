/// @file tests/cpu/semiempirical/test_nddo_api_surface.cpp
/// Public NDDO surface is AM1/RM1 methods + evaluators, not NDDOSolver bag.
#include <fstream>
#include <gtest/gtest.h>
#include <string>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/semiempirical/am1_evaluator.h"
#include "atomiverse/cpu/semiempirical/am1_solver.h"
#include "atomiverse/cpu/semiempirical/rm1_solver.h"

using namespace ATV;
using namespace ATV::cpu::semiempirical;

TEST(NDDOApiSurface, AM1MembersAndRun) {
    cpu::Frame f(2);
    f.span(kAtom_type)[0] = f.span(kAtom_type)[1] = 1;
    f.span(kAtom_x)[0] = 0;
    f.span(kAtom_y)[0] = 0;
    f.span(kAtom_z)[0] = 0;
    f.span(kAtom_x)[1] = 0.74;
    f.span(kAtom_y)[1] = 0;
    f.span(kAtom_z)[1] = 0;
    f.span(kAtom_mass)[0] = f.span(kAtom_mass)[1] = 1.008;
    f.box = cpu::Box::cubic(20.0);

    AM1 am1;
    am1.charge = 0;
    am1.max_iter = 100;
    am1.build(f);
    const auto& r = am1.run(f);
    EXPECT_TRUE(r.converged);
    EXPECT_TRUE(std::isfinite(r.total_energy_ev));
}

TEST(NDDOApiSurface, PublicHeaderTombstonesNDDOSolver) {
    // Grep-style gate: public nddo_solver.h is an #error tombstone.
    std::ifstream in("include/atomiverse/cpu/semiempirical/nddo_solver.h");
    if (!in) {
        // Running from build dir — try source relative path via macro.
        GTEST_SKIP() << "source tree path not available from test cwd";
    }
    std::string line;
    bool has_error = false;
    while (std::getline(in, line)) {
        if (line.find("#error") != std::string::npos)
            has_error = true;
    }
    EXPECT_TRUE(has_error);
}

TEST(NDDOApiSurface, AM1EvaluatorDefaultConstructible) {
    AM1Evaluator ev;
    ev.method().max_iter = 50;
    EXPECT_EQ(ev.method().max_iter, 50);
    EXPECT_TRUE(ev.gradient_enabled());
}
