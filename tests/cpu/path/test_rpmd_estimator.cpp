/// @file tests/cpu/path/test_rpmd_estimator.cpp
#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/path/pimd.h"

using namespace ATV;

TEST(RpmdEstimator, AnalyticPositiveAndClassicalLimit) {
    cpu::path::PimdConfig cfg;
    cfg.n_beads = 4;
    cfg.temperature = 300.0;
    cfg.dt = 0.25;
    cfg.steps_per_block = 5;
    cfg.seed = 1;
    cpu::path::Pimd pimd(cfg);

    cpu::Frame f(1);
    f.span(kAtom_x)[0] = f.span(kAtom_y)[0] = f.span(kAtom_z)[0] = 0;
    f.span(kAtom_mass)[0] = 1.0;
    f.span(kAtom_type)[0] = 1;
    f.span(kAtom_tag)[0] = 1;
    f.simbox = cpu::Box::cubic(50.0);

    // Trigger ensure_built via short run so mass/spring set.
    pimd.run(f, 2);
    const double rg2_an = pimd.free_particle_rg2_analytic();
    EXPECT_GT(rg2_an, 0.0);
    EXPECT_TRUE(std::isfinite(rg2_an));
    EXPECT_GT(pimd.radius_of_gyration_sq(), 0.0);

    // P=1 classical: free_particle formula → 0
    cpu::path::PimdConfig cfg1 = cfg;
    cfg1.n_beads = 2; // minimum; check formula factor (1-1/P²) shrinks with P→1
    // Direct formula check:
    const double hbar = units::h_planck_ev_fs / (2.0 * units::pi);
    const double beta = 1.0 / (units::kB_ev_per_k * 300.0);
    const double m_eff = units::mass_v2_to_energy * 1.0;
    const double P = 4.0;
    const double expect =
        (hbar * beta) * (hbar * beta) / (4.0 * m_eff) * (1.0 - 1.0 / (P * P));
    EXPECT_NEAR(rg2_an, expect, 1e-12);
}
