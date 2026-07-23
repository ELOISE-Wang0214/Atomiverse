/// @file tests/cpu/path/test_pimd.cpp
/// @brief Unit tests for cpu::path::Pimd skeleton.

#include <gtest/gtest.h>

#include <cmath>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/path/pimd.h"

using namespace ATV;

namespace {

cpu::Frame one_particle(double mass = 1.0) {
    cpu::Frame f(1);
    f.span(kAtom_x)[0] = 0.0;
    f.span(kAtom_y)[0] = 0.0;
    f.span(kAtom_z)[0] = 0.0;
    f.span(kAtom_vx)[0] = 0.0;
    f.span(kAtom_vy)[0] = 0.0;
    f.span(kAtom_vz)[0] = 0.0;
    f.span(kAtom_mass)[0] = mass;
    f.span(kAtom_type)[0] = 0;
    f.span(kAtom_tag)[0] = 1;
    f.simbox = cpu::Box::cubic(50.0);
    return f;
}

} // namespace

TEST(Pimd, RejectsSingleBead) {
    cpu::path::PimdConfig cfg;
    cfg.n_beads = 1;
    EXPECT_THROW(cpu::path::Pimd{cfg}, std::invalid_argument);
}

TEST(Pimd, SpringConstantPositive) {
    cpu::path::PimdConfig cfg;
    cfg.n_beads = 4;
    cfg.temperature = 300.0;
    cpu::path::Pimd pimd(cfg);
    auto f = one_particle(1.0);
    pimd.run(f, 1);
    EXPECT_GT(pimd.spring_constant(), 0.0);
    EXPECT_EQ(pimd.n_beads(), 4);
}

TEST(Pimd, FreeParticleCentroidFinite) {
    cpu::path::PimdConfig cfg;
    cfg.n_beads = 4;
    cfg.temperature = 300.0;
    cfg.dt = 0.25;
    cfg.steps_per_block = 5;
    cfg.seed = 11;
    cfg.k_external = 0.0; // free particle
    cpu::path::Pimd pimd(cfg);
    auto f = one_particle(1.0);
    pimd.run(f, 10);
    EXPECT_TRUE(std::isfinite(f.span(kAtom_x)[0]));
    EXPECT_TRUE(std::isfinite(f.span(kAtom_y)[0]));
    EXPECT_TRUE(std::isfinite(f.span(kAtom_z)[0]));
    EXPECT_EQ(f.nlocal(), 1);
}

TEST(Pimd, HarmonicTrapRunsFinite) {
    // External trap is optional scaffolding; skeleton checks it does not
    // crash and leaves a finite centroid (NVE free-particle + trap without
    // thermostat is not a production thermostat path).
    cpu::path::PimdConfig cfg;
    cfg.n_beads = 4;
    cfg.temperature = 50.0;
    cfg.dt = 0.05;
    cfg.steps_per_block = 5;
    cfg.seed = 3;
    cfg.k_external = 1.0; // eV/Å²
    cpu::path::Pimd pimd(cfg);
    auto f = one_particle(1.0);
    pimd.run(f, 5);
    EXPECT_TRUE(std::isfinite(f.span(kAtom_x)[0]));
    EXPECT_TRUE(std::isfinite(f.span(kAtom_y)[0]));
    EXPECT_TRUE(std::isfinite(f.span(kAtom_z)[0]));
}

TEST(Pimd, BeadSpreadGrowsFromIdenticalStart) {
    cpu::path::PimdConfig cfg;
    cfg.n_beads = 8;
    cfg.temperature = 300.0;
    cfg.dt = 0.2;
    cfg.steps_per_block = 10;
    cfg.seed = 99;
    cpu::path::Pimd pimd(cfg);
    auto f = one_particle(1.0);
    pimd.run(f, 5);
    // Thermal velocities + springs → finite polymer size.
    EXPECT_GT(pimd.mean_bead_spread(), 0.0);
    EXPECT_TRUE(std::isfinite(pimd.mean_bead_spread()));
}
