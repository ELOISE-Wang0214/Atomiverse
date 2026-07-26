/// @file tests/cpu/gpw/test_gpw_bomd.cpp
/// @brief Unit gate: MD binds GPWEvaluator forces for short BOMD steps.
///
/// One method surface only — public MD + VelocityVerlet + GPWEvaluator.
/// Full NVE energy-drift pin lives in regressions/gpw/gpw-bomd-water (e2e).

#include <gtest/gtest.h>

#include <cmath>
#include <utility>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/velocity_verlet.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/gpw/gaussian_basis.h"
#include "atomiverse/cpu/gpw/gpw_evaluator.h"
#include "atomiverse/cpu/gpw/gth_potential.h"

using namespace ATV;

namespace {

cpu::Frame make_water_cheap() {
    // Same geometry as gpw-forces-water / CP2K getting-started (Å).
    cpu::Frame f(3);
    const int z[3] = {8, 1, 1};
    const double x[3] = {5.0000, 5.7586, 4.2414};
    const double y[3] = {5.0000, 5.0000, 5.0000};
    const double zc[3] = {5.0000, 5.5043, 5.5043};
    const double mass[3] = {15.999, 1.008, 1.008};
    for (int i = 0; i < 3; ++i) {
        f.span(kAtom_type)[i] = z[i];
        f.span(kAtom_mass)[i] = mass[i];
        f.span(kAtom_tag)[i] = i + 1;
        f.span(kAtom_x)[i] = x[i];
        f.span(kAtom_y)[i] = y[i];
        f.span(kAtom_z)[i] = zc[i];
        f.span(kAtom_vx)[i] = 0.0;
        f.span(kAtom_vy)[i] = 0.0;
        f.span(kAtom_vz)[i] = 0.0;
    }
    // Tiny COM-preserving velocity so KE is non-zero (O≈0 for COM).
    f.span(kAtom_vx)[1] = 0.001;
    f.span(kAtom_vx)[2] = -0.001;
    f.box = cpu::Box::cubic(10.0);
    return f;
}

} // namespace

TEST(GpwBomd, MdBindsEvaluatorForcesTwoSteps) {
    auto frame = make_water_cheap();

    const cpu::gpw::Kind O{
        molrs::Element::O,
        cpu::gpw::GaussianBasis(molrs::Element::O, "DZVP-MOLOPT-GTH"),
        cpu::gpw::GthPotential(molrs::Element::O, "GTH-PBE-q6")};
    const cpu::gpw::Kind H{
        molrs::Element::H,
        cpu::gpw::GaussianBasis(molrs::Element::H, "DZVP-MOLOPT-GTH"),
        cpu::gpw::GthPotential(molrs::Element::H, "GTH-PBE-q1")};

    // Cheap settings: unit test only (e2e pin in regression).
    cpu::gpw::GPWEvaluator eval(
        /*kinds*/ {O, H}, /*cutoff_ry*/ 200.0, /*max_iter*/ 40, /*tol*/ 1e-5,
        /*mix*/ 0.3, /*xc*/ cpu::gpw::XcKind::PBE, /*nonperiodic*/ true,
        /*compute_forces*/ true);

    cpu::md::MD sim;
    sim.set_integrator(cpu::md::VelocityVerlet{0.5});
    sim.on(md::Phase::force).add(std::move(eval));
    sim.enable_timing(false);
    sim.build();
    ASSERT_EQ(sim.run(frame, 2), 2);

    const double pe = frame.meta.get(kEnergy);
    const double ke = frame.meta.get(kKinetic);
    EXPECT_TRUE(std::isfinite(pe)) << "PE=" << pe;
    EXPECT_TRUE(std::isfinite(ke)) << "KE=" << ke;
    EXPECT_GT(ke, 1e-12);

    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(std::isfinite(static_cast<double>(frame.span(kAtom_x)[i])));
        EXPECT_TRUE(std::isfinite(static_cast<double>(frame.span(kAtom_y)[i])));
        EXPECT_TRUE(std::isfinite(static_cast<double>(frame.span(kAtom_z)[i])));
        EXPECT_TRUE(
            std::isfinite(static_cast<double>(frame.span(kAtom_fx)[i])));
        EXPECT_TRUE(
            std::isfinite(static_cast<double>(frame.span(kAtom_fy)[i])));
        EXPECT_TRUE(
            std::isfinite(static_cast<double>(frame.span(kAtom_fz)[i])));
    }
}
