/// @file tests/cpu/md/test_step_timing.cpp
/// @brief Timer subsystem: component classifier + opt-in capture through the
/// Backend's TimerRegistry (ac-001, cpu-md-lammps-perf-assessment).
///
/// Single-rank, bonded-only (no MPI halo) so it is safe under bare ctest.

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "atomiverse/core/mat3.h"
#include "atomiverse/core/timer.h"
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/step_timing.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/harmonic_angle.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
#include "atomiverse/cpu/potentials/multi_harmonic_dihedral.h"

using namespace ATV;
using cpu::md::StepTiming;

namespace {

inline constexpr double kKcalToEv = 1.0 / 23.060547830619026;

cpu::Frame make_chain(int n, double box_len) {
    cpu::Frame frame(n);
    const double theta0 = 109.5 * kPI / 180.0;
    const double half_apex = (kPI - theta0) / 2.0;
    const double dx = 1.53 * std::cos(half_apex);
    const double dy = 1.53 * std::sin(half_apex);
    const double x0 = (box_len - (n - 1) * dx) / 2.0;
    for (int i = 0; i < n; ++i) {
        frame.atoms.x[i] = x0 + i * dx;
        frame.atoms.y[i] = box_len / 2.0 + (i % 2 ? dy : 0.0);
        frame.atoms.z[i] = box_len / 2.0;
        frame.atoms.mass[i] = 14.027;
        frame.atoms.type[i] = 0;
    }
    frame.atoms.map_clear();
    frame.atoms.map_set(0, frame.atoms.nlocal);
    for (int k = 0; k < n - 1; ++k)
        frame.atoms.add_bond(k, k + 1, 0);
    for (int k = 0; k < n - 2; ++k)
        frame.atoms.add_angle(k, k + 1, k + 2, 0);
    for (int k = 0; k < n - 3; ++k)
        frame.atoms.add_dihedral(k, k + 1, k + 2, k + 3, 0);
    frame.atoms.build_special_list();
    frame.box = Mat3::cubic(box_len);
    return frame;
}

void configure_sim(cpu::md::MD& sim) {
    sim.add(cpu::md::VelocityVerlet{0.5});
    sim.on(ATV::md::Phase::force)
        .add(HarmonicBond{{268.0 * kKcalToEv}, {1.53}})
        .add(HarmonicAngle{{60.0 * kKcalToEv}, {109.5 * kPI / 180.0}})
        .add(MultiHarmonicDihedral{{{1.73 * kKcalToEv, -4.49 * kKcalToEv,
                                     0.776 * kKcalToEv, 6.99 * kKcalToEv, 0.0}}});
}

} // namespace

// ── classify: every concrete CPU MD component name maps correctly ──────
TEST(StepTiming, ClassifyKnownComponentNames) {
    EXPECT_EQ(StepTiming::classify("lj"), StepTiming::Pair);
    EXPECT_EQ(StepTiming::classify("morse"), StepTiming::Pair);
    EXPECT_EQ(StepTiming::classify("harmonic_bond"), StepTiming::Bonded);
    EXPECT_EQ(StepTiming::classify("harmonic_angle"), StepTiming::Bonded);
    EXPECT_EQ(StepTiming::classify("multi_harmonic_dihedral"), StepTiming::Bonded);
    EXPECT_EQ(StepTiming::classify("neighbor_list"), StepTiming::Neigh);
    EXPECT_EQ(StepTiming::classify("halo_forward"), StepTiming::Comm);
    EXPECT_EQ(StepTiming::classify("halo_reverse"), StepTiming::Comm);
    EXPECT_EQ(StepTiming::classify("atom_migration"), StepTiming::Comm);
    EXPECT_EQ(StepTiming::classify("rcb_balancer"), StepTiming::Comm);
    EXPECT_EQ(StepTiming::classify("reduce_energy"), StepTiming::Reduce);
    EXPECT_EQ(StepTiming::classify("velocity_verlet"), StepTiming::Other);
    EXPECT_EQ(StepTiming::classify("nose_hoover_npt"), StepTiming::Other);
    EXPECT_EQ(StepTiming::classify("fix_langevin"), StepTiming::Other);
    EXPECT_EQ(StepTiming::classify("kinetic_energy"), StepTiming::Other);
}

// ── opt-in: timing off by default => registry stays empty ──────────────
TEST(StepTiming, OptInDisabledByDefault) {
    auto frame = make_chain(50, 100.0);
    cpu::md::MD sim;
    configure_sim(sim);
    sim.build();
    sim.run(frame, 20);
    const TimerRegistry* reg = sim.timing();
    // Either no backend timing yet, or enabled==false with zero steps.
    if (reg)
        EXPECT_EQ(reg->steps(), 0u);
}

// ── capture: per-category sum within +/-5% of total step-loop wall ─────
TEST(StepTiming, SumWithinFivePercentOfTotal) {
    auto frame = make_chain(500, 200.0);
    cpu::md::MD sim;
    configure_sim(sim);
    sim.enable_timing(true);
    sim.build();

    constexpr int kSteps = 300;
    sim.run(frame, kSteps);

    const TimerRegistry* reg = sim.timing();
    ASSERT_NE(reg, nullptr);
    EXPECT_TRUE(reg->enabled());
    EXPECT_EQ(reg->steps(), static_cast<std::uint64_t>(kSteps));
    EXPECT_GT(reg->loop_seconds(), 0.0);

    const auto roll = reg->rollup_by_category();
    ASSERT_GE(roll.size(), static_cast<std::size_t>(StepTiming::Other) + 1);
    EXPECT_GT(roll[StepTiming::Bonded].total_s, 0.0);
    EXPECT_GT(roll[StepTiming::Other].total_s, 0.0);

    double sum = 0.0;
    for (const auto& s : roll)
        sum += s.total_s;
    const double frac = sum / reg->loop_seconds();
    EXPECT_LE(frac, 1.0 + 1e-6) << "component sum exceeds loop wall";
    EXPECT_GE(frac, 0.95)
        << "instrumentation captured only " << (100.0 * frac) << "%";
}
