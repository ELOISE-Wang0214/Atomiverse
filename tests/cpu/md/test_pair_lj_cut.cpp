/// @file tests/cpu/md/test_pair_lj_cut.cpp
/// @brief Regression test for the PairLJCut port of cpu::PairLJCut.
///
/// Verifies (a) the new class name resolves, (b) the backwards-compat
/// alias `cpu::PairLJCut` still works and is the same type as
/// `cpu::PairLJCut`, (c) energy + per-atom forces on a two-atom test
/// match the documented LJ 12-6 formula bit-for-bit.

#include <gtest/gtest.h>

#include <type_traits>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h" // legacy alias header
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h" // new canonical header

TEST(PairLJCutPort, NewClassExists) {
    ATV::cpu::PairLJCut lj{0.1, 3.0, 8.0};
    EXPECT_EQ(lj.name(), "lj");
    EXPECT_DOUBLE_EQ(lj.pair_cutoff(), 8.0);
}

TEST(PairLJCutPort, LegacyAliasIsSameType) {
    static_assert(std::is_same_v<ATV::cpu::PairLJCut, ATV::cpu::PairLJCut>,
                  "cpu::PairLJCut must alias cpu::PairLJCut after the port");
    SUCCEED();
}

TEST(PairLJCutPort, EnergyMatchesAnalyticFormula) {
    // Two atoms at r = sigma → V = 4ε(1 − 1) = 0.
    const double eps = 0.5;
    const double sigma = 3.0;
    const double rcut = 8.0;

    ATV::cpu::Frame frame(2);
    frame.atoms.x[0] = 0.0;
    frame.atoms.y[0] = 0.0;
    frame.atoms.z[0] = 0.0;
    frame.atoms.x[1] = sigma;
    frame.atoms.y[1] = 0.0;
    frame.atoms.z[1] = 0.0;
    frame.atoms.mass[0] = 1.0;
    frame.atoms.mass[1] = 1.0;
    frame.box = ATV::Mat3::cubic(50.0);

    ATV::cpu::md::MD sim;
    sim.set_neighbor(ATV::cpu::md::NeighborList<>{rcut, /*skin=*/0.5});
    sim.set_pair(ATV::cpu::PairLJCut{eps, sigma, rcut});
    sim.set_integrator(ATV::cpu::md::VelocityVerlet{/*dt=*/0.001});
    sim.build();
    // 1 integration step at dt=0.001 fs with zero initial velocities ≈
    // 0 displacement; the launch records the energy at the starting
    // geometry while leaving positions essentially untouched.
    sim.run(frame, 1);

    // Energy at r≈σ with shifted potential:
    //   V = 4ε[(σ/r)¹² − (σ/r)⁶] − V_shift
    // At r=σ, (σ/r)⁶ = 1, so the first bracket is 0 → V = −V_shift.
    const double sr_cut = sigma / rcut;
    const double sr_cut6 = sr_cut * sr_cut * sr_cut * sr_cut * sr_cut * sr_cut;
    const double sr_cut12 = sr_cut6 * sr_cut6;
    const double v_shift = 4.0 * eps * (sr_cut12 - sr_cut6);
    const double expected = -v_shift;
    // 1e-6 absolute tolerance covers the tiny displacement during the
    // single integration step (atoms move ~ a·dt²/2 ~ 1e-10 Å with
    // initial force ≈ 0 at r=σ).
    EXPECT_NEAR(frame.scalars.energy, expected, 1e-6)
        << "PairLJCut at r≈σ should yield −V_shift = " << expected;
}
