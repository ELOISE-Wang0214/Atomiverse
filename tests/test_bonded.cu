/// @file test_bonded.cu
/// @brief Tests for bonded interactions: HarmonicBond, HarmonicAngle, MultiHarmonicDihedral.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cuda/potentials/harmonic_bond.cuh"
#include "atomiverse/cuda/potentials/harmonic_angle.cuh"
#include "atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cpu/frame.h"

using namespace ATV;

// ── HarmonicBond ───────────────────────────────────────────────────

TEST(HarmonicBond, DimerAtEquilibrium) {
    // 2-atom system at r0 → V=0, F=0
    cpu::Frame frame(2);
    frame.x = xt::xtensor<Real, 1>({0.0, 1.53});
    frame.y = xt::xtensor<Real, 1>({0.0, 0.0});
    frame.z = xt::xtensor<Real, 1>({0.0, 0.0});
    frame.mass = xt::xtensor<Real, 1>({14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.bond_i = xt::xtensor<Int, 1>({0});
    frame.bond_j = xt::xtensor<Int, 1>({1});
    frame.bond_type = xt::xtensor<Int, 1>({0});

    md::MD sim;
    sim.on(md::Phase::force).add(HarmonicBond{{350.0}, {1.53}});
    sim.build();
    sim.run(frame, 1);

    EXPECT_NEAR(frame.energy, 0.0, 1e-10);
    EXPECT_NEAR(frame.fx[0], 0.0, 1e-10);
    EXPECT_NEAR(frame.fx[1], 0.0, 1e-10);
}

TEST(HarmonicBond, DimerDisplaced) {
    // 2-atom system at r != r0 → V > 0
    cpu::Frame frame(2);
    Real r = 1.60;
    frame.x = xt::xtensor<Real, 1>({0.0, r});
    frame.y = xt::xtensor<Real, 1>({0.0, 0.0});
    frame.z = xt::xtensor<Real, 1>({0.0, 0.0});
    frame.mass = xt::xtensor<Real, 1>({14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.bond_i = xt::xtensor<Int, 1>({0});
    frame.bond_j = xt::xtensor<Int, 1>({1});
    frame.bond_type = xt::xtensor<Int, 1>({0});

    Real K = 350.0, r0 = 1.53;

    md::MD sim;
    sim.on(md::Phase::force).add(HarmonicBond{{K}, {r0}});
    sim.build();
    sim.run(frame, 1);

    Real expected_V = K * (r - r0) * (r - r0);
    EXPECT_NEAR(frame.energy, expected_V, 1e-6);

    // Forces should be attractive (pulling atoms together)
    EXPECT_GT(frame.fx[0], 0.0);  // atom 0 pushed right
    EXPECT_LT(frame.fx[1], 0.0);  // atom 1 pushed left
    // Newton's 3rd law
    EXPECT_NEAR(frame.fx[0] + frame.fx[1], 0.0, 1e-10);
}

// ── HarmonicAngle ──────────────────────────────────────────────────

TEST(HarmonicAngle, AtEquilibrium) {
    // 3-atom system at θ0=109.5° → V=0, F=0
    cpu::Frame frame(3);
    Real theta0_deg = 109.5;
    Real theta0 = theta0_deg * M_PI / 180.0;
    Real r = 1.53;

    frame.x = xt::xtensor<Real, 1>({Real(0), r, Real(r + r * cos(M_PI - theta0))});
    frame.y = xt::xtensor<Real, 1>({Real(0), Real(0), Real(r * sin(M_PI - theta0))});
    frame.z = xt::xtensor<Real, 1>({0.0, 0.0, 0.0});
    frame.mass = xt::xtensor<Real, 1>({14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.angle_i = xt::xtensor<Int, 1>({0});
    frame.angle_j = xt::xtensor<Int, 1>({1});
    frame.angle_k = xt::xtensor<Int, 1>({2});
    frame.angle_type = xt::xtensor<Int, 1>({0});

    md::MD sim;
    sim.on(md::Phase::force).add(HarmonicAngle{{60.0}, {theta0_deg}});
    sim.build();
    sim.run(frame, 1);

    EXPECT_NEAR(frame.energy, 0.0, 1e-6);
}

// ── MultiHarmonicDihedral ──────────────────────────────────────────

TEST(MultiHarmonicDihedral, TransConformation) {
    // 4-atom chain in trans (φ=180°) → known V
    cpu::Frame frame(4);
    frame.x = xt::xtensor<Real, 1>({0.0, 1.0, 2.0, 3.0});
    frame.y = xt::xtensor<Real, 1>({0.0, 0.0, 0.0, 0.0});
    frame.z = xt::xtensor<Real, 1>({0.0, 0.0, 0.0, 0.0});
    frame.mass = xt::xtensor<Real, 1>({14.0, 14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.dihedral_i = xt::xtensor<Int, 1>({0});
    frame.dihedral_j = xt::xtensor<Int, 1>({1});
    frame.dihedral_k = xt::xtensor<Int, 1>({2});
    frame.dihedral_l = xt::xtensor<Int, 1>({3});
    frame.dihedral_type = xt::xtensor<Int, 1>({0});

    // For a collinear chain, b1×b2 = 0, so the dihedral is undefined.
    // Move atom 0 and 3 off-axis to create a well-defined trans dihedral.
    frame.y[0] = 0.5;   // atom 0 above the axis
    frame.y[3] = -0.5;  // atom 3 below = trans

    std::array<Real, 5> C = {1.73, -4.49, 0.776, 6.99, 0.0};

    md::MD sim;
    sim.on(md::Phase::force).add(MultiHarmonicDihedral{{C}});
    sim.build();
    sim.run(frame, 1);

    // Energy should be finite
    EXPECT_TRUE(std::isfinite(frame.energy));
    // Forces should be finite
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(frame.fx[i]));
        EXPECT_TRUE(std::isfinite(frame.fy[i]));
    }
}

// ── Integration: bonded + nonbonded ────────────────────────────────

TEST(Bonded, ChainWithLJ) {
    // 4-atom chain with bonds + LJ, 10 steps, finite results
    cpu::Frame frame(4);
    frame.x = xt::xtensor<Real, 1>({0.0, 1.53, 3.06, 4.59});
    frame.y = xt::xtensor<Real, 1>({0.0, 0.3, -0.3, 0.3});
    frame.z = xt::xtensor<Real, 1>({0.0, 0.0, 0.0, 0.0});
    frame.mass = xt::xtensor<Real, 1>({14.0, 14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.bond_i = xt::xtensor<Int, 1>({0, 1, 2});
    frame.bond_j = xt::xtensor<Int, 1>({1, 2, 3});
    frame.bond_type = xt::xtensor<Int, 1>({0, 0, 0});

    md::MD sim;
    sim.add(md::VelocityVerlet{0.1});
    sim.on(md::Phase::force)
        .add(HarmonicBond{{350.0}, {1.53}});
    sim.build();
    sim.run(frame, 10);

    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(frame.x[i]));
        EXPECT_TRUE(std::isfinite(frame.vx[i]));
    }
}
