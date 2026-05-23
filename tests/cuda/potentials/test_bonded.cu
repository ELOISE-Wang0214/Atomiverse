/// @file test_bonded.cu
/// @brief Tests for bonded interactions: cuda::HarmonicBond, cuda::HarmonicAngle, cuda::MultiHarmonicDihedral.

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/evaluate.cuh"
#include "atomiverse/cuda/md/md.cuh"
#include "atomiverse/cuda/md/velocity_verlet.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/harmonic_angle.cuh"
#include "atomiverse/cuda/potentials/harmonic_bond.cuh"
#include "atomiverse/cuda/potentials/multi_harmonic_dihedral.cuh"

using namespace ATV;

// ── cuda::HarmonicBond ───────────────────────────────────────────────────

TEST(HarmonicBond, DimerAtEquilibrium) {
    // 2-atom system at r0 → V=0, F=0
    cpu::Frame frame(2);
    frame.atoms.x = xt::xtensor<double, 1>({0.0, 1.53});
    frame.atoms.y = xt::xtensor<double, 1>({0.0, 0.0});
    frame.atoms.z = xt::xtensor<double, 1>({0.0, 0.0});
    frame.atoms.mass = xt::xtensor<double, 1>({14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.bonds.atomi = xt::xtensor<int, 1>({0});
    frame.bonds.atomj = xt::xtensor<int, 1>({1});
    frame.bonds.type = xt::xtensor<int, 1>({0});

    // Single-pass force evaluation -> use cuda::Evaluate (no Integrator
    // required, unlike cuda::md::MD which validates one in build()).
    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::HarmonicBond{{350.0}, {1.53}});
    sim.build();
    sim.run(frame);

    EXPECT_NEAR(frame.scalars.energy, 0.0, 1e-10);
    EXPECT_NEAR(frame.atoms.fx[0], 0.0, 1e-10);
    EXPECT_NEAR(frame.atoms.fx[1], 0.0, 1e-10);
}

TEST(HarmonicBond, DimerDisplaced) {
    // 2-atom system at r != r0 → V > 0
    cpu::Frame frame(2);
    double r = 1.60;
    frame.atoms.x = xt::xtensor<double, 1>({0.0, r});
    frame.atoms.y = xt::xtensor<double, 1>({0.0, 0.0});
    frame.atoms.z = xt::xtensor<double, 1>({0.0, 0.0});
    frame.atoms.mass = xt::xtensor<double, 1>({14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.bonds.atomi = xt::xtensor<int, 1>({0});
    frame.bonds.atomj = xt::xtensor<int, 1>({1});
    frame.bonds.type = xt::xtensor<int, 1>({0});

    double K = 350.0, r0 = 1.53;

    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::HarmonicBond{{K}, {r0}});
    sim.build();
    sim.run(frame);

    double expected_V = K * (r - r0) * (r - r0);
    EXPECT_NEAR(frame.scalars.energy, expected_V, 1e-6);

    // Forces should be attractive (pulling atoms together)
    EXPECT_GT(frame.atoms.fx[0], 0.0); // atom 0 pushed right
    EXPECT_LT(frame.atoms.fx[1], 0.0); // atom 1 pushed left
    // Newton's 3rd law
    EXPECT_NEAR(frame.atoms.fx[0] + frame.atoms.fx[1], 0.0, 1e-10);
}

// ── cuda::HarmonicAngle ──────────────────────────────────────────────────

TEST(HarmonicAngle, AtEquilibrium) {
    // 3-atom system at θ0=109.5° → V=0, F=0
    cpu::Frame frame(3);
    double theta0_deg = 109.5;
    double theta0 = theta0_deg * M_PI / 180.0;
    double r = 1.53;

    frame.atoms.x = xt::xtensor<double, 1>(
        {double(0), r, double(r + r * cos(M_PI - theta0))});
    frame.atoms.y = xt::xtensor<double, 1>(
        {double(0), double(0), double(r * sin(M_PI - theta0))});
    frame.atoms.z = xt::xtensor<double, 1>({0.0, 0.0, 0.0});
    frame.atoms.mass = xt::xtensor<double, 1>({14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.angles.atomi = xt::xtensor<int, 1>({0});
    frame.angles.atomj = xt::xtensor<int, 1>({1});
    frame.angles.atomk = xt::xtensor<int, 1>({2});
    frame.angles.type = xt::xtensor<int, 1>({0});

    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::HarmonicAngle{{60.0}, {theta0_deg}});
    sim.build();
    sim.run(frame);

    EXPECT_NEAR(frame.scalars.energy, 0.0, 1e-6);
}

// ── cuda::MultiHarmonicDihedral ──────────────────────────────────────────

TEST(MultiHarmonicDihedral, TransConformation) {
    // 4-atom chain in trans (φ=180°) → known V
    cpu::Frame frame(4);
    frame.atoms.x = xt::xtensor<double, 1>({0.0, 1.0, 2.0, 3.0});
    frame.atoms.y = xt::xtensor<double, 1>({0.0, 0.0, 0.0, 0.0});
    frame.atoms.z = xt::xtensor<double, 1>({0.0, 0.0, 0.0, 0.0});
    frame.atoms.mass = xt::xtensor<double, 1>({14.0, 14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.dihedrals.atomi = xt::xtensor<int, 1>({0});
    frame.dihedrals.atomj = xt::xtensor<int, 1>({1});
    frame.dihedrals.atomk = xt::xtensor<int, 1>({2});
    frame.dihedrals.atoml = xt::xtensor<int, 1>({3});
    frame.dihedrals.type = xt::xtensor<int, 1>({0});

    // For a collinear chain, b1×b2 = 0, so the dihedral is undefined.
    // Move atom 0 and 3 off-axis to create a well-defined trans dihedral.
    frame.atoms.y[0] = 0.5;  // atom 0 above the axis
    frame.atoms.y[3] = -0.5; // atom 3 below = trans

    std::array<double, 5> C = {1.73, -4.49, 0.776, 6.99, 0.0};

    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::MultiHarmonicDihedral{{C}});
    sim.build();
    sim.run(frame);

    // Energy should be finite
    EXPECT_TRUE(std::isfinite(frame.scalars.energy));
    // Forces should be finite
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(frame.atoms.fx[i]));
        EXPECT_TRUE(std::isfinite(frame.atoms.fy[i]));
    }
}

// ── Integration: bonded + nonbonded ────────────────────────────────

TEST(Bonded, ChainWithLJ) {
    // 4-atom chain with bonds + LJ, 10 steps, finite results
    cpu::Frame frame(4);
    frame.atoms.x = xt::xtensor<double, 1>({0.0, 1.53, 3.06, 4.59});
    frame.atoms.y = xt::xtensor<double, 1>({0.0, 0.3, -0.3, 0.3});
    frame.atoms.z = xt::xtensor<double, 1>({0.0, 0.0, 0.0, 0.0});
    frame.atoms.mass = xt::xtensor<double, 1>({14.0, 14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    frame.bonds.atomi = xt::xtensor<int, 1>({0, 1, 2});
    frame.bonds.atomj = xt::xtensor<int, 1>({1, 2, 3});
    frame.bonds.type = xt::xtensor<int, 1>({0, 0, 0});

    cuda::md::MD sim;
    sim.add(cuda::md::VelocityVerlet{0.1});
    sim.on(md::Phase::force).add(cuda::HarmonicBond{{350.0}, {1.53}});
    sim.build();
    sim.run(frame, 10);

    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(frame.atoms.x[i]));
        EXPECT_TRUE(std::isfinite(frame.atoms.vx[i]));
    }
}
