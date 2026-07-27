/// @file test_bonded.cu
/// @brief Tests for bonded interactions: cuda::BondHarmonic, cuda::AngleHarmonic, cuda::DihedralMultiHarmonic.

#include <cmath>
#include <gtest/gtest.h>

#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cuda/dynamics/md.cuh"
#include "atomiverse/cuda/dynamics/velocity_verlet.cuh"
#include "atomiverse/cuda/evaluate/evaluate.cuh"
#include "atomiverse/cuda/neighbor/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/angle_harmonic.cuh"
#include "atomiverse/cuda/potentials/bond_harmonic.cuh"
#include "atomiverse/cuda/potentials/dihedral_multi_harmonic.cuh"

using namespace ATV;

// ── cuda::BondHarmonic ───────────────────────────────────────────────────

TEST(BondHarmonic, DimerAtEquilibrium) {
    // 2-atom system at r0 → V=0, F=0
    cpu::Frame frame(2);
    {
        auto x = frame.span(kAtom_x);
        auto y = frame.span(kAtom_y);
        auto z = frame.span(kAtom_z);
        auto mass = frame.span(kAtom_mass);
        auto tag = frame.span(kAtom_tag);
        x[0] = 0.0;
        x[1] = 1.53;
        y[0] = 0.0;
        y[1] = 0.0;
        z[0] = 0.0;
        z[1] = 0.0;
        mass[0] = 14.0;
        mass[1] = 14.0;
        tag[0] = 0;
        tag[1] = 1;
    }
    frame.box = cpu::Box::cubic(20.0);

    frame.clear_tag_index();
    frame.rebuild_tag_index(0, frame.nlocal());
    frame.add_bond(0, 1, 0);

    // Single-pass force evaluation -> use cuda::Evaluate (no Integrator
    // required, unlike cuda::md::MD which validates one in build()).
    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::BondHarmonic{{350.0}, {1.53}});
    sim.build();
    sim.run(frame);

    // CUDA kernels compute in real_t (single precision), so "zero" force
    // at equilibrium carries float round-off (~1e-5), not the double-grade
    // 1e-10 the CPU sibling can assert. Use a float-appropriate tolerance.
    EXPECT_NEAR(frame.meta.get(kEnergy), 0.0, 1e-6);
    EXPECT_NEAR(frame.span(kAtom_fx)[0], 0.0, 1e-3);
    EXPECT_NEAR(frame.span(kAtom_fx)[1], 0.0, 1e-3);
}

TEST(BondHarmonic, DimerDisplaced) {
    // 2-atom system at r != r0 → V > 0
    cpu::Frame frame(2);
    double r = 1.60;
    {
        auto x = frame.span(kAtom_x);
        auto y = frame.span(kAtom_y);
        auto z = frame.span(kAtom_z);
        auto mass = frame.span(kAtom_mass);
        auto tag = frame.span(kAtom_tag);
        x[0] = 0.0;
        x[1] = r;
        y[0] = 0.0;
        y[1] = 0.0;
        z[0] = 0.0;
        z[1] = 0.0;
        mass[0] = 14.0;
        mass[1] = 14.0;
        tag[0] = 0;
        tag[1] = 1;
    }
    frame.box = cpu::Box::cubic(20.0);

    frame.clear_tag_index();
    frame.rebuild_tag_index(0, frame.nlocal());
    frame.add_bond(0, 1, 0);

    double K = 350.0, r0 = 1.53;

    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::BondHarmonic{{K}, {r0}});
    sim.build();
    sim.run(frame);

    double expected_V = K * (r - r0) * (r - r0);
    // Single-precision kernel: ~1e-6 relative round-off on a V≈1.7 energy.
    EXPECT_NEAR(frame.meta.get(kEnergy), expected_V, 1e-4);

    // Forces should be attractive (pulling atoms together)
    EXPECT_GT(frame.span(kAtom_fx)[0], 0.0); // atom 0 pushed right
    EXPECT_LT(frame.span(kAtom_fx)[1], 0.0); // atom 1 pushed left
    // Newton's 3rd law
    EXPECT_NEAR(frame.span(kAtom_fx)[0] + frame.span(kAtom_fx)[1], 0.0, 1e-10);
}

// ── cuda::AngleHarmonic ──────────────────────────────────────────────────

TEST(AngleHarmonic, AtEquilibrium) {
    // 3-atom system at θ0=109.5° → V=0, F=0
    cpu::Frame frame(3);
    double theta0_deg = 109.5;
    double theta0 = theta0_deg * M_PI / 180.0;
    double r = 1.53;

    {
        auto x = frame.span(kAtom_x);
        auto y = frame.span(kAtom_y);
        auto z = frame.span(kAtom_z);
        auto mass = frame.span(kAtom_mass);
        auto tag = frame.span(kAtom_tag);
        x[0] = 0.0;
        x[1] = r;
        x[2] = r + r * cos(M_PI - theta0);
        y[0] = 0.0;
        y[1] = 0.0;
        y[2] = r * sin(M_PI - theta0);
        z[0] = 0.0;
        z[1] = 0.0;
        z[2] = 0.0;
        mass[0] = 14.0;
        mass[1] = 14.0;
        mass[2] = 14.0;
        tag[0] = 0;
        tag[1] = 1;
        tag[2] = 2;
    }
    frame.box = cpu::Box::cubic(20.0);

    frame.clear_tag_index();
    frame.rebuild_tag_index(0, frame.nlocal());
    frame.add_angle(0, 1, 2, 0);

    // AngleHarmonic API contract: theta0 is in RADIANS, matching the CPU
    // sibling (cpu::AngleHarmonic uses theta0_[type] directly in
    // `theta - theta0_[type]`) and the polyethylene_tension benchmark
    // (passes `109.5° · π/180`). Pass `theta0` in radians here.
    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate)
        .add(cuda::AngleHarmonic{{60.0}, {theta0}});
    sim.build();
    sim.run(frame);

    EXPECT_NEAR(frame.meta.get(kEnergy), 0.0, 1e-6);
}

// ── cuda::DihedralMultiHarmonic ──────────────────────────────────────────

TEST(DihedralMultiHarmonic, PerpendicularProducesForce) {
    // φ = 90° dihedral — chosen deliberately to be OFF the V(φ) extrema so
    // BOTH the energy and the torsional force are decisively nonzero. (The
    // earlier "trans" geometry sat at φ=180°, an extremum where dV/dφ=0 and
    // V≈0.006; a force of ~0 is correct there, so an isfinite-only check could
    // not tell "kernel ran" from "kernel never launched" — which is exactly
    // how the __launch_bounds__(128)-vs-256 launch bug hid.)
    //   atom0=(0,1,0) atom1=(0,0,0) atom2=(1.5,0,0) atom3=(1.5,0,1)
    // Planes (0,1,2) and (1,2,3) have normals along -z and -y → φ=90°.
    cpu::Frame frame(4);
    {
        auto x = frame.span(kAtom_x);
        auto y = frame.span(kAtom_y);
        auto z = frame.span(kAtom_z);
        auto mass = frame.span(kAtom_mass);
        auto tag = frame.span(kAtom_tag);
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 1.5;
        x[3] = 1.5;
        y[0] = 1.0;
        y[1] = 0.0;
        y[2] = 0.0;
        y[3] = 0.0;
        z[0] = 0.0;
        z[1] = 0.0;
        z[2] = 0.0;
        z[3] = 1.0;
        mass[0] = 14.0;
        mass[1] = 14.0;
        mass[2] = 14.0;
        mass[3] = 14.0;
        tag[0] = 0;
        tag[1] = 1;
        tag[2] = 2;
        tag[3] = 3;
    }
    frame.box = cpu::Box::cubic(20.0);

    frame.clear_tag_index();
    frame.rebuild_tag_index(0, frame.nlocal());
    frame.add_dihedral(0, 1, 2, 3, 0);

    std::array<double, 5> C = {1.73, -4.49, 0.776, 6.99, 0.0};

    cuda::Evaluate sim;
    sim.on(evaluate::Phase::evaluate).add(cuda::DihedralMultiHarmonic{{C}});
    sim.build();
    sim.run(frame);

    // At φ=90°, cos φ=0 ⇒ V = C0 = 1.73 eV, and dV/dφ = -C1 = 4.49 ≠ 0 ⇒ a
    // real torsional force. Both being nonzero proves the kernel executed.
    EXPECT_TRUE(std::isfinite(frame.meta.get(kEnergy)));
    EXPECT_NEAR(frame.meta.get(kEnergy), 1.73, 0.05)
        << "dihedral energy not ~C0 -> kernel likely never launched";
    double fsum_x = 0.0, fsum_y = 0.0, fsum_z = 0.0, fmax = 0.0;
    auto fx = frame.span(kAtom_fx);
    auto fy = frame.span(kAtom_fy);
    auto fz = frame.span(kAtom_fz);
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(fx[i]));
        EXPECT_TRUE(std::isfinite(fy[i]));
        EXPECT_TRUE(std::isfinite(fz[i]));
        fsum_x += fx[i];
        fsum_y += fy[i];
        fsum_z += fz[i];
        fmax =
            std::max({fmax, std::abs(fx[i]), std::abs(fy[i]), std::abs(fz[i])});
    }
    // A torsional force exists (not all zero) and obeys Newton's 3rd law.
    EXPECT_GT(fmax, 1e-2)
        << "all dihedral forces ~0 -> kernel likely never launched";
    EXPECT_NEAR(fsum_x, 0.0, 1e-3);
    EXPECT_NEAR(fsum_y, 0.0, 1e-3);
    EXPECT_NEAR(fsum_z, 0.0, 1e-3);
}

// ── Integration: bonded + nonbonded ────────────────────────────────

TEST(Bonded, ChainWithLJ) {
    // 4-atom chain with bonds + LJ, 10 steps, finite results
    cpu::Frame frame(4);
    {
        auto x = frame.span(kAtom_x);
        auto y = frame.span(kAtom_y);
        auto z = frame.span(kAtom_z);
        auto mass = frame.span(kAtom_mass);
        auto tag = frame.span(kAtom_tag);
        x[0] = 0.0;
        x[1] = 1.53;
        x[2] = 3.06;
        x[3] = 4.59;
        y[0] = 0.0;
        y[1] = 0.3;
        y[2] = -0.3;
        y[3] = 0.3;
        z[0] = 0.0;
        z[1] = 0.0;
        z[2] = 0.0;
        z[3] = 0.0;
        mass[0] = 14.0;
        mass[1] = 14.0;
        mass[2] = 14.0;
        mass[3] = 14.0;
        tag[0] = 0;
        tag[1] = 1;
        tag[2] = 2;
        tag[3] = 3;
    }
    frame.box = cpu::Box::cubic(20.0);

    frame.clear_tag_index();
    frame.rebuild_tag_index(0, frame.nlocal());
    frame.add_bond(0, 1, 0);
    frame.add_bond(1, 2, 0);
    frame.add_bond(2, 3, 0);

    cuda::md::MD sim;
    sim.set_integrator(cuda::md::VelocityVerlet{0.1});
    sim.on(md::Phase::force).add(cuda::BondHarmonic{{350.0}, {1.53}});
    sim.build();
    sim.run(frame, 10);

    auto x = frame.span(kAtom_x);
    auto vx = frame.span(kAtom_vx);
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(x[i]));
        EXPECT_TRUE(std::isfinite(vx[i]));
    }
}
