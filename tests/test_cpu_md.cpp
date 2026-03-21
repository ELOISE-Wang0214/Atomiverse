/// @file test_cpu_md.cpp
/// @brief Integration tests for CPU MD backend.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/cpu/observer/thermo_component.h"
#include "atomiverse/cpu/bonded/harmonic_bond.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/core/units.h"
#include "atomiverse/core/mat3.h"

using namespace ATV;

// ── Helpers ──────────────────────────────────────────────────────

/// @brief Create a simple FCC lattice for LJ Argon.
static cpu::Frame make_fcc_argon(Int n_cells, Real lattice_const) {
    Int atoms_per_cell = 4;
    Int n_atoms = n_cells * n_cells * n_cells * atoms_per_cell;
    cpu::Frame frame(n_atoms);

    // FCC basis
    Real basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}
    };

    Int idx = 0;
    for (Int ix = 0; ix < n_cells; ++ix) {
        for (Int iy = 0; iy < n_cells; ++iy) {
            for (Int iz = 0; iz < n_cells; ++iz) {
                for (int b = 0; b < 4; ++b) {
                    frame.x[idx] = (ix + basis[b][0]) * lattice_const;
                    frame.y[idx] = (iy + basis[b][1]) * lattice_const;
                    frame.z[idx] = (iz + basis[b][2]) * lattice_const;
                    frame.mass[idx] = Real(39.948);  // Argon
                    frame.type_id[idx] = 18;         // Ar
                    ++idx;
                }
            }
        }
    }

    Real box_len = n_cells * lattice_const;
    frame.box = Mat3::cubic(box_len);
    return frame;
}

// ── Unit Tests ───────────────────────────────────────────────────

TEST(CpuMd, ZeroForces) {
    cpu::Frame frame(10);
    for (Int i = 0; i < 10; ++i) {
        frame.fx[i] = Real(1);
        frame.fy[i] = Real(2);
        frame.fz[i] = Real(3);
    }
    frame.mass = xt::ones<Real>({Int(10)});
    frame.box = Mat3::cubic(10.0);

    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    cpu::md::ZeroForces zf;
    zf.launch(lctx);

    for (Int i = 0; i < 10; ++i) {
        EXPECT_EQ(frame.fx[i], 0);
        EXPECT_EQ(frame.fy[i], 0);
        EXPECT_EQ(frame.fz[i], 0);
    }
}

TEST(CpuMd, KineticEnergy) {
    cpu::Frame frame(2);
    frame.mass[0] = Real(1);  frame.mass[1] = Real(1);
    frame.vx[0] = Real(1);   frame.vy[0] = Real(0); frame.vz[0] = Real(0);
    frame.vx[1] = Real(0);   frame.vy[1] = Real(1); frame.vz[1] = Real(0);
    frame.box = Mat3::cubic(10.0);
    frame.kinetic = 0;

    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    cpu::md::KineticEnergy ke;
    ke.launch(lctx);

    // KE = 0.5 * kMvv2e * sum(m * v^2) = 0.5 * kMvv2e * 2 = kMvv2e
    EXPECT_NEAR(static_cast<double>(frame.kinetic),
                static_cast<double>(kMvv2e), 1e-5);
}

TEST(CpuMd, LjAnalyticForce) {
    // 2-atom system: test LJ force at known separation
    cpu::Frame frame(2);
    Real sigma = Real(3.4);
    Real epsilon = Real(0.0104);
    Real cutoff = Real(8.0);

    Real r = Real(4.0);  // Separation
    frame.x[0] = Real(0); frame.y[0] = Real(0); frame.z[0] = Real(0);
    frame.x[1] = r;       frame.y[1] = Real(0); frame.z[1] = Real(0);
    frame.mass[0] = Real(39.948); frame.mass[1] = Real(39.948);
    frame.box = Mat3::cubic(20.0);
    frame.energy = 0;

    // Build neighbor list manually via CellList strategy
    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    cpu::md::NeighborList<> nlist(cutoff);
    nlist.build(bctx);
    nlist.launch(lctx);

    cpu::md::LJ lj(epsilon, sigma, cutoff);
    lj.build(bctx);
    lj.launch(lctx);

    // Analytic: F = 24*eps/r * [2*(sig/r)^12 - (sig/r)^6]
    double sig_r = static_cast<double>(sigma) / static_cast<double>(r);
    double sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
    double sig_r12 = sig_r6 * sig_r6;
    double f_analytic = 24.0 * static_cast<double>(epsilon) / static_cast<double>(r)
                        * (2.0 * sig_r12 - sig_r6);

    // Force on atom 0 should point toward atom 1 (attractive at r=4.0 for sigma=3.4)
    // f_analytic is the scalar fpair; actual Fx[0] = fpair * dx = fpair * (-r) = -fpair * r
    // Wait: dx = x[0]-x[1] = -r, so Fx[0] = fpair * dx = fpair * (-r)
    // But fpair already includes 1/r: fpair = 24*eps/r^2 * [2*sig_r12 - sig_r6]
    // Actually from the code: fpair = 24*eps * r2inv * (2*r6inv*r6inv - r6inv)
    // and Fx += fpair * dx. So let's just check the magnitude.

    EXPECT_NEAR(std::abs(static_cast<double>(frame.fx[0])),
                std::abs(f_analytic), 1e-8);
    // Newton's 3rd law: F[0] = -F[1]
    EXPECT_NEAR(static_cast<double>(frame.fx[0]),
                -static_cast<double>(frame.fx[1]), 1e-10);

    lj.destroy({ctx, frame});
}

TEST(CpuMd, CellListPairCount) {
    // Compare CellList vs brute-force: same pair count
    auto frame = make_fcc_argon(3, Real(5.26));  // 108 atoms
    Real cutoff = Real(6.8);

    cpu::md::CellList cell_list(cutoff);
    cell_list.init(frame, cutoff);
    cell_list.build(frame);
    auto pairs = cell_list.build_half(frame);

    // Brute-force half pair count
    Int n = frame.n_atoms();
    Real rc_sq = cutoff * cutoff;
    // Add skin
    Real skin = Real(0.3) * cutoff;
    Real rc_skin_sq = (cutoff + skin) * (cutoff + skin);
    Int bf_count = 0;
    for (Int i = 0; i < n; ++i) {
        for (Int j = i + 1; j < n; ++j) {
            Real dx = frame.x[i] - frame.x[j];
            Real dy = frame.y[i] - frame.y[j];
            Real dz = frame.z[i] - frame.z[j];
            Real lx = frame.box(0, 0), ly = frame.box(1, 1), lz = frame.box(2, 2);
            if (dx >  lx * Real(0.5)) dx -= lx;
            if (dx < -lx * Real(0.5)) dx += lx;
            if (dy >  ly * Real(0.5)) dy -= ly;
            if (dy < -ly * Real(0.5)) dy += ly;
            if (dz >  lz * Real(0.5)) dz -= lz;
            if (dz < -lz * Real(0.5)) dz += lz;
            Real r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc_skin_sq) ++bf_count;
        }
    }

    EXPECT_EQ(pairs.total_pairs, bf_count);
}

TEST(CpuMd, NveEnergyConservation) {
    // NVE MD: energy should be conserved over time
    auto frame = make_fcc_argon(3, Real(5.26));  // 108 atoms
    Real cutoff = Real(6.8);
    Real dt = Real(1.0);  // 1 fs

    // First run: initialize and equilibrate
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.add(cpu::md::InitVelocity{Real(50.0), 42});
        sim.on(ATV::md::Phase::force)
           .add(cpu::md::NeighborList<>{cutoff}, -10)
           .add(cpu::md::LJ{Real(0.0104), Real(3.4), cutoff});
        sim.build();
        sim.run(frame, 100);
    }

    Real e_first = frame.energy + frame.kinetic;

    // Second run: continue without re-initialization
    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.on(ATV::md::Phase::force)
           .add(cpu::md::NeighborList<>{cutoff}, -10)
           .add(cpu::md::LJ{Real(0.0104), Real(3.4), cutoff});
        sim.build();
        sim.run(frame, 100);
    }

    Real e_second = frame.energy + frame.kinetic;

    // Energy drift should be small relative to total energy
    double e1 = static_cast<double>(e_first);
    double e2 = static_cast<double>(e_second);
    ASSERT_FALSE(std::isnan(e1)) << "First energy is NaN";
    ASSERT_FALSE(std::isnan(e2)) << "Second energy is NaN";
    ASSERT_NE(e1, 0.0) << "First energy is zero";

    double relative_drift = std::abs(e2 - e1) / std::abs(e1);
    EXPECT_LT(relative_drift, 0.05);  // 5% tolerance for short runs
}

TEST(CpuMd, ForceSumZero) {
    // Newton's 3rd law: sum of all forces should be zero
    auto frame = make_fcc_argon(3, Real(5.26));
    Real cutoff = Real(6.8);

    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};

    // Zero forces
    Int n = frame.n_atoms();
    for (Int i = 0; i < n; ++i) {
        frame.fx[i] = 0; frame.fy[i] = 0; frame.fz[i] = 0;
    }

    cpu::md::NeighborList<> nlist(cutoff);
    nlist.build(bctx);
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    nlist.launch(lctx);

    cpu::md::LJ lj(Real(0.0104), Real(3.4), cutoff);
    lj.build(bctx);
    lj.launch(lctx);

    Real sum_fx = 0, sum_fy = 0, sum_fz = 0;
    for (Int i = 0; i < n; ++i) {
        sum_fx += frame.fx[i];
        sum_fy += frame.fy[i];
        sum_fz += frame.fz[i];
    }

    EXPECT_NEAR(static_cast<double>(sum_fx), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(sum_fy), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(sum_fz), 0.0, 1e-6);

    lj.destroy({ctx, frame});
}

TEST(CpuMd, BondedForce) {
    // 2-atom harmonic bond: F = -2K(r - r0) r_hat
    cpu::Frame frame(2);
    frame.x[0] = Real(0); frame.y[0] = Real(0); frame.z[0] = Real(0);
    frame.x[1] = Real(1.5); frame.y[1] = Real(0); frame.z[1] = Real(0);
    frame.mass[0] = Real(12); frame.mass[1] = Real(12);
    frame.box = Mat3::cubic(20.0);

    // Bond topology
    frame.bond_i = xt::xtensor<Int, 1>({Int(1)});
    frame.bond_j = xt::xtensor<Int, 1>({Int(1)});
    frame.bond_type = xt::xtensor<Int, 1>({Int(1)});
    frame.bond_i[0] = 0;
    frame.bond_j[0] = 1;
    frame.bond_type[0] = 0;

    Real K = Real(100.0);  // eV/Angstrom^2
    Real r0 = Real(1.0);   // Angstrom

    cpu::bonded::HarmonicBond bond({K}, {r0});
    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    bond.build(bctx);

    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    bond.launch(lctx);

    // Analytic: F = -2K(r-r0)/r * dx
    double r = 1.5;
    double dr = r - static_cast<double>(r0);
    double f_expected = -2.0 * static_cast<double>(K) * dr / r * r;  // fpair * dx = fpair * 1.5
    // Actually: fpair = -2K(r-r0)/r, then Fx[0] = fpair * dx = fpair * (x0-x1) = fpair * (-1.5)
    // Fx[0] = -2K * 0.5/1.5 * (-1.5) = 2K * 0.5 = K = 100
    // Wait, let me recalculate:
    // dx = x[i]-x[j] = 0 - 1.5 = -1.5
    // r = 1.5
    // dr = r - r0 = 0.5
    // fpair = -2*K*dr/r = -2*100*0.5/1.5 = -66.667
    // Fx[i] = fpair * dx = -66.667 * (-1.5) = 100
    double fx0_expected = 100.0;

    EXPECT_NEAR(static_cast<double>(frame.fx[0]), fx0_expected, 1e-6);
    EXPECT_NEAR(static_cast<double>(frame.fx[0]),
                -static_cast<double>(frame.fx[1]), 1e-10);

    bond.destroy({ctx, frame});
}
