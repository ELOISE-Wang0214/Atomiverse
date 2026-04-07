/// @file tests/cpu/md/test_md.cpp
/// @brief Unit and integration tests for CPU MD backend components.

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

// ── Fixture: FCC Argon system ───────────────────────────────────────

class CpuMdTest : public ::testing::Test {
protected:
    static cpu::Frame make_fcc_argon(int n_cells, double lattice_const) {
        int atoms_per_cell = 4;
        int n_atoms = n_cells * n_cells * n_cells * atoms_per_cell;
        cpu::Frame frame(n_atoms);

        double basis[4][3] = {
            {0.0, 0.0, 0.0},
            {0.5, 0.5, 0.0},
            {0.5, 0.0, 0.5},
            {0.0, 0.5, 0.5}
        };

        int idx = 0;
        for (int ix = 0; ix < n_cells; ++ix)
            for (int iy = 0; iy < n_cells; ++iy)
                for (int iz = 0; iz < n_cells; ++iz)
                    for (int b = 0; b < 4; ++b) {
                        frame.atoms.x[idx] = (ix + basis[b][0]) * lattice_const;
                        frame.atoms.y[idx] = (iy + basis[b][1]) * lattice_const;
                        frame.atoms.z[idx] = (iz + basis[b][2]) * lattice_const;
                        frame.atoms.mass[idx] = double(39.948);
                        frame.atoms.type[idx] = 18;
                        ++idx;
                    }

        double box_len = n_cells * lattice_const;
        frame.box = Mat3::cubic(box_len);
        return frame;
    }

    // LJ parameters for Argon
    static constexpr double kEpsilon = double(0.0104);
    static constexpr double kSigma = double(3.4);
    static constexpr double kCutoff = double(6.8);
};

// ── ZeroForces ──────────────────────────────────────────────────────

TEST_F(CpuMdTest, ZeroForces) {
    cpu::Frame frame(10);
    for (int i = 0; i < 10; ++i) {
        frame.atoms.fx[i] = double(1);
        frame.atoms.fy[i] = double(2);
        frame.atoms.fz[i] = double(3);
    }
    frame.atoms.mass = xt::ones<double>({int(10)});
    frame.box = Mat3::cubic(10.0);

    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    cpu::md::ZeroForces zf;
    zf.launch(lctx);

    for (int i = 0; i < 10; ++i) {
        EXPECT_EQ(frame.atoms.fx[i], 0);
        EXPECT_EQ(frame.atoms.fy[i], 0);
        EXPECT_EQ(frame.atoms.fz[i], 0);
    }
}

// ── KineticEnergy ───────────────────────────────────────────────────

TEST_F(CpuMdTest, KineticEnergy) {
    cpu::Frame frame(2);
    frame.atoms.mass[0] = double(1);  frame.atoms.mass[1] = double(1);
    frame.atoms.vx[0] = double(1);   frame.atoms.vy[0] = double(0); frame.atoms.vz[0] = double(0);
    frame.atoms.vx[1] = double(0);   frame.atoms.vy[1] = double(1); frame.atoms.vz[1] = double(0);
    frame.box = Mat3::cubic(10.0);
    frame.kinetic = 0;

    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    cpu::md::KineticEnergy ke;
    ke.launch(lctx);

    EXPECT_NEAR(static_cast<double>(frame.kinetic),
                static_cast<double>(kMvv2e), 1e-5);
}

// ── LJ force ────────────────────────────────────────────────────────

TEST_F(CpuMdTest, LjAnalyticForce) {
    cpu::Frame frame(2);
    double sigma = double(3.4);
    double epsilon = double(0.0104);
    double cutoff = double(8.0);

    double r = double(4.0);
    frame.atoms.x[0] = double(0); frame.atoms.y[0] = double(0); frame.atoms.z[0] = double(0);
    frame.atoms.x[1] = r;       frame.atoms.y[1] = double(0); frame.atoms.z[1] = double(0);
    frame.atoms.mass[0] = double(39.948); frame.atoms.mass[1] = double(39.948);
    frame.box = Mat3::cubic(20.0);
    frame.energy = 0;

    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    cpu::LaunchContext lctx{ctx, frame, 0, 0};

    cpu::md::NeighborList<> nlist(cutoff);
    nlist.build(bctx);
    nlist.launch(lctx);

    cpu::md::LJ lj(epsilon, sigma, cutoff);
    lj.build(bctx);
    lj.launch(lctx);

    double sig_r = static_cast<double>(sigma) / static_cast<double>(r);
    double sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
    double sig_r12 = sig_r6 * sig_r6;
    double f_analytic = 24.0 * static_cast<double>(epsilon) / static_cast<double>(r)
                        * (2.0 * sig_r12 - sig_r6);

    EXPECT_NEAR(std::abs(static_cast<double>(frame.atoms.fx[0])),
                std::abs(f_analytic), 1e-8);
    EXPECT_NEAR(static_cast<double>(frame.atoms.fx[0]),
                -static_cast<double>(frame.atoms.fx[1]), 1e-10);

    lj.destroy({ctx, frame});
}

TEST_F(CpuMdTest, ForceSumZero) {
    auto frame = make_fcc_argon(3, double(5.26));

    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};

    int n = frame.atoms.size();
    for (int i = 0; i < n; ++i) {
        frame.atoms.fx[i] = 0; frame.atoms.fy[i] = 0; frame.atoms.fz[i] = 0;
    }

    cpu::md::NeighborList<> nlist(kCutoff);
    nlist.build(bctx);
    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    nlist.launch(lctx);

    cpu::md::LJ lj(kEpsilon, kSigma, kCutoff);
    lj.build(bctx);
    lj.launch(lctx);

    double sum_fx = 0, sum_fy = 0, sum_fz = 0;
    for (int i = 0; i < n; ++i) {
        sum_fx += frame.atoms.fx[i];
        sum_fy += frame.atoms.fy[i];
        sum_fz += frame.atoms.fz[i];
    }

    EXPECT_NEAR(static_cast<double>(sum_fx), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(sum_fy), 0.0, 1e-6);
    EXPECT_NEAR(static_cast<double>(sum_fz), 0.0, 1e-6);

    lj.destroy({ctx, frame});
}

// ── CellList ────────────────────────────────────────────────────────

TEST_F(CpuMdTest, CellListPairCount) {
    auto frame = make_fcc_argon(3, double(5.26));

    cpu::md::CellList cell_list(kCutoff);
    cell_list.init(frame, kCutoff);
    cell_list.build(frame);
    auto pairs = cell_list.build_half(frame);

    int n = frame.atoms.size();
    double skin = double(0.3) * kCutoff;
    double rc_skin_sq = (kCutoff + skin) * (kCutoff + skin);
    int bf_count = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = frame.atoms.x[i] - frame.atoms.x[j];
            double dy = frame.atoms.y[i] - frame.atoms.y[j];
            double dz = frame.atoms.z[i] - frame.atoms.z[j];
            double lx = frame.box(0, 0), ly = frame.box(1, 1), lz = frame.box(2, 2);
            if (dx >  lx * double(0.5)) dx -= lx;
            if (dx < -lx * double(0.5)) dx += lx;
            if (dy >  ly * double(0.5)) dy -= ly;
            if (dy < -ly * double(0.5)) dy += ly;
            if (dz >  lz * double(0.5)) dz -= lz;
            if (dz < -lz * double(0.5)) dz += lz;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc_skin_sq) ++bf_count;
        }
    }

    EXPECT_EQ(pairs.total_pairs, bf_count);
}

// ── NVE energy conservation ─────────────────────────────────────────

TEST_F(CpuMdTest, NveEnergyConservation) {
    auto frame = make_fcc_argon(3, double(5.26));
    double dt = double(1.0);

    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.add(cpu::md::InitVelocity{double(50.0), 42});
        sim.on(ATV::md::Phase::force)
           .add(cpu::md::NeighborList<>{kCutoff}, -10)
           .add(cpu::md::LJ{kEpsilon, kSigma, kCutoff});
        sim.build();
        sim.run(frame, 100);
    }

    double e_first = frame.energy + frame.kinetic;

    {
        cpu::md::MD sim;
        sim.add(cpu::md::VelocityVerlet{dt});
        sim.on(ATV::md::Phase::force)
           .add(cpu::md::NeighborList<>{kCutoff}, -10)
           .add(cpu::md::LJ{kEpsilon, kSigma, kCutoff});
        sim.build();
        sim.run(frame, 100);
    }

    double e_second = frame.energy + frame.kinetic;

    double e1 = static_cast<double>(e_first);
    double e2 = static_cast<double>(e_second);
    ASSERT_FALSE(std::isnan(e1)) << "First energy is NaN";
    ASSERT_FALSE(std::isnan(e2)) << "Second energy is NaN";
    ASSERT_NE(e1, 0.0) << "First energy is zero";

    double relative_drift = std::abs(e2 - e1) / std::abs(e1);
    EXPECT_LT(relative_drift, 0.05);
}

// ── Bonded forces ───────────────────────────────────────────────────

TEST_F(CpuMdTest, BondedForce) {
    cpu::Frame frame(2);
    frame.atoms.x[0] = double(0); frame.atoms.y[0] = double(0); frame.atoms.z[0] = double(0);
    frame.atoms.x[1] = double(1.5); frame.atoms.y[1] = double(0); frame.atoms.z[1] = double(0);
    frame.atoms.mass[0] = double(12); frame.atoms.mass[1] = double(12);
    frame.box = Mat3::cubic(20.0);

    frame.bonds.atomi = xt::xtensor<int, 1>({int(1)});
    frame.bonds.atomj = xt::xtensor<int, 1>({int(1)});
    frame.bonds.type = xt::xtensor<int, 1>({int(1)});
    frame.bonds.atomi[0] = 0;
    frame.bonds.atomj[0] = 1;
    frame.bonds.type[0] = 0;

    double K = double(100.0);
    double r0 = double(1.0);

    cpu::bonded::HarmonicBond bond({K}, {r0});
    cpu::Context ctx;
    cpu::BuildContext bctx{ctx, frame};
    bond.build(bctx);

    cpu::LaunchContext lctx{ctx, frame, 0, 0};
    bond.launch(lctx);

    double fx0_expected = 100.0;

    EXPECT_NEAR(static_cast<double>(frame.atoms.fx[0]), fx0_expected, 1e-6);
    EXPECT_NEAR(static_cast<double>(frame.atoms.fx[0]),
                -static_cast<double>(frame.atoms.fx[1]), 1e-10);

    bond.destroy({ctx, frame});
}
