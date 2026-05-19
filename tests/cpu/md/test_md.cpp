/// @file tests/cpu/md/test_md.cpp
/// @brief Unit and integration tests for CPU MD backend components.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/cpu/md/fix_langevin.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/md/velocity_verlet.h"
#include "atomiverse/cpu/md/lj.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/nose_hoover.h"
#include "atomiverse/cpu/md/nose_hoover_npt.h"
#include "atomiverse/cpu/md/nve_limit.h"
#include "atomiverse/cpu/md/neighbor_list.h"
#include "atomiverse/cpu/md/pair_list.h"
#include "atomiverse/cpu/md/virial_stress.h"
#include "atomiverse/cpu/recorder/thermo.h"
#include "atomiverse/cpu/potentials/harmonic_bond.h"
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

/// Reproducer: when `update()` is called with a smaller `n_active` than
/// the previous `build()`, the stale `binned_atoms_` entries from rows
/// that are no longer "active" must NOT leak into the new pair list --
/// otherwise `build_half` will index `frame.atoms.x[]` past the active
/// count and crash. This is the failure mode triggered by AtomMigration
/// shrinking `dd->n_local + dd->n_ghost` between MD steps under MPI.
TEST_F(CpuMdTest, CellListShrinkingNActive) {
    // 4x4x4 FCC = 256 atoms in 24 Å box.
    auto frame = make_fcc_argon(4, /*lattice_const=*/6.0);
    const int n_total = frame.atoms.size();
    ASSERT_EQ(n_total, 256);

    cpu::md::CellList cell_list(kCutoff);
    cell_list.init(frame, kCutoff);

    // First build with all atoms active.
    cell_list.build(frame, n_total);
    auto pairs_full = cell_list.build_half(frame, n_total);

    // Now shrink: pretend half the atoms migrated away. Subsequent
    // `update()` should NOT report any pair index beyond [0, n_active).
    const int n_shrink = n_total / 2;
    // Force movement so needs_rebuild based on displacement *would* still
    // skip the rebuild — only the n change should trigger it.
    cell_list.update(frame, n_shrink);
    auto pairs_shrink = cell_list.build_half(frame, n_shrink);

    EXPECT_LE(pairs_shrink.inum, n_shrink) << "inum exceeds n_active";

    for (int idx = 0; idx < pairs_shrink.total_pairs; ++idx) {
        const int j = pairs_shrink.neighbors[idx];
        ASSERT_GE(j, 0)        << "negative neighbor index";
        ASSERT_LT(j, n_shrink) << "neighbor index " << j
                                << " >= n_active=" << n_shrink
                                << " (stale binned_atoms_ leaked)";
    }
}

/// Reproducer #2: NPT-style affine box-and-atom rescale. Pair count must
/// be invariant under affine scaling (atoms preserve relative distances
/// scaled by `scale`, so the same pairs are within `cutoff*scale` of each
/// other). If the cell grid is sized for the OLD box, atoms that drift
/// into cells that no longer span `cutoff` will produce false NEGATIVES
/// and miss pairs -- the eventual symptom is an LJ blow-up when one such
/// "missed" close pair finally registers.
TEST_F(CpuMdTest, CellListExtendedGridBoxChange) {
    auto frame = make_fcc_argon(4, /*a=*/6.0);
    const int n_total = frame.atoms.size();

    cpu::md::CellList cl(kCutoff);
    int dims[3] = {2, 1, 1};
    cl.init_extended(frame, kCutoff, dims, /*capacity=*/n_total);

    cl.build(frame, n_total);
    auto pairs0 = cl.build_half(frame, n_total);
    ASSERT_GT(pairs0.total_pairs, 0);
    const int pairs_before = pairs0.total_pairs;

    // NPT-like: expand the box by 2% and displace atoms affinely.
    const double scale = 1.02;
    frame.box.scale_isotropic(scale);
    for (int i = 0; i < n_total; ++i) {
        frame.atoms.x[i] *= scale;
        frame.atoms.y[i] *= scale;
        frame.atoms.z[i] *= scale;
    }

    cl.update(frame, n_total);
    auto pairs1 = cl.build_half(frame, n_total);

    EXPECT_LE(pairs1.inum, n_total);
    for (int idx = 0; idx < pairs1.total_pairs; ++idx) {
        const int j = pairs1.neighbors[idx];
        ASSERT_GE(j, 0);
        ASSERT_LT(j, n_total) << "neighbor index OOB after box change";
    }
    // Affine rescale by `scale` < 1 + skin/cutoff preserves the pair set:
    // the smallest pair distance was at most `cutoff_skin`; after the 2%
    // expansion (with skin = 0.3*cutoff), the same pair is still within
    // 1.02 * cutoff_skin < cutoff_skin * (1 + skin_fraction) so it must
    // still be picked. Allow a small downward drift (~5%) but reject the
    // catastrophic "missed everything" failure mode.
    EXPECT_GE(pairs1.total_pairs, pairs_before * 9 / 10)
        << "stale cell grid missed pairs after box rescale: "
        << "before=" << pairs_before << " after=" << pairs1.total_pairs;
}

/// Reproducer #3: full multi-rank sequence -- shrink n_active (migration
/// out), then grow it (ghosts arriving), then shrink again. This is the
/// exact pattern NoseHooverNPT provokes each step.
TEST_F(CpuMdTest, CellListMultiRankMigrationCycle) {
    auto frame = make_fcc_argon(4, /*a=*/6.0);
    const int N = frame.atoms.size();

    cpu::md::CellList cl(kCutoff);
    int dims[3] = {2, 1, 1};
    cl.init_extended(frame, kCutoff, dims, /*capacity=*/N);

    const int n_owners = N / 2;
    const int n_with_ghost = n_owners + N / 4;

    // Launch 1: owners only.
    cl.build(frame, n_owners);
    auto p1 = cl.build_half(frame, n_owners);
    for (int k = 0; k < p1.total_pairs; ++k) {
        ASSERT_LT(p1.neighbors[k], n_owners) << "launch 1 OOB";
    }

    // Launch 2: owners + ghosts (grew). Mimic halo exchange by writing
    // new atom positions into slots [n_owners, n_with_ghost) -- the real
    // framework's HaloExchangeForward does this.
    for (int i = n_owners; i < n_with_ghost; ++i) {
        frame.atoms.x[i] = frame.atoms.x[i - n_owners] - frame.box(0,0);
        frame.atoms.y[i] = frame.atoms.y[i - n_owners];
        frame.atoms.z[i] = frame.atoms.z[i - n_owners];
    }
    cl.update(frame, n_with_ghost);
    auto p2 = cl.build_half(frame, n_with_ghost);
    for (int k = 0; k < p2.total_pairs; ++k) {
        ASSERT_LT(p2.neighbors[k], n_with_ghost) << "launch 2 OOB";
    }

    // Launch 3: migration shrank owners again. The previously-written
    // ghost slots now carry stale data; if the rebuild guard doesn't
    // fire, binned_atoms_ can still hand out those stale indices.
    const int n_shrunk = n_owners - 10;
    cl.update(frame, n_shrunk);
    auto p3 = cl.build_half(frame, n_shrunk);
    for (int k = 0; k < p3.total_pairs; ++k) {
        ASSERT_LT(p3.neighbors[k], n_shrunk) << "launch 3 OOB";
    }
}

/// PBC edge pair: with cutoff=6.0 and a box of 30 Å, two atoms placed
/// at x=0.5 and x=29.5 are 1 Å apart through the periodic image. With
/// `init` (single-rank, all axes wrap during binning), these two atoms
/// should still be reported as one pair via the minimum-image distance
/// check. The half-stencil only spans cell offsets ±1 per axis, so this
/// also implicitly checks that the cell list either uses periodic
/// stencil offsets, has cells big enough that two opposite-edge atoms
/// land in adjacent cells, or relies on something else for PBC.
TEST_F(CpuMdTest, CellListPbcPairAcrossEdge) {
    cpu::Frame frame(2);
    frame.atoms.x[0] = 0.5;  frame.atoms.y[0] = 1.0; frame.atoms.z[0] = 1.0;
    frame.atoms.x[1] = 29.5; frame.atoms.y[1] = 1.0; frame.atoms.z[1] = 1.0;
    for (int i = 0; i < 2; ++i) {
        frame.atoms.mass[i] = 1.0;
    }
    frame.box = Mat3::cubic(30.0);

    cpu::md::CellList cl(/*cutoff=*/6.0);
    cl.init(frame, 6.0);
    cl.build(frame, 2);
    auto pairs = cl.build_half(frame, 2);
    EXPECT_EQ(pairs.total_pairs, 1)
        << "PBC pair across box edge not detected by cell list "
        << "(cell list misses pairs at periodic boundaries)";
}

/// Reproducer #4: the bin_atoms scatter phase assumes `binned_atoms_`
/// slots beyond the current scatter are never READ. But if n grew
/// between launches (no rebuild required by displacement), `update`
/// skips rebuild, so any read into `binned_atoms_` could return stale
/// larger-index data from the previous build. Verify that an n_active
/// that GROWS between calls still yields valid neighbors.
TEST_F(CpuMdTest, CellListGrowingNActive) {
    auto frame = make_fcc_argon(4, /*a=*/6.0);
    const int N = frame.atoms.size();

    cpu::md::CellList cl(kCutoff);
    cl.init(frame, kCutoff);

    const int n_small = N / 2;
    cl.build(frame, n_small);
    auto p_small = cl.build_half(frame, n_small);
    for (int k = 0; k < p_small.total_pairs; ++k) {
        ASSERT_LT(p_small.neighbors[k], n_small);
    }

    // Grow to full set without any displacement.
    cl.update(frame, N);
    auto p_full = cl.build_half(frame, N);
    EXPECT_LE(p_full.inum, N);
    for (int k = 0; k < p_full.total_pairs; ++k) {
        ASSERT_LT(p_full.neighbors[k], N);
    }
    // A grown launch must actually see more atoms than the half launch.
    EXPECT_GT(p_full.total_pairs, p_small.total_pairs)
        << "growing n_active produced fewer pairs -- cell list must not "
           "have rebuilt";
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

    HarmonicBond bond({K}, {r0});
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

// ── NVELimit ────────────────────────────────────────────────────────

TEST_F(CpuMdTest, NveLimitCapsHotVelocity) {
    cpu::Frame frame(2);
    frame.atoms.x[0] = 0; frame.atoms.y[0] = 0; frame.atoms.z[0] = 0;
    frame.atoms.x[1] = 5; frame.atoms.y[1] = 0; frame.atoms.z[1] = 0;
    frame.atoms.mass[0] = 1.0; frame.atoms.mass[1] = 1.0;
    // pre-existing huge velocity in +x: |v| = 10 Å/fs
    frame.atoms.vx[0] = 10.0;
    frame.box = Mat3::cubic(20.0);

    const double dt = 0.5;
    const double xmax = 0.05;     // → vmax = 0.1 Å/fs
    cpu::md::NVELimit limit(dt, xmax);

    cpu::Context ctx;
    cpu::LaunchContext lctx_init{ctx, frame, static_cast<int>(ATV::md::Phase::initial_integrate), 0};
    limit.launch(lctx_init);

    const double vmax = xmax / dt;
    const double v_after = std::sqrt(frame.atoms.vx[0] * frame.atoms.vx[0]
                                   + frame.atoms.vy[0] * frame.atoms.vy[0]
                                   + frame.atoms.vz[0] * frame.atoms.vz[0]);
    EXPECT_NEAR(v_after, vmax, 1e-10) << "NVELimit failed to cap |v| <= xmax/dt";
}

TEST_F(CpuMdTest, NveLimitNoCapBelowThreshold) {
    cpu::Frame frame(1);
    frame.atoms.mass[0] = 1.0;
    frame.atoms.vx[0] = 0.05;     // |v| = 0.05, vmax = 0.1 -> no cap
    frame.box = Mat3::cubic(20.0);

    cpu::md::NVELimit limit(0.5, 0.05);
    cpu::Context ctx;
    cpu::LaunchContext lctx{ctx, frame, static_cast<int>(ATV::md::Phase::final_integrate), 0};
    limit.launch(lctx);
    EXPECT_NEAR(static_cast<double>(frame.atoms.vx[0]), 0.05, 1e-12);
}

// ── FixLangevin ─────────────────────────────────────────────────────
//
// Statistical correctness: under fix_langevin in NVE thermostatted by the
// drag+random forces, the equilibrium <KE> should approach (3N/2)*kT for
// large N, matching the canonical ensemble. We seed deterministically so
// the test is reproducible.

TEST_F(CpuMdTest, FixLangevinReachesEquipartition) {
    const int n_atoms = 256;
    const double T_target = 300.0;
    const double damp     = 50.0;        // fs
    const double dt       = 1.0;         // fs
    const double mass     = 16.0;        // amu

    auto frame = make_fcc_argon(4, 4.0); // 256 atoms, neutral box
    for (int i = 0; i < n_atoms; ++i) frame.atoms.mass[i] = mass;

    cpu::md::MD sim;
    sim.add(cpu::md::VelocityVerlet{dt});
    sim.add(cpu::md::FixLangevin{T_target, damp, dt, /*seed=*/12345});
    sim.build();

    sim.run(frame, 4000);  // ~80 damping times -- well-equilibrated

    // Average KE from current velocities (no time integration, just sample).
    double KE = 0;
    for (int i = 0; i < n_atoms; ++i) {
        const double v2 = frame.atoms.vx[i]*frame.atoms.vx[i]
                        + frame.atoms.vy[i]*frame.atoms.vy[i]
                        + frame.atoms.vz[i]*frame.atoms.vz[i];
        KE += 0.5 * mass * v2;
    }
    KE *= kMvv2e;

    const double T_inst = (2.0 * KE) / (3.0 * (n_atoms - 1) * kBoltzmann);
    // Single-snapshot temperature has stochastic spread; 256 atoms gives
    // a relative width of ~sqrt(2/(3N)) ≈ 5%, so allow ±15% for safety.
    EXPECT_NEAR(T_inst, T_target, 0.15 * T_target)
        << "fix_langevin failed to thermostat to within 15% of target";
}

// ── NoseHooverNPT ───────────────────────────────────────────────────

TEST_F(CpuMdTest, NoseHooverNptRelaxesPressure) {
    // Slightly compressed FCC argon -- LJ pair virial drives P > 0.
    auto frame = make_fcc_argon(4, 5.1);  // a < equilibrium for argon
    const int n = frame.atoms.size();
    for (int i = 0; i < n; ++i) frame.atoms.mass[i] = 39.948;

    const double dt    = 1.0;
    const double T_tgt = 100.0;
    const double tau_T = 100.0;
    const double tau_P = 1000.0;

    cpu::md::MD sim;
    sim.add(cpu::md::InitVelocity{T_tgt, /*seed=*/777});
    sim.add(cpu::md::NoseHooverNPT{{T_tgt, tau_T, /*P_target=*/0.0, tau_P, dt}});
    sim.on(ATV::md::Phase::force).add(cpu::md::NeighborList<>{kCutoff}, -10);
    sim.on(ATV::md::Phase::force).add(cpu::md::LJ{kEpsilon, kSigma, kCutoff});
    sim.add(cpu::md::VirialStress{});
    sim.build();

    const double V0 = frame.box(0,0) * frame.box(1,1) * frame.box(2,2);
    sim.run(frame, 2000);
    const double V1 = frame.box(0,0) * frame.box(1,1) * frame.box(2,2);

    // No NaN explosion.
    for (int i = 0; i < n; ++i) {
        ASSERT_FALSE(std::isnan(frame.atoms.x[i])) << "atom " << i << " NaN";
    }
    // System started COMPRESSED; NPT must EXPAND to relieve pressure.
    EXPECT_GT(V1, V0) << "NPT did not grow over-compressed box: V0=" << V0
                      << " V1=" << V1;
}

// ── HaloExchangeForward::compute_bonded_extension ────────────────────
//
// Pure CPU unit tests for the topology-scan helper that tells
// HaloExchangeForward how much extra ghost skin to carve beyond the pair
// cutoff so every atom referenced by a locally-owned bonded tuple can be
// resolved on the owner rank.
//
// @see .claude/specs/parallel-md-engine.md § A.5 "Bonded interactions
// across rank boundaries": Ghost skin =
//     cutoff + δ_skin + bonded_extension
// with bonded_extension = max(bond_len, angle_span, dihedral_span).

TEST_F(CpuMdTest, HaloForwardBondedExtensionEmptyTopology) {
    // No bonds / angles / dihedrals → no extra skin required.
    cpu::Frame frame(10);
    const double ext =
        cpu::md::mpi::HaloExchangeForward::compute_bonded_extension(frame);
    EXPECT_DOUBLE_EQ(ext, 0.0);
}

TEST_F(CpuMdTest, HaloForwardBondedExtensionBondLength) {
    // Diatomic: single bond of length 1.5 Å. Extension = bond length.
    cpu::Frame frame(2);
    frame.atoms.x[0] = 0.0; frame.atoms.y[0] = 0.0; frame.atoms.z[0] = 0.0;
    frame.atoms.x[1] = 1.5; frame.atoms.y[1] = 0.0; frame.atoms.z[1] = 0.0;
    frame.bonds.atomi = xt::zeros<int>({1});
    frame.bonds.atomj = xt::zeros<int>({1});
    frame.bonds.type  = xt::zeros<int>({1});
    frame.bonds.atomi[0] = 0;
    frame.bonds.atomj[0] = 1;
    const double ext =
        cpu::md::mpi::HaloExchangeForward::compute_bonded_extension(frame);
    EXPECT_NEAR(ext, 1.5, 1e-12);
}

TEST_F(CpuMdTest, HaloForwardBondedExtensionAngleSpan) {
    // 3 atoms in a line, spaced 1.0 apart → angle-tuple span = 2.0.
    cpu::Frame frame(3);
    for (int i = 0; i < 3; ++i) {
        frame.atoms.x[i] = double(i);
        frame.atoms.y[i] = 0.0;
        frame.atoms.z[i] = 0.0;
    }
    frame.angles.atomi = xt::zeros<int>({1});
    frame.angles.atomj = xt::zeros<int>({1});
    frame.angles.atomk = xt::zeros<int>({1});
    frame.angles.type  = xt::zeros<int>({1});
    frame.angles.atomi[0] = 0;
    frame.angles.atomj[0] = 1;
    frame.angles.atomk[0] = 2;
    const double ext =
        cpu::md::mpi::HaloExchangeForward::compute_bonded_extension(frame);
    EXPECT_NEAR(ext, 2.0, 1e-12);
}

TEST_F(CpuMdTest, HaloForwardBondedExtensionDihedralSpan) {
    // 4 atoms in a trans-zigzag, bond = 1.529, angle = 112.7°.
    // Dihedral 1-to-4 straight-line span ≈ bond * sqrt(9 cos² + sin²) ≈ 3.9.
    cpu::Frame frame(4);
    const double bond = 1.529;
    const double half_apex = (kPI - 112.7 * kPI / 180.0) / 2.0;
    const double dx = bond * std::cos(half_apex);
    const double dy = bond * std::sin(half_apex);
    for (int i = 0; i < 4; ++i) {
        frame.atoms.x[i] = i * dx;
        frame.atoms.y[i] = (i % 2) ? dy : 0.0;
        frame.atoms.z[i] = 0.0;
    }
    frame.dihedrals.atomi = xt::zeros<int>({1});
    frame.dihedrals.atomj = xt::zeros<int>({1});
    frame.dihedrals.atomk = xt::zeros<int>({1});
    frame.dihedrals.atoml = xt::zeros<int>({1});
    frame.dihedrals.type  = xt::zeros<int>({1});
    frame.dihedrals.atomi[0] = 0;
    frame.dihedrals.atomj[0] = 1;
    frame.dihedrals.atomk[0] = 2;
    frame.dihedrals.atoml[0] = 3;
    const double ext =
        cpu::md::mpi::HaloExchangeForward::compute_bonded_extension(frame);
    // Expect the i-l straight-line span, roughly 3.9 Å for these params.
    EXPECT_GT(ext, 3.5) << "dihedral span too small: " << ext;
    EXPECT_LT(ext, 4.2) << "dihedral span too large: " << ext;
    // The widest pair in a tuple is always i-l, not i-j / i-k / j-k / j-l / k-l.
    const double expected = bond * std::sqrt(9.0 * std::cos(half_apex) *
                                              std::cos(half_apex)
                                              + std::sin(half_apex) *
                                                std::sin(half_apex));
    EXPECT_NEAR(ext, expected, 1e-10);
}

TEST_F(CpuMdTest, HaloForwardBondedExtensionMaxAcrossTupleKinds) {
    // Mixed topology: bond 1.0 + angle 3.0 + dihedral 2.0 → returns 3.0.
    cpu::Frame frame(4);
    frame.atoms.x = xt::xtensor<double, 1>({0.0, 1.0, 1.5, 3.5});
    frame.atoms.y = xt::zeros<double>({4});
    frame.atoms.z = xt::zeros<double>({4});

    // One bond (0-1, length 1.0)
    frame.bonds.atomi = xt::zeros<int>({1});
    frame.bonds.atomj = xt::zeros<int>({1});
    frame.bonds.type  = xt::zeros<int>({1});
    frame.bonds.atomi[0] = 0;
    frame.bonds.atomj[0] = 1;

    // One angle (0-1-3, i-to-k straight-line = 3.5; but we want angle max = 3.0 case)
    // Relax: atoms at 0, 1, 1.5 → angle 0-1-2 i-k = 1.5; use atoms 0-1-3 span 3.5
    // Actually we want the angle-max > dihedral-max to prove max is taken.
    // Redo atoms: 0, 1, 4, 2.5 → bond 0-1 = 1.0, angle 0-1-2 span 4.0, dihedral 0-1-2-3 span 2.5.
    frame.atoms.x = xt::xtensor<double, 1>({0.0, 1.0, 4.0, 2.5});
    // bond 0-1 length = 1.0
    // angle: atoms 0-1-2, i-k span = |x[0] - x[2]| = 4.0
    // dihedral: atoms 0-1-2-3, i-l span = |x[0] - x[3]| = 2.5
    frame.angles.atomi = xt::zeros<int>({1});
    frame.angles.atomj = xt::zeros<int>({1});
    frame.angles.atomk = xt::zeros<int>({1});
    frame.angles.type  = xt::zeros<int>({1});
    frame.angles.atomi[0] = 0;
    frame.angles.atomj[0] = 1;
    frame.angles.atomk[0] = 2;

    frame.dihedrals.atomi = xt::zeros<int>({1});
    frame.dihedrals.atomj = xt::zeros<int>({1});
    frame.dihedrals.atomk = xt::zeros<int>({1});
    frame.dihedrals.atoml = xt::zeros<int>({1});
    frame.dihedrals.type  = xt::zeros<int>({1});
    frame.dihedrals.atomi[0] = 0;
    frame.dihedrals.atomj[0] = 1;
    frame.dihedrals.atomk[0] = 2;
    frame.dihedrals.atoml[0] = 3;

    const double ext =
        cpu::md::mpi::HaloExchangeForward::compute_bonded_extension(frame);
    EXPECT_NEAR(ext, 4.0, 1e-12)
        << "expected max(bond=1, angle=4, dihedral=2.5) = 4, got " << ext;
}

TEST_F(CpuMdTest, HaloForwardBondedExtensionPeriodicBondWrap) {
    // Two atoms bonded across the periodic boundary: x=0.5 and x=L-0.5
    // look 79 Å apart without PBC, but the real bond is 1.0 Å. Without the
    // minimum-image fix in pair_dist_sq the bonded extension picks up the
    // unwrapped span and inflates the halo cutoff until every owned atom
    // ships as a ghost on every face (n_ghost explodes past capacity).
    // Regression guard for the debug_amorphous 2/4-rank crash.
    cpu::Frame frame(2);
    frame.box = Mat3::cubic(80.0);
    frame.atoms.x[0] = 0.5;  frame.atoms.y[0] = 0.0; frame.atoms.z[0] = 0.0;
    frame.atoms.x[1] = 79.5; frame.atoms.y[1] = 0.0; frame.atoms.z[1] = 0.0;
    frame.bonds.atomi = xt::zeros<int>({1});
    frame.bonds.atomj = xt::zeros<int>({1});
    frame.bonds.type  = xt::zeros<int>({1});
    frame.bonds.atomi[0] = 0;
    frame.bonds.atomj[0] = 1;
    const double ext =
        cpu::md::mpi::HaloExchangeForward::compute_bonded_extension(frame);
    EXPECT_NEAR(ext, 1.0, 1e-12)
        << "periodic bond picked up unwrapped span " << ext << " (box=80)";
}
