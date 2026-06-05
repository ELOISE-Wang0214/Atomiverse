/// @file tests/cuda/md/test_multi_gpu.cu
/// @brief cuda::md::MultiGpu — partition a global frame across GPUs + gather
///        back; atom conservation and position round-trip.

#include <gtest/gtest.h>
#include <cmath>

#include "atomiverse/core/mat3.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/md/multi_gpu.cuh"

using namespace ATV;

namespace {
// 4x4x4 grid of atoms in a cubic box, deterministic tags 0..63.
cpu::Frame make_grid(double L) {
    const int s = 4, n = s * s * s;
    cpu::Frame f(n);
    f.box = Mat3::cubic(L);
    double sp = L / s;
    int idx = 0;
    for (int ix = 0; ix < s; ++ix)
        for (int iy = 0; iy < s; ++iy)
            for (int iz = 0; iz < s; ++iz) {
                f.atoms.x[idx] = (ix + 0.5) * sp;
                f.atoms.y[idx] = (iy + 0.5) * sp;
                f.atoms.z[idx] = (iz + 0.5) * sp;
                f.atoms.vx[idx] = 0.01 * idx;
                f.atoms.mass[idx] = 1.0;
                f.atoms.type[idx] = 0;
                f.atoms.tag[idx] = idx;
                ++idx;
            }
    return f;
}
} // namespace

TEST(MultiGpu, PartitionConservesAtoms) {
    auto g = make_grid(12.0);
    cuda::md::MultiGpu mg;
    mg.partition(g);
    int total = 0;
    for (int s = 0; s < mg.count(); ++s)
        total += mg.n_local(s);
    EXPECT_EQ(total, g.atoms.size()); // every atom owned by exactly one GPU
}

TEST(MultiGpu, PartitionGatherRoundTrip) {
    auto g = make_grid(12.0);
    cuda::md::MultiGpu mg;
    mg.partition(g);

    cpu::Frame out(g.atoms.size());
    for (int i = 0; i < g.atoms.size(); ++i)
        out.atoms.tag[i] = i;
    mg.gather(out);

    // cuda::Frame stores per-atom data as real_t (float), so the
    // partition->materialize->gather round trip is exact only to float
    // precision (e.g. vx = 0.01*idx is not representable in float).
    for (int i = 0; i < g.atoms.size(); ++i) {
        EXPECT_NEAR(out.atoms.x[i], g.atoms.x[i], 1e-5) << "atom " << i;
        EXPECT_NEAR(out.atoms.y[i], g.atoms.y[i], 1e-5);
        EXPECT_NEAR(out.atoms.z[i], g.atoms.z[i], 1e-5);
        EXPECT_NEAR(out.atoms.vx[i], g.atoms.vx[i], 1e-5);
    }
}

TEST(MultiGpu, HaloExchangeXSplit) {
    // 2-GPU x-split (dims {2,1,1}); a single forward halo pass is complete.
    cuda::md::MultiGpu mg(2);
    if (mg.count() < 2)
        GTEST_SKIP() << "needs >= 2 GPUs";

    auto g = make_grid(12.0); // grid coords 1.5/4.5/7.5/10.5
    mg.partition(g);          // creates the DomainDecomposition (dd_)
    ASSERT_EQ(mg.domain().dims(),
              (std::array<int, 3>{2, 1, 1})); // boundary at x = 6
    mg.exchange_halo(3.0);                    // cutoff 3.0

    // Each slab owns 32 atoms; boundary atoms within 3.0 of x=0/6/12 are sent
    // both ways, so each side receives the other's near-boundary atoms.
    EXPECT_GT(mg.n_ghost(0), 0);
    EXPECT_EQ(mg.n_ghost(0), mg.n_ghost(1)); // symmetric exchange

    // Ghosts on slab 0 must be slab-1 atoms (x in [6,12)); on slab 1, x in [0,6).
    for (int slot = 0; slot < 2; ++slot) {
        const auto& f = mg.frame(slot);
        int ng = mg.n_ghost(slot);
        std::vector<float> gx(ng);
        cudaSetDevice(mg.devices().device_id(slot));
        cudaMemcpy(gx.data(), f.atoms.x + mg.n_local(slot), ng * sizeof(float),
                   cudaMemcpyDeviceToHost);
        for (int k = 0; k < ng; ++k) {
            if (slot == 0)
                EXPECT_GE(gx[k], 6.0f) << "slab0 ghost should be slab1 atom";
            else
                EXPECT_LT(gx[k], 6.0f) << "slab1 ghost should be slab0 atom";
        }
    }
}

namespace {
// Shift the x-coordinate of every owned atom by `delta` (test helper for
// forward_comm: emulates one integration step moving the owners).
__global__ void shift_owned_x(float* x, int n_local, float delta) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n_local)
        x[i] += delta;
}
} // namespace

TEST(MultiGpu, ForwardCommRefreshesGhostPositions) {
    // build_halo establishes the ghost SET + a gather plan; forward_comm then
    // refreshes only ghost POSITIONS from the current owners. Moving every
    // owned atom by a fixed delta in x must move every ghost by the same delta
    // (a ghost is a copy of some owner), with y/z and the ghost count unchanged.
    cuda::md::MultiGpu mg(2);
    if (mg.count() < 2)
        GTEST_SKIP() << "needs >= 2 GPUs";
    auto g = make_grid(12.0);
    mg.partition(g);
    mg.build_halo(3.0);

    auto read_ghosts = [&](int slot, std::vector<float>& gx,
                           std::vector<float>& gy, std::vector<float>& gz) {
        int ng = mg.n_ghost(slot);
        gx.resize(ng);
        gy.resize(ng);
        gz.resize(ng);
        const auto& f = mg.frame(slot);
        int off = mg.n_local(slot);
        cudaSetDevice(mg.devices().device_id(slot));
        cudaMemcpy(gx.data(), f.atoms.x + off, ng * sizeof(float),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(gy.data(), f.atoms.y + off, ng * sizeof(float),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(gz.data(), f.atoms.z + off, ng * sizeof(float),
                   cudaMemcpyDeviceToHost);
    };

    std::vector<std::vector<float>> bx(2), by(2), bz(2);
    for (int s = 0; s < 2; ++s)
        read_ghosts(s, bx[s], by[s], bz[s]);

    const float delta = 0.123f;
    for (int s = 0; s < 2; ++s) {
        int nl = mg.n_local(s);
        cudaSetDevice(mg.devices().device_id(s));
        int block = 128, grid = (nl + block - 1) / block;
        shift_owned_x<<<grid, block>>>(mg.frame(s).atoms.x, nl, delta);
    }
    mg.forward_comm();

    for (int s = 0; s < 2; ++s) {
        std::vector<float> ax, ay, az;
        read_ghosts(s, ax, ay, az);
        ASSERT_EQ(ax.size(), bx[s].size()) << "ghost count changed on slot " << s;
        for (size_t k = 0; k < ax.size(); ++k) {
            EXPECT_NEAR(ax[k], bx[s][k] + delta, 1e-5) << "slot " << s << " k " << k;
            EXPECT_NEAR(ay[k], by[s][k], 1e-5);
            EXPECT_NEAR(az[k], bz[s][k], 1e-5);
        }
    }
}

TEST(MultiGpu, LjForceParityVsSingleGpu) {
    // The crux of multi-GPU correctness: an owned atom must feel the SAME LJ
    // force whether the system runs on one GPU or is split across GPUs with a
    // halo. Same brute-force kernel both ways; only the atom set in each owned
    // atom's loop differs (single: all; multi: owned+ghost). With halo cutoff
    // == LJ cutoff, the ghost set captures every cross-boundary neighbor, so
    // the per-atom force sums are identical (up to float summation order).
    auto g = make_grid(12.0); // spacing 3.0
    const int n = g.atoms.size();
    // Break the perfect-lattice symmetry (else neighbor forces cancel to ~0
    // and the parity check is vacuous). Deterministic per-atom jitter.
    for (int i = 0; i < n; ++i) {
        g.atoms.x[i] += 0.4 * std::sin(0.7 * i + 0.1);
        g.atoms.y[i] += 0.4 * std::sin(1.3 * i + 0.2);
        g.atoms.z[i] += 0.4 * std::sin(1.9 * i + 0.3);
    }
    const double cut = 3.5, eps = 1.0, sig = 2.5;

    cpu::Frame ref(n);
    for (int i = 0; i < n; ++i)
        ref.atoms.tag[i] = i;
    {
        cuda::md::MultiGpu mg1(1); // whole system on one GPU
        mg1.partition(g);
        mg1.exchange_halo(cut); // no ghosts (self-neighbors skipped)
        mg1.compute_lj_bruteforce(cut, eps, sig);
        mg1.gather_forces(ref);
    }

    cuda::md::MultiGpu mg(2);
    if (mg.count() < 2)
        GTEST_SKIP() << "needs >= 2 GPUs";
    mg.partition(g);
    mg.exchange_halo(cut);
    mg.compute_lj_bruteforce(cut, eps, sig);
    cpu::Frame mf(n);
    for (int i = 0; i < n; ++i)
        mf.atoms.tag[i] = i;
    mg.gather_forces(mf);

    double maxabs = 0;
    for (int i = 0; i < n; ++i)
        maxabs = std::max(maxabs, std::abs(ref.atoms.fx[i]));
    EXPECT_GT(maxabs, 1e-3) << "reference forces should be nonzero";
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(mf.atoms.fx[i], ref.atoms.fx[i], 1e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.fy[i], ref.atoms.fy[i], 1e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.fz[i], ref.atoms.fz[i], 1e-3) << "atom " << i;
    }
}

namespace {
// s x s x s grid in a cubic box of edge L, deterministic tags.
cpu::Frame make_grid_n(int s, double L) {
    const int n = s * s * s;
    cpu::Frame f(n);
    f.box = Mat3::cubic(L);
    double sp = L / s;
    int idx = 0;
    for (int ix = 0; ix < s; ++ix)
        for (int iy = 0; iy < s; ++iy)
            for (int iz = 0; iz < s; ++iz) {
                f.atoms.x[idx] = (ix + 0.5) * sp;
                f.atoms.y[idx] = (iy + 0.5) * sp;
                f.atoms.z[idx] = (iz + 0.5) * sp;
                f.atoms.mass[idx] = 1.0;
                f.atoms.type[idx] = 0;
                f.atoms.tag[idx] = idx;
                ++idx;
            }
    return f;
}
} // namespace

TEST(MultiGpu, LjForceParityDims2Overlap) {
    // Regression for the dims==2 double-send bug: along an axis split into
    // exactly 2 subdomains, both faces point to the SAME neighbor (periodic
    // self-pair). An atom in the overlap band (within cutoff of BOTH shared
    // boundaries) used to be sent TWICE, and the min-image force kernel then
    // counted it twice -> doubled force. Here spacing 2.0 with cutoff 4.5 puts
    // many atoms in the overlap band AND makes them genuine cross-boundary
    // neighbors, so the bug would break parity. The dedup (mode=2 single send)
    // restores it.
    auto g = make_grid_n(6, 16.0); // box/2 = 8; slab width 8, overlap band [3.5,4.5]
    const int n = g.atoms.size();  // spacing 2.667 > sigma -> gentle forces
    for (int i = 0; i < n; ++i) {
        g.atoms.x[i] += 0.3 * std::sin(0.7 * i + 0.1);
        g.atoms.y[i] += 0.3 * std::sin(1.3 * i + 0.2);
        g.atoms.z[i] += 0.3 * std::sin(1.9 * i + 0.3);
    }
    const double cut = 4.5, eps = 1.0, sig = 2.5;

    cpu::Frame ref(n);
    for (int i = 0; i < n; ++i)
        ref.atoms.tag[i] = i;
    {
        cuda::md::MultiGpu mg1(1);
        mg1.partition(g, 4.0); // generous ghost headroom for the dense band
        mg1.build_halo(cut);
        mg1.compute_lj_bruteforce(cut, eps, sig);
        mg1.gather_forces(ref);
    }

    cuda::md::MultiGpu mg(2);
    if (mg.count() < 2)
        GTEST_SKIP() << "needs >= 2 GPUs";
    mg.partition(g, 4.0);
    ASSERT_EQ(mg.domain().dims(), (std::array<int, 3>{2, 1, 1}));
    mg.build_halo(cut);
    cpu::Frame mf(n);
    for (int i = 0; i < n; ++i)
        mf.atoms.tag[i] = i;
    mg.compute_lj_bruteforce(cut, eps, sig);
    mg.gather_forces(mf);

    double maxabs = 0;
    for (int i = 0; i < n; ++i)
        maxabs = std::max(maxabs, std::abs(ref.atoms.fx[i]));
    EXPECT_GT(maxabs, 1e-3) << "reference forces should be nonzero";
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(mf.atoms.fx[i], ref.atoms.fx[i], 2e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.fy[i], ref.atoms.fy[i], 2e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.fz[i], ref.atoms.fz[i], 2e-3) << "atom " << i;
    }
}

TEST(MultiGpu, LjForceParityVsSingleGpu4Way) {
    // The 4-GPU decomposition is 2x2x1, so atoms near the (x=6, y=6) edge have
    // neighbors owned by the DIAGONAL subdomain. A single forward halo pass
    // would never deliver those corner atoms; only the 3-axis cascade (x-pass
    // then y-pass forwarding x-ghosts) does. This test fails without the
    // cascade and passes with it — owned-atom LJ forces must match single-GPU.
    cuda::md::MultiGpu mg(4);
    if (mg.count() < 4)
        GTEST_SKIP() << "needs 4 GPUs";
    ASSERT_EQ(mg.devices().count(), 4);

    auto g = make_grid(12.0);
    const int n = g.atoms.size();
    for (int i = 0; i < n; ++i) {
        g.atoms.x[i] += 0.4 * std::sin(0.7 * i + 0.1);
        g.atoms.y[i] += 0.4 * std::sin(1.3 * i + 0.2);
        g.atoms.z[i] += 0.4 * std::sin(1.9 * i + 0.3);
    }
    const double cut = 3.5, eps = 1.0, sig = 2.5;

    cpu::Frame ref(n);
    for (int i = 0; i < n; ++i)
        ref.atoms.tag[i] = i;
    {
        cuda::md::MultiGpu mg1(1);
        mg1.partition(g);
        mg1.build_halo(cut);
        mg1.compute_lj_bruteforce(cut, eps, sig);
        mg1.gather_forces(ref);
    }

    mg.partition(g); // creates the DomainDecomposition (dd_)

    // Confirm the decomposition genuinely needs the cascade (>1 axis split).
    const auto dims = mg.domain().dims();
    int split_axes = (dims[0] > 1) + (dims[1] > 1) + (dims[2] > 1);
    ASSERT_GE(split_axes, 2) << "4-GPU layout should split >=2 axes";

    mg.build_halo(cut);
    mg.compute_lj_bruteforce(cut, eps, sig);
    cpu::Frame mf(n);
    for (int i = 0; i < n; ++i)
        mf.atoms.tag[i] = i;
    mg.gather_forces(mf);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(mf.atoms.fx[i], ref.atoms.fx[i], 1e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.fy[i], ref.atoms.fy[i], 1e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.fz[i], ref.atoms.fz[i], 1e-3) << "atom " << i;
    }
}

TEST(MultiGpu, NveTrajectoryParity) {
    // Run a short NVE trajectory single-GPU vs split-across-GPUs and require the
    // owned positions to agree — validating the full step loop (halo + force +
    // integrate) end to end. Small dt / few steps so atoms stay in their slab
    // (no migration) and float chaos stays bounded.
    auto g = make_grid(12.0);
    const int n = g.atoms.size();
    for (int i = 0; i < n; ++i) {
        g.atoms.x[i] += 0.4 * std::sin(0.7 * i + 0.1);
        g.atoms.y[i] += 0.4 * std::sin(1.3 * i + 0.2);
        g.atoms.z[i] += 0.4 * std::sin(1.9 * i + 0.3);
        g.atoms.vx[i] = 0.02 * std::sin(0.5 * i);
        g.atoms.vy[i] = 0.02 * std::sin(0.9 * i);
        g.atoms.vz[i] = 0.02 * std::sin(1.1 * i);
        g.atoms.mass[i] = 1.0;
    }
    const double cut = 3.5, eps = 1.0, sig = 2.5, dt = 0.002;
    const int K = 10;

    cpu::Frame ref(n);
    for (int i = 0; i < n; ++i)
        ref.atoms.tag[i] = i;
    {
        cuda::md::MultiGpu mg1(1);
        mg1.partition(g);
        mg1.run_nve(K, dt, cut, eps, sig);
        mg1.gather(ref);
    }

    cuda::md::MultiGpu mg(2);
    if (mg.count() < 2)
        GTEST_SKIP() << "needs >= 2 GPUs";
    mg.partition(g);
    mg.run_nve(K, dt, cut, eps, sig);
    cpu::Frame mf(n);
    for (int i = 0; i < n; ++i)
        mf.atoms.tag[i] = i;
    mg.gather(mf);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(mf.atoms.x[i], ref.atoms.x[i], 5e-3) << "atom " << i;
        EXPECT_NEAR(mf.atoms.y[i], ref.atoms.y[i], 5e-3);
        EXPECT_NEAR(mf.atoms.z[i], ref.atoms.z[i], 5e-3);
    }
}

TEST(MultiGpu, EachAtomOnOwnerSlab) {
    auto g = make_grid(12.0);
    cuda::md::MultiGpu mg;
    mg.partition(g);
    if (mg.count() < 2)
        GTEST_SKIP() << "single GPU: trivially one slab";
    // Owner count per slab should match DomainDecomposition::owner_of.
    const auto& dd = mg.domain();
    std::vector<int> expect(mg.count(), 0);
    for (int i = 0; i < g.atoms.size(); ++i)
        expect[dd.owner_of(g.atoms.x[i], g.atoms.y[i], g.atoms.z[i])]++;
    for (int s = 0; s < mg.count(); ++s)
        EXPECT_EQ(mg.n_local(s), expect[s]) << "slab " << s;
}
