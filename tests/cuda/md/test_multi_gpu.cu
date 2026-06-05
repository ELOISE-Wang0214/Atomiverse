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
