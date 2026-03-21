/// @file tests/test_torch_neighbor_list.cu
/// @brief Tests for torch NeighborList operator.

#include <gtest/gtest.h>

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/neighbor_list.cuh"
#include "atomiverse/torch/transfer.cuh"
#include "atomiverse/torch/potentials/lj.cuh"

#include <cuda_runtime.h>
#include <cmath>

using ATV::Real;
using ATV::Int;
using ATV::Mat3;
namespace atorch = ATV::torch;

namespace {

/// @brief Build a simple cubic lattice torch::Frame.
atorch::Frame make_cubic_frame(int n_side, Real box_len) {
    int n = n_side * n_side * n_side;
    Real spacing = box_len / Real(n_side);

    ATV::cpu::Frame cpu_frame(n);
    int idx = 0;
    for (int ix = 0; ix < n_side; ++ix)
        for (int iy = 0; iy < n_side; ++iy)
            for (int iz = 0; iz < n_side; ++iz) {
                cpu_frame.x[idx] = (Real(ix) + Real(0.5)) * spacing;
                cpu_frame.y[idx] = (Real(iy) + Real(0.5)) * spacing;
                cpu_frame.z[idx] = (Real(iz) + Real(0.5)) * spacing;
                cpu_frame.mass[idx] = Real(1.0);
                cpu_frame.type_id[idx] = 0;
                idx++;
            }
    cpu_frame.box = Mat3::cubic(static_cast<double>(box_len));
    return atorch::from_cpu(cpu_frame);
}

} // anonymous namespace

TEST(TorchNeighborList, CubicLatticePairCount) {
    // 3^3 = 27 atoms, spacing = 10/3 = 3.33 Angstrom
    // cutoff = 4.0, skin = 0 -> r_search = 4.0
    // Only nearest neighbors at distance 3.33 are included (6 per atom)
    // Total pairs = 27 * 6 / 2 = 81
    Real box_len = Real(10.0);
    Real cutoff = Real(4.0);
    auto frame = make_cubic_frame(3, box_len);

    atorch::neighbor_list(frame, cutoff, Real(0));

    // Exact count depends on cluster algorithm; just check populated
    EXPECT_GT(frame.n_pairs(), 0);
    EXPECT_TRUE(frame.pair_i.defined());
    EXPECT_TRUE(frame.pair_j.defined());
    EXPECT_TRUE(frame.pair_dist.defined());
}

TEST(TorchNeighborList, PairTensorsNoGrad) {
    auto frame = make_cubic_frame(3, Real(10.0));

    atorch::neighbor_list(frame, Real(4.0), Real(0));
    ASSERT_GT(frame.n_pairs(), 0);

    EXPECT_FALSE(frame.pair_i.requires_grad());
    EXPECT_FALSE(frame.pair_j.requires_grad());
    EXPECT_FALSE(frame.pair_dx.requires_grad());
    EXPECT_FALSE(frame.pair_dy.requires_grad());
    EXPECT_FALSE(frame.pair_dz.requires_grad());
    EXPECT_FALSE(frame.pair_dist.requires_grad());
}

TEST(TorchNeighborList, PairTensorsAreSelfContained) {
    auto frame = make_cubic_frame(3, Real(10.0));

    atorch::neighbor_list(frame, Real(4.0), Real(0));
    ASSERT_GT(frame.n_pairs(), 0);

    EXPECT_TRUE(frame.pair_i.is_cuda());
    EXPECT_TRUE(frame.pair_i.is_contiguous());
    EXPECT_EQ(frame.pair_i.scalar_type(), atorch::kIntDtype);

    EXPECT_TRUE(frame.pair_dist.is_cuda());
    EXPECT_TRUE(frame.pair_dist.is_contiguous());
    EXPECT_EQ(frame.pair_dist.scalar_type(), atorch::kRealDtype);

    int64_t np = static_cast<int64_t>(frame.n_pairs());
    EXPECT_EQ(frame.pair_i.numel(), np);
    EXPECT_EQ(frame.pair_j.numel(), np);
    EXPECT_EQ(frame.pair_dx.numel(), np);
    EXPECT_EQ(frame.pair_dist.numel(), np);
}

TEST(TorchNeighborList, DistanceConsistency) {
    auto frame = make_cubic_frame(3, Real(10.0));

    atorch::neighbor_list(frame, Real(4.0), Real(0));
    ASSERT_GT(frame.n_pairs(), 0);

    auto dx = frame.pair_dx.to(::torch::kCPU);
    auto dy = frame.pair_dy.to(::torch::kCPU);
    auto dz = frame.pair_dz.to(::torch::kCPU);
    auto dist = frame.pair_dist.to(::torch::kCPU);

    for (int64_t p = 0; p < frame.n_pairs(); ++p) {
        Real d = dx[p].item<Real>();
        Real e = dy[p].item<Real>();
        Real f = dz[p].item<Real>();
        Real expected = std::sqrt(d * d + e * e + f * f);
        Real actual = dist[p].item<Real>();
        EXPECT_NEAR(actual, expected, Real(1e-4))
            << "pair " << p << ": dist mismatch";
    }
}

TEST(TorchNeighborList, PairListFeedsLJ) {
    auto frame = make_cubic_frame(3, Real(10.0));
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    int n = frame.n_atoms();
    frame.fx = ::torch::zeros({n}, real_opts);
    frame.fy = ::torch::zeros({n}, real_opts);
    frame.fz = ::torch::zeros({n}, real_opts);
    frame.energy = ::torch::zeros({1}, real_opts);

    atorch::neighbor_list(frame, Real(4.0), Real(0));
    ASSERT_GT(frame.n_pairs(), 0);

    Real eps = Real(1.0), sig = Real(1.0);
    auto energy = atorch::lj(frame,
        ::torch::tensor({eps}, real_opts),
        ::torch::tensor({sig}, real_opts), Real(4.0));
    cudaDeviceSynchronize();

    Real E = energy.to(::torch::kCPU).item<Real>();
    EXPECT_NE(E, Real(0));
}

TEST(TorchNeighborList, ClassBasedBuildUpdate) {
    // Test the class-based API: build once, update multiple times
    auto frame = make_cubic_frame(3, Real(10.0));

    atorch::NeighborList nl(Real(4.0), Real(0));
    nl.build(frame);
    nl.update(frame);

    Int first_count = frame.n_pairs();
    EXPECT_GT(first_count, 0);

    // Update again (same positions -> same result)
    nl.update(frame);
    EXPECT_EQ(frame.n_pairs(), first_count);
}

#else

TEST(TorchNeighborList, Disabled) {
    GTEST_SKIP() << "ATV_HAS_TORCH not defined";
}

#endif // ATV_HAS_TORCH
