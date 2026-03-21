#include <gtest/gtest.h>
#include <cuda_runtime.h>

#include "atomiverse/cuda/backend.cuh"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/facade/stage.h"

using namespace ATV;

TEST(Frame, HostCreate) {
    cpu::Frame frame(100);
    EXPECT_EQ(frame.n_atoms(), 100);
    EXPECT_EQ(frame.x.size(), 100u);
    EXPECT_EQ(frame.mass.size(), 100u);
}

TEST(Frame, Mat3Box) {
    cpu::Frame frame(2);
    frame.box = Mat3::cubic(10.0);
    EXPECT_DOUBLE_EQ(frame.box(0, 0), 10.0);
    EXPECT_DOUBLE_EQ(frame.box(1, 1), 10.0);
    EXPECT_DOUBLE_EQ(frame.box(2, 2), 10.0);
    EXPECT_DOUBLE_EQ(frame.box(0, 1), 0.0);
}

TEST(Frame, ToDeviceAndBack) {
    cpu::Frame frame(4);
    frame.x = xt::xtensor<Real, 1>({1.0, 2.0, 3.0, 4.0});
    frame.mass = xt::xtensor<Real, 1>({1.0, 1.0, 1.0, 1.0});
    frame.box = Mat3::cubic(10.0);

    auto d_frame = cpu::to_cuda(frame);
    EXPECT_NE(d_frame.d_n_atoms, nullptr);

    frame.x.fill(Real(0));
    cuda::to_cpu(d_frame, frame);
    EXPECT_DOUBLE_EQ(frame.x[0], 1.0);
    EXPECT_DOUBLE_EQ(frame.x[3], 4.0);

    cuda::destroy(d_frame);
}

TEST(Frame, TopologyRoundTrip) {
    cpu::Frame frame(4);
    frame.x = xt::xtensor<Real, 1>({0.0, 1.53, 3.06, 4.59});
    frame.y = xt::xtensor<Real, 1>({0.0, 0.0, 0.0, 0.0});
    frame.z = xt::xtensor<Real, 1>({0.0, 0.0, 0.0, 0.0});
    frame.mass = xt::xtensor<Real, 1>({14.0, 14.0, 14.0, 14.0});
    frame.box = Mat3::cubic(20.0);

    // 3 bonds: 0-1, 1-2, 2-3
    frame.bond_i = xt::xtensor<Int, 1>({0, 1, 2});
    frame.bond_j = xt::xtensor<Int, 1>({1, 2, 3});
    frame.bond_type = xt::xtensor<Int, 1>({0, 0, 0});

    // 2 angles: 0-1-2, 1-2-3
    frame.angle_i = xt::xtensor<Int, 1>({0, 1});
    frame.angle_j = xt::xtensor<Int, 1>({1, 2});
    frame.angle_k = xt::xtensor<Int, 1>({2, 3});
    frame.angle_type = xt::xtensor<Int, 1>({0, 0});

    // 1 dihedral: 0-1-2-3
    frame.dihedral_i = xt::xtensor<Int, 1>({0});
    frame.dihedral_j = xt::xtensor<Int, 1>({1});
    frame.dihedral_k = xt::xtensor<Int, 1>({2});
    frame.dihedral_l = xt::xtensor<Int, 1>({3});
    frame.dihedral_type = xt::xtensor<Int, 1>({0});

    auto d_frame = cpu::to_cuda(frame);

    // Verify device pointers allocated
    EXPECT_NE(d_frame.bonds.i, nullptr);
    EXPECT_NE(d_frame.bonds.n_bonds, nullptr);
    EXPECT_NE(d_frame.angles.i, nullptr);
    EXPECT_NE(d_frame.angles.n_angles, nullptr);
    EXPECT_NE(d_frame.dihedrals.i, nullptr);
    EXPECT_NE(d_frame.dihedrals.n_dihedrals, nullptr);
    EXPECT_NE(d_frame.stress, nullptr);

    // Verify counts via D2H
    Int nb = 0, na = 0, nd = 0;
    cudaMemcpy(&nb, d_frame.bonds.n_bonds, sizeof(Int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&na, d_frame.angles.n_angles, sizeof(Int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&nd, d_frame.dihedrals.n_dihedrals, sizeof(Int), cudaMemcpyDeviceToHost);
    EXPECT_EQ(nb, 3);
    EXPECT_EQ(na, 2);
    EXPECT_EQ(nd, 1);

    // Verify bond indices round-trip
    std::vector<Int> host_bond_i(3), host_bond_j(3);
    cudaMemcpy(host_bond_i.data(), d_frame.bonds.i, 3 * sizeof(Int), cudaMemcpyDeviceToHost);
    cudaMemcpy(host_bond_j.data(), d_frame.bonds.j, 3 * sizeof(Int), cudaMemcpyDeviceToHost);
    EXPECT_EQ(host_bond_i[0], 0); EXPECT_EQ(host_bond_j[0], 1);
    EXPECT_EQ(host_bond_i[1], 1); EXPECT_EQ(host_bond_j[1], 2);
    EXPECT_EQ(host_bond_i[2], 2); EXPECT_EQ(host_bond_j[2], 3);

    // Verify host-cached metadata in cuda::Frame
    d_frame.n_atoms_ = 4;
    d_frame.n_bonds_ = 3;
    d_frame.n_angles_ = 2;
    d_frame.n_dihedrals_ = 1;
    EXPECT_EQ(d_frame.n_atoms(), 4);
    EXPECT_EQ(d_frame.n_bonds(), 3);
    EXPECT_EQ(d_frame.n_angles(), 2);
    EXPECT_EQ(d_frame.n_dihedrals(), 1);

    cuda::destroy(d_frame);
}

TEST(Frame, EmptyTopologyRoundTrip) {
    // cpu::Frame with no topology should still work
    cpu::Frame frame(4);
    frame.x = xt::xtensor<Real, 1>({1.0, 2.0, 3.0, 4.0});
    frame.mass = xt::xtensor<Real, 1>({1.0, 1.0, 1.0, 1.0});
    frame.box = Mat3::cubic(10.0);

    auto d_frame = cpu::to_cuda(frame);
    // Topology pointers should be nullptr
    EXPECT_EQ(d_frame.bonds.i, nullptr);
    EXPECT_EQ(d_frame.angles.i, nullptr);
    EXPECT_EQ(d_frame.dihedrals.i, nullptr);

    // Round-trip positions still work
    frame.x.fill(Real(0));
    cuda::to_cpu(d_frame, frame);
    EXPECT_DOUBLE_EQ(frame.x[0], 1.0);
    EXPECT_DOUBLE_EQ(frame.x[3], 4.0);

    cuda::destroy(d_frame);
}

TEST(Frame, ZeroForces) {
    cpu::Frame frame(32);
    frame.mass.fill(Real(1.0));
    frame.box = Mat3::cubic(10.0);

    // Use Driver-level lifecycle
    Backend backend;
    backend.init();

    backend.materialize(frame);
    auto& d_frame = backend.frame();

    // Set forces to 1 via host copy
    std::vector<Real> ones(32, Real(1));
    cudaMemcpy(d_frame.fx, ones.data(), 32 * sizeof(Real), cudaMemcpyHostToDevice);

    // Build + launch ZeroForces via Plan
    md::ZeroForces zf;
    PlanBuilder builder;
    builder.add({0, 0, &zf});
    auto plan = builder.build(static_cast<int>(md::Phase::COUNT));
    backend.load(plan);

    BuildContext bctx{backend.context(), d_frame};
    zf.build(bctx);

    backend.step(0);
    backend.sync();

    DestroyContext dctx{backend.context(), d_frame};
    zf.destroy(dctx);

    cuda::to_cpu(d_frame, frame);
    for (int i = 0; i < 32; ++i)
        EXPECT_DOUBLE_EQ(frame.fx[i], 0.0);

    backend.shutdown();
}
