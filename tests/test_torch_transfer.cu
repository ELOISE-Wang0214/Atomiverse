/// @file tests/test_torch_transfer.cu
/// @brief Tests for torch transfer functions.

#include <gtest/gtest.h>

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/transfer.cuh"
#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/torch/stream_guard.cuh"

#include <cuda_runtime.h>

using ATV::Real;
using ATV::Int;
using ATV::Mat3;
namespace atorch = ATV::torch;

class TorchTransferTest : public ::testing::Test {
protected:
    void SetUp() override {
        cpu_frame_ = std::make_unique<ATV::cpu::Frame>(4);
        cpu_frame_->x = xt::xtensor<Real, 1>({1.0, 2.0, 3.0, 4.0});
        cpu_frame_->y = xt::xtensor<Real, 1>({0.1, 0.2, 0.3, 0.4});
        cpu_frame_->z = xt::xtensor<Real, 1>({0.0, 0.0, 0.0, 0.0});
        cpu_frame_->mass = xt::xtensor<Real, 1>({12.0, 12.0, 12.0, 12.0});
        cpu_frame_->type_id = xt::xtensor<Int, 1>({0, 0, 1, 1});
        cpu_frame_->box = Mat3::cubic(10.0);
        cpu_frame_->energy = Real(1.5);
    }

    std::unique_ptr<ATV::cpu::Frame> cpu_frame_;
};

TEST_F(TorchTransferTest, FromCpuBasic) {
    auto tframe = atorch::from_cpu(*cpu_frame_);
    EXPECT_EQ(tframe.n_atoms(), 4);
    EXPECT_TRUE(tframe.x.is_cuda());
    EXPECT_EQ(tframe.x.scalar_type(), atorch::kRealDtype);
    EXPECT_EQ(tframe.x.numel(), 4);

    auto x_cpu = tframe.x.to(::torch::kCPU);
    EXPECT_NEAR(x_cpu[0].item<Real>(), Real(1.0), Real(1e-6));
    EXPECT_NEAR(x_cpu[3].item<Real>(), Real(4.0), Real(1e-6));
}

TEST_F(TorchTransferTest, ToCpuRoundTrip) {
    auto tframe = atorch::from_cpu(*cpu_frame_);

    ATV::cpu::Frame dst(4);
    atorch::to_cpu(tframe, dst);

    EXPECT_NEAR(dst.x[0], Real(1.0), Real(1e-6));
    EXPECT_NEAR(dst.x[3], Real(4.0), Real(1e-6));
    EXPECT_NEAR(dst.y[1], Real(0.2), Real(1e-6));
    EXPECT_NEAR(dst.mass[0], Real(12.0), Real(1e-6));
}

TEST_F(TorchTransferTest, ToCudaZeroCopy) {
    auto cuda_frame = ATV::cpu::to_cuda(*cpu_frame_);

    auto tframe = atorch::to_torch(cuda_frame);
    EXPECT_EQ(tframe.n_atoms(), 4);
    EXPECT_TRUE(tframe.x.is_cuda());

    // Zero-copy: tensor data_ptr should equal cuda::Frame pointer
    EXPECT_EQ(tframe.x.data_ptr<Real>(), cuda_frame.x);
    EXPECT_EQ(tframe.y.data_ptr<Real>(), cuda_frame.y);
    EXPECT_EQ(tframe.fx.data_ptr<Real>(), cuda_frame.fx);

    ATV::cuda::destroy(cuda_frame);
}

TEST_F(TorchTransferTest, ToCudaExtractPointers) {
    auto tframe = atorch::from_cpu(*cpu_frame_);

    ATV::cuda::Frame cframe;
    atorch::to_cuda(tframe, cframe);

    EXPECT_EQ(cframe.x, tframe.x.data_ptr<Real>());
    EXPECT_EQ(cframe.n_atoms(), 4);
}

TEST_F(TorchTransferTest, DestroyResetsFrame) {
    auto tframe = atorch::from_cpu(*cpu_frame_);
    EXPECT_TRUE(tframe.x.defined());

    atorch::destroy(tframe);
    EXPECT_FALSE(tframe.x.defined());
    EXPECT_EQ(tframe.n_atoms(), 0);
}

TEST_F(TorchTransferTest, StreamGuardRestores) {
    cudaStream_t s;
    cudaStreamCreate(&s);

    auto before = ::c10::cuda::getCurrentCUDAStream();
    {
        atorch::StreamGuard guard(static_cast<void*>(s));
        auto during = ::c10::cuda::getCurrentCUDAStream();
        EXPECT_EQ(during.stream(), s);
    }
    auto after = ::c10::cuda::getCurrentCUDAStream();
    EXPECT_EQ(after.stream(), before.stream());

    cudaStreamDestroy(s);
}

TEST_F(TorchTransferTest, TopologyRoundTrip) {
    // Add 2 bonds to the cpu::Frame
    cpu_frame_->bond_i = xt::xtensor<Int, 1>({0, 1});
    cpu_frame_->bond_j = xt::xtensor<Int, 1>({1, 2});
    cpu_frame_->bond_type = xt::xtensor<Int, 1>({0, 1});

    auto tframe = atorch::from_cpu(*cpu_frame_);
    EXPECT_EQ(tframe.n_bonds(), 2);

    // Verify bond_i
    EXPECT_TRUE(tframe.bond_i.defined());
    EXPECT_TRUE(tframe.bond_i.is_cuda());
    EXPECT_EQ(tframe.bond_i.numel(), 2);
    auto bi_cpu = tframe.bond_i.to(::torch::kCPU);
    EXPECT_EQ(bi_cpu[0].item<Int>(), Int(0));
    EXPECT_EQ(bi_cpu[1].item<Int>(), Int(1));

    // Verify bond_j
    auto bj_cpu = tframe.bond_j.to(::torch::kCPU);
    EXPECT_EQ(bj_cpu[0].item<Int>(), Int(1));
    EXPECT_EQ(bj_cpu[1].item<Int>(), Int(2));

    // Verify bond_type
    auto bt_cpu = tframe.bond_type.to(::torch::kCPU);
    EXPECT_EQ(bt_cpu[0].item<Int>(), Int(0));
    EXPECT_EQ(bt_cpu[1].item<Int>(), Int(1));
}

TEST_F(TorchTransferTest, ToCudaRejectsWrongDevice) {
    // Create a CPU tensor and try to extract via to_cuda
    auto tframe = atorch::from_cpu(*cpu_frame_);
    tframe.x = tframe.x.to(::torch::kCPU);  // Move x to CPU

    ATV::cuda::Frame cframe;
    EXPECT_THROW(atorch::to_cuda(tframe, cframe), std::invalid_argument);
}

TEST_F(TorchTransferTest, ToCudaRejectsWrongDtype) {
    auto tframe = atorch::from_cpu(*cpu_frame_);
    // Replace x with wrong dtype (Int instead of Real)
    auto int_opts = ::torch::TensorOptions().dtype(atorch::kIntDtype).device(::torch::kCUDA);
    tframe.x = ::torch::zeros({4}, int_opts);

    ATV::cuda::Frame cframe;
    EXPECT_THROW(atorch::to_cuda(tframe, cframe), std::invalid_argument);
}

TEST_F(TorchTransferTest, ToCudaRejectsNonContiguous) {
    auto tframe = atorch::from_cpu(*cpu_frame_);
    // Create a non-contiguous tensor via stride manipulation
    auto big = ::torch::zeros({8}, ::torch::TensorOptions()
        .dtype(atorch::kRealDtype).device(::torch::kCUDA));
    tframe.x = big.slice(0, 0, 8, 2);  // stride=2, non-contiguous, 4 elements
    EXPECT_FALSE(tframe.x.is_contiguous());

    ATV::cuda::Frame cframe;
    EXPECT_THROW(atorch::to_cuda(tframe, cframe), std::invalid_argument);
}

TEST_F(TorchTransferTest, BoxRoundTrip) {
    // Set a non-cubic box on the source frame
    cpu_frame_->box = Mat3::diagonal(Real(12.5), Real(13.0), Real(14.5));

    auto tframe = atorch::from_cpu(*cpu_frame_);

    // Verify torch::Frame box tensor has correct values
    EXPECT_TRUE(tframe.box.defined());
    EXPECT_EQ(tframe.box.numel(), 9);
    auto box_cpu = tframe.box.to(::torch::kCPU);
    EXPECT_NEAR(box_cpu[0].item<Real>(), Real(12.5), Real(1e-6)); // [0,0]
    EXPECT_NEAR(box_cpu[4].item<Real>(), Real(13.0), Real(1e-6)); // [1,1]
    EXPECT_NEAR(box_cpu[8].item<Real>(), Real(14.5), Real(1e-6)); // [2,2]
    // Off-diagonal should be zero for diagonal box
    EXPECT_NEAR(box_cpu[1].item<Real>(), Real(0), Real(1e-6));    // [0,1]
    EXPECT_NEAR(box_cpu[3].item<Real>(), Real(0), Real(1e-6));    // [1,0]

    // Round-trip: transfer back to cpu::Frame
    ATV::cpu::Frame dst(4);
    atorch::to_cpu(tframe, dst);

    EXPECT_NEAR(dst.box(0, 0), Real(12.5), Real(1e-6));
    EXPECT_NEAR(dst.box(1, 1), Real(13.0), Real(1e-6));
    EXPECT_NEAR(dst.box(2, 2), Real(14.5), Real(1e-6));
    EXPECT_NEAR(dst.box(0, 1), Real(0),    Real(1e-6));
    EXPECT_NEAR(dst.box(1, 0), Real(0),    Real(1e-6));
    EXPECT_NEAR(dst.box(2, 0), Real(0),    Real(1e-6));
}

#else

TEST(TorchTransfer, Disabled) {
    GTEST_SKIP() << "ATV_HAS_TORCH not defined";
}

#endif // ATV_HAS_TORCH
