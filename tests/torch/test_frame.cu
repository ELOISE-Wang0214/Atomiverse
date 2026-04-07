/// @file tests/test_torch_frame.cu
/// @brief Tests for torch::Frame type.

#include <gtest/gtest.h>

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/frame.cuh"

namespace atorch = ATV::torch;

TEST(TorchFrame, PrecisionMapping) {
    if constexpr (std::is_same_v<ATV::cuda::real_t, float>) {
        EXPECT_EQ(atorch::kRealDtype, ::at::ScalarType::Float);
    } else {
        EXPECT_EQ(atorch::kRealDtype, ::at::ScalarType::Double);
    }
    // kIntDtype is always Int (int -> ScalarType::Int)
    EXPECT_EQ(atorch::kIntDtype, ::at::ScalarType::Int);
}

TEST(TorchFrame, DefaultConstruct) {
    atorch::Frame frame;
    EXPECT_EQ(frame.n_atoms(), 0);
    EXPECT_EQ(frame.n_bonds(), 0);
    EXPECT_EQ(frame.n_pairs(), 0);
    EXPECT_FALSE(frame.atoms.x.defined());
}

TEST(TorchFrame, ManualConstruct) {
    int64_t n = 10;
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);
    auto int_opts = ::torch::TensorOptions().dtype(atorch::kIntDtype).device(::torch::kCUDA);

    atorch::Frame frame;
    frame.n_atoms_ = static_cast<int>(n);
    frame.atoms.x = ::torch::zeros({n}, real_opts);
    frame.atoms.y = ::torch::zeros({n}, real_opts);
    frame.atoms.z = ::torch::zeros({n}, real_opts);
    frame.atoms.fx = ::torch::zeros({n}, real_opts);
    frame.atoms.fy = ::torch::zeros({n}, real_opts);
    frame.atoms.fz = ::torch::zeros({n}, real_opts);
    frame.atoms.mass = ::torch::ones({n}, real_opts);
    frame.atoms.type = ::torch::zeros({n}, int_opts);
    frame.energy = ::torch::zeros({1}, real_opts);
    frame.kinetic = ::torch::zeros({1}, real_opts);
    frame.stress = ::torch::zeros({6}, real_opts);
    frame.box = ::torch::zeros({9}, real_opts);

    EXPECT_EQ(frame.n_atoms(), static_cast<int>(n));
    EXPECT_TRUE(frame.atoms.x.is_cuda());
    EXPECT_EQ(frame.atoms.x.scalar_type(), atorch::kRealDtype);
    EXPECT_EQ(frame.atoms.x.numel(), n);
    EXPECT_TRUE(frame.atoms.x.is_contiguous());
}

TEST(TorchFrame, RequiresGrad) {
    auto opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);
    auto x = ::torch::randn({5}, opts).requires_grad_(true);
    EXPECT_TRUE(x.requires_grad());

    auto int_opts = ::torch::TensorOptions().dtype(atorch::kIntDtype).device(::torch::kCUDA);
    auto type_id = ::torch::zeros({5}, int_opts);
    EXPECT_FALSE(type_id.requires_grad());
}

#else

TEST(TorchFrame, Disabled) {
    GTEST_SKIP() << "ATV_HAS_TORCH not defined";
}

#endif // ATV_HAS_TORCH
