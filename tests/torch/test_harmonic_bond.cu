/// @file tests/test_torch_harmonic_bond.cu
/// @brief Tests for differentiable HarmonicBond operator in the torch backend.

#include <gtest/gtest.h>

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/potentials/harmonic_bond.cuh"

#include <cuda_runtime.h>
#include <cmath>

using ATV::Real;
using ATV::Int;
namespace atorch = ATV::torch;

namespace {

Real analytical_bond_energy(Real r, Real K, Real r0) {
    Real dr = r - r0;
    return K * dr * dr;
}

atorch::Frame make_bond_frame(Real x0, Real x1) {
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);
    auto int_opts = ::torch::TensorOptions().dtype(atorch::kIntDtype).device(::torch::kCUDA);

    atorch::Frame frame;
    frame.n_atoms_ = 2;
    frame.n_bonds_ = 1;
    frame.atoms.x = ::torch::tensor({x0, x1}, real_opts);
    frame.atoms.y = ::torch::zeros({2}, real_opts);
    frame.atoms.z = ::torch::zeros({2}, real_opts);
    frame.atoms.fx = ::torch::zeros({2}, real_opts);
    frame.atoms.fy = ::torch::zeros({2}, real_opts);
    frame.atoms.fz = ::torch::zeros({2}, real_opts);
    frame.energy = ::torch::zeros({1}, real_opts);
    frame.bonds.atomi = ::torch::tensor({Int(0)}, int_opts);
    frame.bonds.atomj = ::torch::tensor({Int(1)}, int_opts);
    frame.bonds.type = ::torch::zeros({1}, int_opts);

    Real box_data[9] = {100, 0, 0, 0, 100, 0, 0, 0, 100};
    auto cpu_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype);
    frame.box = ::torch::from_blob(box_data, {9}, cpu_opts).to(::torch::kCUDA);
    return frame;
}

} // anonymous namespace

TEST(TorchHarmonicBond, SingleBondAnalytical) {
    Real K_val = Real(100.0), r0_val = Real(1.54), x1 = Real(1.8);
    auto frame = make_bond_frame(Real(0), x1);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    auto energy = atorch::harmonic_bond(frame,
        ::torch::tensor({K_val}, real_opts),
        ::torch::tensor({r0_val}, real_opts));
    cudaDeviceSynchronize();

    Real computed = energy.to(::torch::kCPU).item<Real>();
    Real expected = analytical_bond_energy(x1, K_val, r0_val);
    EXPECT_NEAR(computed, expected, Real(1e-4));
}

TEST(TorchHarmonicBond, GradcheckK) {
    Real K_val = Real(100.0), r0_val = Real(1.54), x1 = Real(1.8);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    Real h = Real(0.01);
    Real numerical = (analytical_bond_energy(x1, K_val + h, r0_val)
                    - analytical_bond_energy(x1, K_val - h, r0_val)) / (Real(2) * h);

    auto frame = make_bond_frame(Real(0), x1);
    auto K = ::torch::tensor({K_val}, real_opts).requires_grad_(true);
    auto r0 = ::torch::tensor({r0_val}, real_opts);

    auto energy = atorch::harmonic_bond(frame, K, r0);
    energy.backward();
    cudaDeviceSynchronize();

    Real autograd = K.grad().to(::torch::kCPU).item<Real>();
    EXPECT_NEAR(autograd, numerical, Real(1e-2));
}

TEST(TorchHarmonicBond, GradcheckR0) {
    Real K_val = Real(100.0), r0_val = Real(1.54), x1 = Real(1.8);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    Real h = Real(1e-4);
    Real numerical = (analytical_bond_energy(x1, K_val, r0_val + h)
                    - analytical_bond_energy(x1, K_val, r0_val - h)) / (Real(2) * h);

    auto frame = make_bond_frame(Real(0), x1);
    auto K = ::torch::tensor({K_val}, real_opts);
    auto r0 = ::torch::tensor({r0_val}, real_opts).requires_grad_(true);

    auto energy = atorch::harmonic_bond(frame, K, r0);
    energy.backward();
    cudaDeviceSynchronize();

    Real autograd = r0.grad().to(::torch::kCPU).item<Real>();
    EXPECT_NEAR(autograd, numerical, Real(1e-2));
}

TEST(TorchHarmonicBond, ZeroBonds) {
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    atorch::Frame frame;
    frame.n_atoms_ = 2;
    frame.n_bonds_ = 0;
    frame.atoms.x = ::torch::zeros({2}, real_opts);
    frame.atoms.y = ::torch::zeros({2}, real_opts);
    frame.atoms.z = ::torch::zeros({2}, real_opts);
    frame.atoms.fx = ::torch::zeros({2}, real_opts);
    frame.atoms.fy = ::torch::zeros({2}, real_opts);
    frame.atoms.fz = ::torch::zeros({2}, real_opts);
    frame.energy = ::torch::zeros({1}, real_opts);

    auto energy = atorch::harmonic_bond(frame,
        ::torch::tensor({Real(100.0)}, real_opts),
        ::torch::tensor({Real(1.54)}, real_opts));
    cudaDeviceSynchronize();

    EXPECT_NEAR(energy.to(::torch::kCPU).item<Real>(), Real(0), Real(1e-6));
}

TEST(TorchHarmonicBond, ForceSymmetry) {
    // Newton's 3rd law: for a stretched bond, fx[0] = -fx[1]
    Real K_val = Real(100.0), r0_val = Real(1.54);
    Real x0 = Real(0), x1 = Real(1.8);
    auto frame = make_bond_frame(x0, x1);
    auto real_opts = ::torch::TensorOptions().dtype(atorch::kRealDtype).device(::torch::kCUDA);

    auto energy = atorch::harmonic_bond(frame,
        ::torch::tensor({K_val}, real_opts),
        ::torch::tensor({r0_val}, real_opts));
    cudaDeviceSynchronize();

    auto fx_cpu = frame.atoms.fx.to(::torch::kCPU);
    Real fx0 = fx_cpu[0].item<Real>();
    Real fx1 = fx_cpu[1].item<Real>();

    // Forces must be opposite: fx[0] = -fx[1]
    EXPECT_NEAR(fx0, -fx1, Real(1e-6));

    // Forces should be non-zero (bond is stretched beyond r0)
    EXPECT_GT(std::abs(fx0), Real(1e-6));

    // fy and fz should be zero (bond is along x only)
    auto fy_cpu = frame.atoms.fy.to(::torch::kCPU);
    auto fz_cpu = frame.atoms.fz.to(::torch::kCPU);
    EXPECT_NEAR(fy_cpu[0].item<Real>(), Real(0), Real(1e-6));
    EXPECT_NEAR(fy_cpu[1].item<Real>(), Real(0), Real(1e-6));
    EXPECT_NEAR(fz_cpu[0].item<Real>(), Real(0), Real(1e-6));
    EXPECT_NEAR(fz_cpu[1].item<Real>(), Real(0), Real(1e-6));
}

#else

TEST(TorchHarmonicBond, Disabled) {
    GTEST_SKIP() << "ATV_HAS_TORCH not defined";
}

#endif // ATV_HAS_TORCH
