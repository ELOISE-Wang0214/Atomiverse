/// @file torch/bonded/harmonic_bond.cu
/// @brief Differentiable harmonic bond -- forward reuses cuda kernel, backward analytical.

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/potentials/harmonic_bond.cuh"
#include "atomiverse/cuda/potentials/harmonic_bond_kernel.cuh"
#include "atomiverse/cuda/launch.cuh"

#include <torch/torch.h>
#include <c10/cuda/CUDAStream.h>
#include <cuda_runtime.h>

namespace ATV {
namespace torch {

namespace {

/// @brief Backward kernel for harmonic bond parameter gradients.
///
/// Computes dE/dK[t] = sum_bonds_of_type_t (r - r0[t])^2
///          dE/dr0[t] = sum_bonds_of_type_t -2*K[t]*(r - r0[t])
/// Also computes dE/dx for each atom.
__global__ void harmonic_bond_backward(
    const Real* __restrict__ x, const Real* __restrict__ y, const Real* __restrict__ z,
    Real* __restrict__ grad_x, Real* __restrict__ grad_y, Real* __restrict__ grad_z,
    Real* __restrict__ grad_K, Real* __restrict__ grad_r0,
    const Int* __restrict__ bi, const Int* __restrict__ bj, const Int* __restrict__ bt,
    const Real* __restrict__ K, const Real* __restrict__ r0,
    const Real* __restrict__ box,
    Real grad_energy, Int n_bonds)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_bonds) return;

    Int i = bi[idx], j = bj[idx], t = bt[idx];

    Real dx = x[j] - x[i];
    Real dy = y[j] - y[i];
    Real dz = z[j] - z[i];

    // Minimum image convention (orthorhombic)
    Real bx = box[0], by = box[4], bz = box[8];
    dx -= bx * round(dx / bx);
    dy -= by * round(dy / by);
    dz -= bz * round(dz / bz);

    Real r = sqrt(dx * dx + dy * dy + dz * dz);
    Real dr = r - r0[t];

    // dE/dx: F = -2K(r-r0) * r_hat, dE/dx = -F
    Real f_over_r = Real(2) * K[t] * dr / r;
    Real dEdx_j = f_over_r * dx * grad_energy;
    Real dEdy_j = f_over_r * dy * grad_energy;
    Real dEdz_j = f_over_r * dz * grad_energy;

    atomicAdd(&grad_x[j], dEdx_j);
    atomicAdd(&grad_y[j], dEdy_j);
    atomicAdd(&grad_z[j], dEdz_j);
    atomicAdd(&grad_x[i], -dEdx_j);
    atomicAdd(&grad_y[i], -dEdy_j);
    atomicAdd(&grad_z[i], -dEdz_j);

    // dE/dK[t] = (r - r0)^2
    atomicAdd(&grad_K[t], dr * dr * grad_energy);

    // dE/dr0[t] = -2*K*(r - r0)
    atomicAdd(&grad_r0[t], Real(-2) * K[t] * dr * grad_energy);
}

struct HarmonicBondFunction : public ::torch::autograd::Function<HarmonicBondFunction> {
    static ::torch::autograd::variable_list forward(
        ::torch::autograd::AutogradContext* ctx,
        Tensor x, Tensor y, Tensor z,
        Tensor fx, Tensor fy, Tensor fz,
        Tensor energy, Tensor box,
        Tensor bond_i, Tensor bond_j, Tensor bond_type,
        Tensor K, Tensor r0,
        int64_t n_bonds)
    {
        auto stream = ::c10::cuda::getCurrentCUDAStream().stream();
        int block = cuda::kDefaultBlock;
        int grid = cuda::grid_size(static_cast<int>(n_bonds));

        cuda::kernels::harmonic_bond_k<<<grid, block, block * sizeof(Real), stream>>>(
            x.data_ptr<Real>(), y.data_ptr<Real>(), z.data_ptr<Real>(),
            fx.data_ptr<Real>(), fy.data_ptr<Real>(), fz.data_ptr<Real>(),
            energy.data_ptr<Real>(),
            bond_i.data_ptr<Int>(), bond_j.data_ptr<Int>(), bond_type.data_ptr<Int>(),
            K.data_ptr<Real>(), r0.data_ptr<Real>(),
            box.data_ptr<Real>(), static_cast<Int>(n_bonds));

        ctx->save_for_backward({x, y, z, box, bond_i, bond_j, bond_type, K, r0});
        ctx->saved_data["n_bonds"] = n_bonds;

        return {energy};
    }

    static ::torch::autograd::variable_list backward(
        ::torch::autograd::AutogradContext* ctx,
        ::torch::autograd::variable_list grad_outputs)
    {
        auto saved = ctx->get_saved_variables();
        auto x = saved[0], y = saved[1], z = saved[2], box = saved[3];
        auto bond_i = saved[4], bond_j = saved[5], bond_type = saved[6];
        auto K = saved[7], r0 = saved[8];

        int64_t n_bonds = ctx->saved_data["n_bonds"].toInt();
        Real ge = grad_outputs[0].item<Real>();

        int64_t n_atoms = x.numel();
        int64_t n_types = K.numel();

        auto opts = ::torch::TensorOptions().dtype(kRealDtype).device(x.device());
        auto grad_x = ::torch::zeros({n_atoms}, opts);
        auto grad_y = ::torch::zeros({n_atoms}, opts);
        auto grad_z = ::torch::zeros({n_atoms}, opts);
        auto grad_K = ::torch::zeros({n_types}, opts);
        auto grad_r0 = ::torch::zeros({n_types}, opts);

        auto stream = ::c10::cuda::getCurrentCUDAStream().stream();
        int block = cuda::kDefaultBlock;
        int grid = cuda::grid_size(static_cast<int>(n_bonds));

        harmonic_bond_backward<<<grid, block, 0, stream>>>(
            x.data_ptr<Real>(), y.data_ptr<Real>(), z.data_ptr<Real>(),
            grad_x.data_ptr<Real>(), grad_y.data_ptr<Real>(), grad_z.data_ptr<Real>(),
            grad_K.data_ptr<Real>(), grad_r0.data_ptr<Real>(),
            bond_i.data_ptr<Int>(), bond_j.data_ptr<Int>(), bond_type.data_ptr<Int>(),
            K.data_ptr<Real>(), r0.data_ptr<Real>(),
            box.data_ptr<Real>(), ge, static_cast<Int>(n_bonds));

        // Return grads: x, y, z, fx, fy, fz, energy, box,
        //               bond_i, bond_j, bond_type, K, r0, n_bonds
        return {grad_x, grad_y, grad_z,
                Tensor(), Tensor(), Tensor(), Tensor(), Tensor(), // fx,fy,fz,energy,box
                Tensor(), Tensor(), Tensor(), // bond_i, bond_j, bond_type
                grad_K, grad_r0, Tensor()};  // K, r0, n_bonds
    }
};

} // anonymous namespace

Tensor harmonic_bond(Frame& frame, const Tensor& K, const Tensor& r0) {
    Int nb = frame.n_bonds();
    if (nb == 0) return frame.energy;

    frame.energy.zero_();

    auto result = HarmonicBondFunction::apply(
        frame.x, frame.y, frame.z,
        frame.fx, frame.fy, frame.fz,
        frame.energy, frame.box,
        frame.bond_i, frame.bond_j, frame.bond_type,
        K, r0,
        static_cast<int64_t>(nb));

    return result[0];
}

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
