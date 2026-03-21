/// @file torch/md/lj.cu
/// @brief Differentiable LJ operator -- forward reuses cuda kernel, backward is analytical.

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/potentials/lj.cuh"
#include "atomiverse/cuda/potentials/lj_kernel.cuh"
#include "atomiverse/cuda/launch.cuh"

#include <torch/torch.h>
#include <c10/cuda/CUDAStream.h>
#include <cuda_runtime.h>

namespace ATV {
namespace torch {

namespace {

/// @brief Backward kernel: compute dE/d(epsilon), dE/d(sigma) from pair distances.
///
/// dE/d(epsilon) = sum_pairs 4 * (sr12 - sr6) - v_shift/epsilon
/// dE/d(sigma) = sum_pairs 4*epsilon * (-12*sigma^11/r^12 + 6*sigma^5/r^6)
__global__ void lj_param_grad(
    Real* __restrict__ grad_eps_out,
    Real* __restrict__ grad_sig_out,
    const Real* __restrict__ pdist,
    const Int* __restrict__ n_pairs_ptr,
    Real epsilon, Real sigma, Real cutoff, Real v_shift)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    Int n_pairs = *n_pairs_ptr;

    Real de = Real(0);
    Real ds = Real(0);
    if (idx < n_pairs) {
        Real r = pdist[idx];
        if (r < cutoff && r >= cuda::kernels::kMinPairDist) {
            Real inv_r = Real(1) / r;
            Real sr = sigma * inv_r;
            Real sr2 = sr * sr;
            Real sr6 = sr2 * sr2 * sr2;
            Real sr12 = sr6 * sr6;

            // dV/d(epsilon) = 4*(sr12 - sr6) - v_shift/epsilon
            de = Real(4) * (sr12 - sr6) - v_shift / epsilon;

            // dV/d(sigma): d/dsigma [4*eps*(sigma/r)^12 - 4*eps*(sigma/r)^6]
            //   = 4*eps * (12*sigma^11/r^12 - 6*sigma^5/r^6)
            Real inv_sig = Real(1) / sigma;
            ds = Real(4) * epsilon * (Real(12) * sr12 * inv_sig - Real(6) * sr6 * inv_sig);
        }
    }

    // Block reduction for grad_eps
    extern __shared__ Real sdata[];
    Real* s_eps = sdata;
    Real* s_sig = sdata + blockDim.x;
    s_eps[threadIdx.x] = de;
    s_sig[threadIdx.x] = ds;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s_eps[threadIdx.x] += s_eps[threadIdx.x + s];
            s_sig[threadIdx.x] += s_sig[threadIdx.x + s];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        if (s_eps[0] != Real(0)) atomicAdd(grad_eps_out, s_eps[0]);
        if (s_sig[0] != Real(0)) atomicAdd(grad_sig_out, s_sig[0]);
    }
}

/// @brief autograd::Function for LJ forward+backward.
struct LJFunction : public ::torch::autograd::Function<LJFunction> {
    static ::torch::autograd::variable_list forward(
        ::torch::autograd::AutogradContext* ctx,
        Tensor x, Tensor y, Tensor z,
        Tensor fx, Tensor fy, Tensor fz,
        Tensor energy,
        Tensor pair_i, Tensor pair_j,
        Tensor pair_dx, Tensor pair_dy, Tensor pair_dz,
        Tensor pair_dist, Tensor n_pairs_t,
        Tensor epsilon, Tensor sigma,
        double cutoff, int64_t n_atoms, int64_t max_pairs)
    {
        Real eps_val = epsilon.item<Real>();
        Real sig_val = sigma.item<Real>();
        Real cut_val = static_cast<Real>(cutoff);
        Real v_shift = cuda::kernels::lj_shift(eps_val, sig_val, cut_val);

        auto stream = ::c10::cuda::getCurrentCUDAStream().stream();
        int block = cuda::kDefaultBlock;
        int grid = cuda::grid_size(static_cast<int>(max_pairs));

        cuda::kernels::lj_force<<<grid, block, block * sizeof(Real), stream>>>(
            fx.data_ptr<Real>(), fy.data_ptr<Real>(), fz.data_ptr<Real>(),
            energy.data_ptr<Real>(),
            pair_i.data_ptr<Int>(), pair_j.data_ptr<Int>(),
            pair_dx.data_ptr<Real>(), pair_dy.data_ptr<Real>(),
            pair_dz.data_ptr<Real>(), pair_dist.data_ptr<Real>(),
            eps_val, sig_val, cut_val, v_shift,
            n_pairs_t.data_ptr<Int>());

        // Save for backward
        ctx->save_for_backward({pair_i, pair_j, pair_dx, pair_dy, pair_dz,
                                pair_dist, n_pairs_t, epsilon, sigma,
                                fx, fy, fz});
        ctx->saved_data["cutoff"] = cutoff;
        ctx->saved_data["n_atoms"] = n_atoms;
        ctx->saved_data["max_pairs"] = max_pairs;
        ctx->saved_data["v_shift"] = static_cast<double>(v_shift);

        return {energy};
    }

    static ::torch::autograd::variable_list backward(
        ::torch::autograd::AutogradContext* ctx,
        ::torch::autograd::variable_list grad_outputs)
    {
        auto saved = ctx->get_saved_variables();
        auto pair_i = saved[0], pair_j = saved[1];
        auto pair_dx = saved[2], pair_dy = saved[3], pair_dz = saved[4];
        auto pair_dist = saved[5], n_pairs_t = saved[6];
        auto epsilon = saved[7], sigma = saved[8];
        auto fx = saved[9], fy = saved[10], fz = saved[11];

        auto grad_energy = grad_outputs[0]; // scalar
        Real ge = grad_energy.item<Real>();

        double cutoff = ctx->saved_data["cutoff"].toDouble();
        int64_t n_atoms = ctx->saved_data["n_atoms"].toInt();
        int64_t max_pairs = ctx->saved_data["max_pairs"].toInt();
        Real v_shift = static_cast<Real>(ctx->saved_data["v_shift"].toDouble());

        // dE/dx = -fx (force is negative gradient of energy)
        // These are already computed in forward; scale by grad_energy
        Tensor grad_x = Tensor();
        Tensor grad_y = Tensor();
        Tensor grad_z = Tensor();
        if (ge != Real(0)) {
            grad_x = -fx * ge;
            grad_y = -fy * ge;
            grad_z = -fz * ge;
        }

        // dE/d(epsilon) and dE/d(sigma) via custom kernel
        Tensor grad_eps = Tensor();
        Tensor grad_sig = Tensor();

        Real eps_val = epsilon.item<Real>();
        Real sig_val = sigma.item<Real>();

        if (epsilon.requires_grad() || sigma.requires_grad()) {
            auto opts = ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
            auto d_grad_eps = ::torch::zeros({1}, opts);
            auto d_grad_sig = ::torch::zeros({1}, opts);

            auto stream = ::c10::cuda::getCurrentCUDAStream().stream();
            int block = cuda::kDefaultBlock;
            int grid = cuda::grid_size(static_cast<int>(max_pairs));
            lj_param_grad<<<grid, block, 2 * block * sizeof(Real), stream>>>(
                d_grad_eps.data_ptr<Real>(),
                d_grad_sig.data_ptr<Real>(),
                pair_dist.data_ptr<Real>(),
                n_pairs_t.data_ptr<Int>(),
                eps_val, sig_val, static_cast<Real>(cutoff), v_shift);

            if (epsilon.requires_grad())
                grad_eps = d_grad_eps * ge;
            if (sigma.requires_grad())
                grad_sig = d_grad_sig * ge;
        }

        // Return grads in same order as forward args:
        // x, y, z, fx, fy, fz, energy, pair_i, pair_j,
        // pair_dx, pair_dy, pair_dz, pair_dist, n_pairs_t,
        // epsilon, sigma, cutoff, n_atoms, max_pairs
        return {grad_x, grad_y, grad_z,
                Tensor(), Tensor(), Tensor(), Tensor(), // fx, fy, fz, energy
                Tensor(), Tensor(),   // pair_i, pair_j
                Tensor(), Tensor(), Tensor(), Tensor(), Tensor(), // pair_dx..n_pairs
                grad_eps, grad_sig,
                Tensor(), Tensor(), Tensor()}; // cutoff, n_atoms, max_pairs
    }
};

} // anonymous namespace

Tensor lj(Frame& frame, const Tensor& epsilon, const Tensor& sigma, Real cutoff) {
    Int na = frame.n_atoms();
    Int np = frame.n_pairs_;

    // Estimate max pairs for grid sizing
    int64_t max_pairs = np > 0 ? static_cast<int64_t>(np)
                               : static_cast<int64_t>(na) * 50;

    // n_pairs as device tensor (for kernel)
    Tensor n_pairs_t;
    if (frame.pair_i.defined() && frame.pair_i.numel() > 0) {
        // If n_pairs_ is set, create a device scalar
        auto opts = ::torch::TensorOptions().dtype(kIntDtype).device(::torch::kCUDA);
        n_pairs_t = ::torch::full({1}, static_cast<int64_t>(np), opts);
    } else {
        auto opts = ::torch::TensorOptions().dtype(kIntDtype).device(::torch::kCUDA);
        n_pairs_t = ::torch::zeros({1}, opts);
    }

    // Zero energy before accumulation
    frame.energy.zero_();

    auto result = LJFunction::apply(
        frame.x, frame.y, frame.z,
        frame.fx, frame.fy, frame.fz,
        frame.energy,
        frame.pair_i, frame.pair_j,
        frame.pair_dx, frame.pair_dy, frame.pair_dz,
        frame.pair_dist, n_pairs_t,
        epsilon, sigma,
        static_cast<double>(cutoff), static_cast<int64_t>(na), max_pairs);

    return result[0];
}

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
