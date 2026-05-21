/// @file torch/neighbor_list.cu
/// @brief NeighborList for torch backend -- holds ATV::NeighborList.

#ifdef ATV_HAS_TORCH

#include "atomiverse/cpu/transfer.h"
#include "atomiverse/cuda/component.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/transfer.cuh"
#include "atomiverse/torch/neighbor_list.cuh"
#include "atomiverse/torch/transfer.cuh"

#include <c10/cuda/CUDAStream.h>
#include <cuda_runtime.h>
#include <torch/torch.h>

namespace ATV {
namespace torch {

// ── NeighborList class ─────────────────────────────────────────────────

NeighborList::NeighborList(double cutoff, double skin)
    : cutoff_(cutoff), skin_(skin) {}

NeighborList::~NeighborList() {
    if (built_ && nl_) {
        DestroyContext dctx{ctx_, cframe_};
        nl_->destroy(dctx);
        // Free the cuda::Frame arrays allocated by cpu::to_cuda
        // (but NOT the pairs -- already freed by nl_->destroy)
        auto saved_pairs = cframe_.pairs;
        cframe_.pairs = {}; // prevent double-free
        cuda::destroy(cframe_);
        cframe_.pairs = saved_pairs; // already freed
    }
}

NeighborList::NeighborList(NeighborList&&) noexcept = default;
NeighborList& NeighborList::operator=(NeighborList&&) noexcept = default;

void NeighborList::build(Frame& frame) {
    // Set up CUDA context using the current Torch stream
    auto stream = ::c10::cuda::getCurrentCUDAStream();
    ctx_.stream = static_cast<void*>(stream.stream());
    ctx_.defer_stream = ctx_.stream;

    // Transfer torch::Frame -> cpu::Frame -> cuda::Frame (core allocator)
    // This ensures the cuda::Frame uses cudaMalloc, which the NeighborList
    // pipeline expects.
    cpu::Frame cpu_frame(frame.n_atoms());
    to_cpu(frame, cpu_frame);
    cframe_ = cpu::to_cuda(cpu_frame);

    // Create and build the CUDA NeighborList
    nl_ = std::make_unique<ATV::NeighborList>(cutoff_, skin_);
    BuildContext bctx{ctx_, cframe_};
    nl_->build(bctx);

    built_ = true;
}

void NeighborList::update(Frame& frame) {
    if (!built_) {
        build(frame);
    }

    // Copy current positions from torch tensors to cuda::Frame arrays
    // (the cuda::Frame owns its own cudaMalloc memory)
    cudaStream_t stream = static_cast<cudaStream_t>(ctx_.stream);

    int n = frame.n_atoms();
    cudaMemcpyAsync(cframe_.atoms.x, frame.atoms.x.data_ptr<cuda::real_t>(),
                    n * sizeof(cuda::real_t), cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(cframe_.atoms.y, frame.atoms.y.data_ptr<cuda::real_t>(),
                    n * sizeof(cuda::real_t), cudaMemcpyDeviceToDevice, stream);
    cudaMemcpyAsync(cframe_.atoms.z, frame.atoms.z.data_ptr<cuda::real_t>(),
                    n * sizeof(cuda::real_t), cudaMemcpyDeviceToDevice, stream);
    if (frame.box.defined() && cframe_.box) {
        cudaMemcpyAsync(cframe_.box, frame.box.data_ptr<cuda::prec_t>(),
                        9 * sizeof(cuda::prec_t), cudaMemcpyDeviceToDevice,
                        stream);
    }

    // Launch pair search
    LaunchContext lctx{ctx_, cframe_, 0, 0};
    nl_->launch(lctx);

    cudaStreamSynchronize(stream);

    // Read n_pairs from device
    int h_n_pairs = 0;
    if (cframe_.pairs.n_pairs) {
        cudaMemcpy(&h_n_pairs, cframe_.pairs.n_pairs, sizeof(int),
                   cudaMemcpyDeviceToHost);
    }
    frame.n_pairs_ = h_n_pairs;

    // Wrap pair list as self-contained Tensors (cloned)
    if (h_n_pairs > 0) {
        auto real_opts =
            ::torch::TensorOptions().dtype(kRealDtype).device(::torch::kCUDA);
        auto int_opts =
            ::torch::TensorOptions().dtype(kIntDtype).device(::torch::kCUDA);
        int64_t np = static_cast<int64_t>(h_n_pairs);

        frame.pair_i =
            ::torch::from_blob(cframe_.pairs.i, {np}, int_opts).clone();
        frame.pair_j =
            ::torch::from_blob(cframe_.pairs.j, {np}, int_opts).clone();
        frame.pair_dx =
            ::torch::from_blob(cframe_.pairs.dx, {np}, real_opts).clone();
        frame.pair_dy =
            ::torch::from_blob(cframe_.pairs.dy, {np}, real_opts).clone();
        frame.pair_dz =
            ::torch::from_blob(cframe_.pairs.dz, {np}, real_opts).clone();
        frame.pair_dist =
            ::torch::from_blob(cframe_.pairs.dist, {np}, real_opts).clone();
    } else {
        frame.pair_i = Tensor();
        frame.pair_j = Tensor();
        frame.pair_dx = Tensor();
        frame.pair_dy = Tensor();
        frame.pair_dz = Tensor();
        frame.pair_dist = Tensor();
    }
}

// ── Convenience free function ──────────────────────────────────────────

void neighbor_list(Frame& frame, double cutoff, double skin) {
    NeighborList nl(cutoff, skin);
    nl.build(frame);
    nl.update(frame);
}

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
