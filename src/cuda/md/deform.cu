/// @file md/deform.cu
/// @brief Deform -- uniaxial constant-strain-rate deformation.

#include "atomiverse/cuda/md/deform.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/facade/stage.h"
#include <cuda_runtime.h>

namespace ATV {
namespace md {

Deform::Deform(int axis, Real strain_rate, Real dt)
    : axis_(axis), strain_rate_(strain_rate), dt_(dt) {}

std::vector<DefaultBinding> Deform::default_bindings() const {
    return {{static_cast<int>(Phase::initial_integrate), -100}};
}

namespace {

__global__ void deform_k(
    Real* __restrict__ pos,       // x, y, or z depending on axis
    Real* __restrict__ box,       // 3x3 row-major
    int axis, Real scale_factor, Int n)
{
    Int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    // Remap atom position proportionally
    pos[idx] *= scale_factor;

    // Thread 0 updates the box diagonal element
    if (idx == 0) {
        box[axis * 3 + axis] *= scale_factor;
    }
}

} // anonymous namespace

void Deform::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    Int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    Real scale = Real(1) + strain_rate_ * dt_;

    Real* pos = nullptr;
    switch (axis_) {
        case 0: pos = f.x; break;
        case 1: pos = f.y; break;
        case 2: pos = f.z; break;
    }

    deform_k<<<grid, block, 0, stream>>>(pos, f.box, axis_, scale, n);
}

} // namespace md
} // namespace ATV
