/// @file md/deform.cu
/// @brief Deform -- uniaxial constant-strain-rate deformation.

#include "atomiverse/cuda/md/deform.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/facade/stage.h"
#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"

namespace ATV {

using cuda::real_t;
using cuda::work_t;
using cuda::prec_t;

namespace md {

Deform::Deform(int axis, double strain_rate, double dt)
    : axis_(axis), strain_rate_(strain_rate), dt_(dt) {}

std::vector<DefaultBinding> Deform::default_bindings() const {
    return {{static_cast<int>(Phase::initial_integrate), -100}};
}

namespace {

__global__ void deform_k(
    real_t* __restrict__ pos,       // x, y, or z depending on axis
    prec_t* __restrict__ box,       // 3x3 row-major
    int axis, work_t scale_factor, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    // Remap atom position proportionally
    pos[idx] = static_cast<real_t>(static_cast<work_t>(pos[idx]) * scale_factor);

    // Thread 0 updates the box diagonal element
    if (idx == 0) {
        int box_idx = axis * 3 + axis;
        box[box_idx] = static_cast<prec_t>(static_cast<work_t>(box[box_idx]) * scale_factor);
    }
}

} // anonymous namespace

void Deform::launch(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int block = cuda::kDefaultBlock;
    int grid = cuda::grid_size(static_cast<int>(n));

    work_t scale = work_t(1) + strain_rate_ * dt_;

    real_t* pos = nullptr;
    switch (axis_) {
        case 0: pos = f.atoms.x; break;
        case 1: pos = f.atoms.y; break;
        case 2: pos = f.atoms.z; break;
    }

    deform_k<<<grid, block, 0, stream>>>(pos, f.box, axis_, scale, n);
}

} // namespace md
} // namespace ATV
