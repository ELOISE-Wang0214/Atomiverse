/// @file neighborlist/export_pairs.cu
/// @brief Cluster flat pairs → atom-level PairList with full triclinic PBC.

#include "nlist_data.cuh"
#include "nlist_error.cuh"
#include "triclinic.cuh"
#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

#include <cuda_runtime.h>

namespace ATV::md::nlist {

using ATV::Int;
using ATV::Real;

__global__ void export_atom_pairs_kernel(
    const int* __restrict__ flat_ci,
    const int* __restrict__ flat_cj,
    const uint64_t* __restrict__ flat_mask,
    const int* __restrict__ num_flat_device,
    const int* __restrict__ original_index,
    const Real* __restrict__ pos_x,
    const Real* __restrict__ pos_y,
    const Real* __restrict__ pos_z,
    bool periodic,
    const Real* __restrict__ cell_matrix,
    const Real* __restrict__ inv_cell,
    Int* __restrict__ out_i,
    Int* __restrict__ out_j,
    Real* __restrict__ out_dx,
    Real* __restrict__ out_dy,
    Real* __restrict__ out_dz,
    Real* __restrict__ out_dist,
    Int* __restrict__ out_count,
    const Int* __restrict__ capacity_ptr) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int n_flat = *num_flat_device;
    if (idx >= n_flat) return;

    int ci = flat_ci[idx];
    int cj = flat_cj[idx];
    uint64_t mask = flat_mask[idx];
    Int cap = *capacity_ptr;

    // Batch slot reservation: count all pairs this thread will emit,
    // reserve slots with a single atomicAdd, then write sequentially.
    // Reduces global atomic contention from O(n_pairs) to O(n_threads).
    int my_count = __popcll(mask);
    Int my_base = 0;
    if (my_count > 0) my_base = atomicAdd(out_count, Int(my_count));

    int local = 0;
    while (mask != 0) {
        int bit = __ffsll(mask) - 1;
        mask &= mask - 1;

        int ai = original_index[ci * kClusterSize + bit / kClusterSize];
        int aj = original_index[cj * kClusterSize + bit % kClusterSize];

        Real dx = pos_x[aj] - pos_x[ai];
        Real dy = pos_y[aj] - pos_y[ai];
        Real dz = pos_z[aj] - pos_z[ai];

        if (periodic) {
            minimum_image(dx, dy, dz, cell_matrix, inv_cell);
        }

        Real r = sqrt(dx * dx + dy * dy + dz * dz);
        Int slot = my_base + local;
        if (slot < cap) {
            out_i[slot] = static_cast<Int>(ai);
            out_j[slot] = static_cast<Int>(aj);
            out_dx[slot] = dx;
            out_dy[slot] = dy;
            out_dz[slot] = dz;
            out_dist[slot] = r;
        }
        local++;
    }
}

void launch_export_atom_pairs(
    const int* flat_ci, const int* flat_cj, const uint64_t* flat_mask,
    const int* num_flat_device, int n_flat, const int* original_index,
    const Real* pos_x, const Real* pos_y, const Real* pos_z,
    bool periodic, const Real* cell_matrix, const Real* inv_cell,
    Int* out_i, Int* out_j, Real* out_dx, Real* out_dy, Real* out_dz,
    Real* out_dist, Int* out_count, const Int* capacity_ptr,
    cudaStream_t stream) {
    if (n_flat == 0) return;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (n_flat + block - 1) / block;
    export_atom_pairs_kernel<<<grid, block, 0, stream>>>(
        flat_ci, flat_cj, flat_mask, num_flat_device, original_index,
        pos_x, pos_y, pos_z, periodic, cell_matrix, inv_cell,
        out_i, out_j, out_dx, out_dy, out_dz, out_dist,
        out_count, capacity_ptr);
    ATV_CUDA_CHECK(cudaGetLastError());
}

} // namespace ATV::md::nlist
