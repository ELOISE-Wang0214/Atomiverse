/// @file fused_search.cu
/// @brief CSR↔Flat conversion, displacement tracking, and flat-mask refresh kernels.
///
/// Supports the GPU-Resident Fused Neighbor List (GRFNL) pipeline.
/// The build uses the proven legacy search pipeline (coarse count → fill →
/// compact) for full warp utilization, then expands CSR to flat pair format
/// for the smart prune path.

#include "nlist_impl.cuh"
#include "triclinic.cuh"

#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

#include <cstdint>

namespace ATV::md::nlist {

// ---------------------------------------------------------------------------
// Kernel: CSR → Flat expansion (one warp per cluster, parallel copy)
// ---------------------------------------------------------------------------

__global__ void csr_to_flat(const int* __restrict__ row_ptr,
                            const int* __restrict__ col_idx,
                            const uint64_t* __restrict__ mask_in,
                            int num_clusters, int* __restrict__ flat_ci,
                            int* __restrict__ flat_cj,
                            uint64_t* __restrict__ flat_mask) {

    int warp_id = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int lane = threadIdx.x % 32;

    if (warp_id >= num_clusters)
        return;

    int ci = warp_id;
    int begin = row_ptr[ci];
    int end = row_ptr[ci + 1];

    for (int p = begin + lane; p < end; p += 32) {
        flat_ci[p] = ci;
        flat_cj[p] = col_idx[p];
        flat_mask[p] = mask_in[p];
    }
}

void launch_csr_to_flat(const int* row_ptr, const int* col_idx,
                        const uint64_t* interaction_mask, int num_clusters,
                        int* flat_ci, int* flat_cj, uint64_t* flat_mask,
                        cudaStream_t stream) {
    if (num_clusters == 0)
        return;

    int threads_needed = num_clusters * 32;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (threads_needed + block - 1) / block;

    csr_to_flat<<<grid, block, 0, stream>>>(row_ptr, col_idx, interaction_mask,
                                            num_clusters, flat_ci, flat_cj,
                                            flat_mask);
    ATV_CUDA_CHECK(cudaGetLastError());
}

// ---------------------------------------------------------------------------
// Kernel: Per-cluster displacement + dirty flags + rebuild flag
// ---------------------------------------------------------------------------

template <typename work_t>
__global__ void cluster_displacement(
    const real_t* __restrict__ cluster_x, const real_t* __restrict__ cluster_y,
    const real_t* __restrict__ cluster_z,
    const real_t* __restrict__ ref_cluster_x,
    const real_t* __restrict__ ref_cluster_y,
    const real_t* __restrict__ ref_cluster_z,
    const int* __restrict__ cluster_count, int num_clusters,
    work_t dirty_threshold_sq, real_t* __restrict__ disp_sq_out,
    int* __restrict__ dirty_flags, int* __restrict__ rebuild_flag,
    work_t half_skin_sq) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int lane = tid % kClusterSize;
    int ci = tid / kClusterSize;

    if (ci >= num_clusters)
        return;

    int count = cluster_count[ci];
    int idx = ci * kClusterSize + lane;

    work_t dx =
        (lane < count) ? (cluster_x[idx] - ref_cluster_x[idx]) : work_t(0);
    work_t dy =
        (lane < count) ? (cluster_y[idx] - ref_cluster_y[idx]) : work_t(0);
    work_t dz =
        (lane < count) ? (cluster_z[idx] - ref_cluster_z[idx]) : work_t(0);
    work_t d2 = dx * dx + dy * dy + dz * dz;

    // Warp-level max reduction within the kClusterSize-thread group
    unsigned mask =
        0xFFu << (((threadIdx.x % 32) / kClusterSize) * kClusterSize);
    for (int offset = kClusterSize / 2; offset >= 1; offset >>= 1) {
        work_t tmp = __shfl_down_sync(mask, d2, offset);
        d2 = fmax(d2, tmp);
    }

    if (lane == 0) {
        disp_sq_out[ci] = d2;
        dirty_flags[ci] = (d2 > dirty_threshold_sq) ? 1 : 0;
        if (d2 > half_skin_sq) {
            atomicMax(rebuild_flag, 1);
        }
    }
}

template <typename work_t>
void launch_cluster_displacement(
    const real_t* cluster_x, const real_t* cluster_y, const real_t* cluster_z,
    const real_t* ref_cluster_x, const real_t* ref_cluster_y,
    const real_t* ref_cluster_z, const int* cluster_count, int num_clusters,
    work_t dirty_threshold_sq, work_t* cluster_disp_sq, int* dirty_flags,
    int* rebuild_flag, work_t half_skin_sq, cudaStream_t stream) {
    if (num_clusters == 0)
        return;

    ATV_CUDA_CHECK(cudaMemsetAsync(rebuild_flag, 0, sizeof(int), stream));

    int threads_needed = num_clusters * kClusterSize;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (threads_needed + block - 1) / block;

    cluster_displacement<work_t><<<grid, block, 0, stream>>>(
        cluster_x, cluster_y, cluster_z, ref_cluster_x, ref_cluster_y,
        ref_cluster_z, cluster_count, num_clusters, dirty_threshold_sq,
        cluster_disp_sq, dirty_flags, rebuild_flag, half_skin_sq);
    ATV_CUDA_CHECK(cudaGetLastError());
}

template void launch_cluster_displacement<float>(const float*, const float*,
                                                 const float*, const float*,
                                                 const float*, const float*,
                                                 const int*, int, float, float*,
                                                 int*, int*, float,
                                                 cudaStream_t);
template void launch_cluster_displacement<double>(const double*, const double*,
                                                  const double*, const double*,
                                                  const double*, const double*,
                                                  const int*, int, double,
                                                  double*, int*, int*, double,
                                                  cudaStream_t);

// ---------------------------------------------------------------------------
// Kernel: Refresh atom-pair masks on the flat pair list.
//
// Correctness note:
// The mask must be recomputed from current coordinates instead of being
// monotonically intersected with the previous mask. Otherwise a pair that
// leaves the cutoff and later re-enters before the next rebuild can be lost
// permanently.
// ---------------------------------------------------------------------------

template <typename work_t>
__global__ __launch_bounds__(ATV::cuda::kDefaultBlock) void smart_rolling_prune(
    const real_t* __restrict__ cluster_x, const real_t* __restrict__ cluster_y,
    const real_t* __restrict__ cluster_z, const int* __restrict__ flat_ci,
    const int* __restrict__ flat_cj,
    const int* __restrict__ num_flat_pairs_device,
    const int* __restrict__ dirty_flags, bool periodic, work_t cutoff_sq,
    const real_t* __restrict__ cell_matrix, const real_t* __restrict__ inv_cell,
    uint64_t* __restrict__ flat_mask, int prune_step, int prune_period) {

    // kClusterSize threads per pair (8 threads compute the mask)
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int pair_idx = tid / kClusterSize;
    int i_lane = tid % kClusterSize;
    int num_flat_pairs = *num_flat_pairs_device;

    if (pair_idx >= num_flat_pairs)
        return;

    int ci = flat_ci[pair_idx];
    int cj = flat_cj[pair_idx];

    work_t xi = cluster_x[ci * kClusterSize + i_lane];
    work_t yi = cluster_y[ci * kClusterSize + i_lane];
    work_t zi = cluster_z[ci * kClusterSize + i_lane];

    work_t jx = cluster_x[cj * kClusterSize + i_lane];
    work_t jy = cluster_y[cj * kClusterSize + i_lane];
    work_t jz = cluster_z[cj * kClusterSize + i_lane];

    uint8_t my_bits = 0;
    int group_in_warp = (threadIdx.x % 32) / kClusterSize;
    unsigned submask = 0xFFu << (group_in_warp * kClusterSize);
    int group_base = group_in_warp * kClusterSize;

    for (int j = 0; j < kClusterSize; ++j) {
        work_t jx_j = __shfl_sync(submask, jx, group_base + j);
        work_t jy_j = __shfl_sync(submask, jy, group_base + j);
        work_t jz_j = __shfl_sync(submask, jz, group_base + j);

        work_t dx = jx_j - xi;
        work_t dy = jy_j - yi;
        work_t dz = jz_j - zi;

        if (periodic) {
            minimum_image(dx, dy, dz, cell_matrix, inv_cell);
        }

        work_t r2 = dx * dx + dy * dy + dz * dz;
        if (r2 < cutoff_sq) {
            my_bits |= (1u << j);
        }
    }

    // Half-Newton exclusion for diagonal pairs
    if (ci == cj) {
        uint8_t upper = static_cast<uint8_t>(0xFFu << (i_lane + 1));
        my_bits &= upper;
    }

    uint64_t new_mask = 0;
    for (int k = 0; k < kClusterSize; ++k) {
        uint8_t other = static_cast<uint8_t>(
            __shfl_sync(submask, static_cast<int>(my_bits), group_base + k));
        new_mask |= (static_cast<uint64_t>(other) << (k * kClusterSize));
    }

    if (i_lane == 0) {
        flat_mask[pair_idx] = new_mask;
    }
}

template <typename work_t>
void launch_smart_rolling_prune(
    const real_t* cluster_x, const real_t* cluster_y, const real_t* cluster_z,
    const int* flat_ci, const int* flat_cj, const int* num_flat_pairs_device,
    int pair_capacity, const int* dirty_flags, bool periodic, work_t cutoff_sq,
    const real_t* cell_matrix, const real_t* inv_cell, uint64_t* flat_mask,
    int prune_step, int prune_period, cudaStream_t stream) {
    if (pair_capacity == 0)
        return;

    int threads_needed = pair_capacity * kClusterSize;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (threads_needed + block - 1) / block;

    smart_rolling_prune<work_t><<<grid, block, 0, stream>>>(
        cluster_x, cluster_y, cluster_z, flat_ci, flat_cj,
        num_flat_pairs_device, dirty_flags, periodic, cutoff_sq, cell_matrix,
        inv_cell, flat_mask, prune_step, prune_period);
    ATV_CUDA_CHECK(cudaGetLastError());
}

template void launch_smart_rolling_prune<float>(
    const float*, const float*, const float*, const int*, const int*,
    const int*, int, const int*, bool, float, const float*, const float*,
    uint64_t*, int, int, cudaStream_t);
template void launch_smart_rolling_prune<double>(
    const double*, const double*, const double*, const int*, const int*,
    const int*, int, const int*, bool, double, const double*, const double*,
    uint64_t*, int, int, cudaStream_t);

} // namespace ATV::md::nlist
