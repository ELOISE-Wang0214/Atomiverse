/// @file pair_search.cu
/// @brief Stencil-based cluster pair search — single-pass with coarse count.
///
/// Pipeline: coarse_count → prefix_sum → fill(AABB) → prefix_sum → compact.
/// The expensive AABB traversal happens only ONCE (in the fill kernel).

#include "nlist_impl.cuh"
#include "triclinic.cuh"

#include <cuda_runtime.h>
#include "atomiverse/cuda/types.cuh"

namespace ATV::md::nlist {

// ---------------------------------------------------------------------------
// Device: Minimum AABB distance² between two clusters
// ---------------------------------------------------------------------------

template <typename work_t>
__device__ work_t min_aabb_dist_sq(work_t ci_min_x, work_t ci_min_y, work_t ci_min_z,
                                 work_t ci_max_x, work_t ci_max_y, work_t ci_max_z,
                                 work_t cj_min_x, work_t cj_min_y, work_t cj_min_z,
                                 work_t cj_max_x, work_t cj_max_y, work_t cj_max_z,
                                 bool periodic,
                                 const real_t* __restrict__ cell_matrix,
                                 const real_t* __restrict__ inv_cell) {
    work_t cx_i = (ci_min_x + ci_max_x) * work_t(0.5);
    work_t cy_i = (ci_min_y + ci_max_y) * work_t(0.5);
    work_t cz_i = (ci_min_z + ci_max_z) * work_t(0.5);
    work_t cx_j = (cj_min_x + cj_max_x) * work_t(0.5);
    work_t cy_j = (cj_min_y + cj_max_y) * work_t(0.5);
    work_t cz_j = (cj_min_z + cj_max_z) * work_t(0.5);

    work_t dx = cx_j - cx_i, dy = cy_j - cy_i, dz = cz_j - cz_i;
    if (periodic) minimum_image(dx, dy, dz, cell_matrix, inv_cell);

    work_t hx_i = (ci_max_x - ci_min_x) * work_t(0.5);
    work_t hy_i = (ci_max_y - ci_min_y) * work_t(0.5);
    work_t hz_i = (ci_max_z - ci_min_z) * work_t(0.5);
    work_t hx_j = (cj_max_x - cj_min_x) * work_t(0.5);
    work_t hy_j = (cj_max_y - cj_min_y) * work_t(0.5);
    work_t hz_j = (cj_max_z - cj_min_z) * work_t(0.5);

    work_t gx = fmax(work_t(0), fabs(dx) - hx_i - hx_j);
    work_t gy = fmax(work_t(0), fabs(dy) - hy_i - hy_j);
    work_t gz = fmax(work_t(0), fabs(dz) - hz_i - hz_j);
    return gx * gx + gy * gy + gz * gz;
}

// ---------------------------------------------------------------------------
// Kernel: Coarse pair count (no AABB check — just count stencil-cell clusters)
//
// One thread per i-cluster.  Much cheaper than the full AABB-based count
// because it only reads cell_start[] values, not AABB coordinates.
// Produces an upper-bound pair_counts for pre-allocation.
// ---------------------------------------------------------------------------

__global__ void coarse_count_pairs(
    const int* __restrict__ cell_start,
    const int* __restrict__ cluster_cell_idx,
    int dims_x, int dims_y, int dims_z,
    const int* __restrict__ stencil_dx,
    const int* __restrict__ stencil_dy,
    const int* __restrict__ stencil_dz,
    const int* __restrict__ num_stencil_ptr, int num_clusters, bool periodic,
    int* __restrict__ pair_counts) {
    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    if (ci >= num_clusters) return;

    int num_stencil = *num_stencil_ptr;
    int home = cluster_cell_idx[ci];
    int ncells_yz = dims_y * dims_z;
    int hx = home / ncells_yz;
    int hy = (home / dims_z) % dims_y;
    int hz = home % dims_z;

    int count = 0;
    for (int s = 0; s < num_stencil; ++s) {
        int nx = hx + stencil_dx[s];
        int ny = hy + stencil_dy[s];
        int nz = hz + stencil_dz[s];
        if (periodic) {
            nx = ((nx % dims_x) + dims_x) % dims_x;
            ny = ((ny % dims_y) + dims_y) % dims_y;
            nz = ((nz % dims_z) + dims_z) % dims_z;
        } else {
            if (nx < 0 || nx >= dims_x || ny < 0 || ny >= dims_y ||
                nz < 0 || nz >= dims_z)
                continue;
        }
        int cell_idx = nx * ncells_yz + ny * dims_z + nz;
        count += cell_start[cell_idx + 1] - cell_start[cell_idx];
    }
    pair_counts[ci] = count;
}

void launch_coarse_count_pairs(const int* cell_start,
                                const int* cluster_cell_idx,
                                int dims_x, int dims_y, int dims_z,
                                const int* stencil_dx,
                                const int* stencil_dy,
                                const int* stencil_dz,
                                const int* num_stencil_ptr, int num_clusters,
                                bool periodic, int* pair_counts,
                                cudaStream_t stream) {
    if (num_clusters == 0) return;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (num_clusters + block - 1) / block;
    coarse_count_pairs<<<grid, block, 0, stream>>>(
        cell_start, cluster_cell_idx, dims_x, dims_y, dims_z,
        stencil_dx, stencil_dy, stencil_dz,
        num_stencil_ptr, num_clusters, periodic, pair_counts);
    ATV_CUDA_CHECK(cudaGetLastError());
}

// ---------------------------------------------------------------------------
// Kernel: Fill cluster pairs (single AABB pass)
//
// Writes pairs at row_ptr_upper[ci] offsets (pre-allocated from coarse count).
// Outputs actual pair_counts[ci] for subsequent exact prefix sum.
// Uses cluster_cell_idx to avoid recomputing home cell (P5).
// ---------------------------------------------------------------------------

template <typename work_t>
__global__ __launch_bounds__(ATV::cuda::kDefaultBlock)
void fill_pairs(
    const int* __restrict__ cell_start,
    const int* __restrict__ cell_clusters,
    const int* __restrict__ cluster_cell_idx,
    int dims_x, int dims_y, int dims_z,
    const int* __restrict__ stencil_dx,
    const int* __restrict__ stencil_dy,
    const int* __restrict__ stencil_dz,
    const int* __restrict__ num_stencil_ptr,
    const real_t* __restrict__ aabb_min_x,
    const real_t* __restrict__ aabb_min_y,
    const real_t* __restrict__ aabb_min_z,
    const real_t* __restrict__ aabb_max_x,
    const real_t* __restrict__ aabb_max_y,
    const real_t* __restrict__ aabb_max_z,
    int num_clusters, bool periodic,
    work_t cutoff_plus_skin_sq,
    const real_t* __restrict__ cell_matrix,
    const real_t* __restrict__ inv_cell,
    const int* __restrict__ row_ptr_upper,
    int* __restrict__ col_idx,
    int* __restrict__ pair_counts) {
    int warp_id = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int lane = threadIdx.x % 32;

    if (warp_id >= num_clusters) return;

    int num_stencil = *num_stencil_ptr;
    int ci = warp_id;

    // Home cell from pre-computed cluster_cell_idx (P5: no cart_to_frac)
    int home = cluster_cell_idx[ci];
    int ncells_yz = dims_y * dims_z;
    int hx = home / ncells_yz;
    int hy = (home / dims_z) % dims_y;
    int hz = home % dims_z;

    int write_offset = row_ptr_upper[ci];
    int count = 0;

    for (int s = 0; s < num_stencil; ++s) {
        int nx = hx + stencil_dx[s];
        int ny = hy + stencil_dy[s];
        int nz = hz + stencil_dz[s];

        if (periodic) {
            nx = ((nx % dims_x) + dims_x) % dims_x;
            ny = ((ny % dims_y) + dims_y) % dims_y;
            nz = ((nz % dims_z) + dims_z) % dims_z;
        } else {
            if (nx < 0 || nx >= dims_x || ny < 0 || ny >= dims_y ||
                nz < 0 || nz >= dims_z)
                continue;
        }

        int cell_idx = nx * ncells_yz + ny * dims_z + nz;
        int start = cell_start[cell_idx];
        int end = cell_start[cell_idx + 1];
        int num_cj = end - start;

        for (int batch = 0; batch < (num_cj + 31) / 32; ++batch) {
            int local = batch * 32 + lane;
            bool passes = false;
            int cj_val = -1;

            if (local < num_cj) {
                cj_val = cell_clusters[start + local];
                if (cj_val >= ci) {
                    work_t dist_sq = min_aabb_dist_sq(
                        aabb_min_x[ci], aabb_min_y[ci], aabb_min_z[ci],
                        aabb_max_x[ci], aabb_max_y[ci], aabb_max_z[ci],
                        aabb_min_x[cj_val], aabb_min_y[cj_val],
                        aabb_min_z[cj_val],
                        aabb_max_x[cj_val], aabb_max_y[cj_val],
                        aabb_max_z[cj_val],
                        periodic, cell_matrix, inv_cell);
                    passes = (dist_sq <= cutoff_plus_skin_sq);
                }
            }

            unsigned ballot = __ballot_sync(0xFFFFFFFF, passes);
            int prefix = __popc(ballot & ((1u << lane) - 1));

            if (passes) {
                col_idx[write_offset + prefix] = cj_val;
            }

            int batch_count = __popc(ballot);
            write_offset += batch_count;
            count += batch_count;
        }
    }

    if (lane == 0) {
        pair_counts[ci] = count;
    }
}

template <typename work_t>
void launch_fill_cluster_pairs(const int* cell_start,
                               const int* cell_clusters,
                               const int* cluster_cell_idx,
                               int dims_x, int dims_y, int dims_z,
                               const int* stencil_dx,
                               const int* stencil_dy,
                               const int* stencil_dz,
                               const int* num_stencil_ptr,
                               const real_t* aabb_min_x,
                               const real_t* aabb_min_y,
                               const real_t* aabb_min_z,
                               const real_t* aabb_max_x,
                               const real_t* aabb_max_y,
                               const real_t* aabb_max_z,
                               int num_clusters, bool periodic,
                               work_t cutoff_plus_skin_sq,
                               const real_t* cell_matrix,
                               const real_t* inv_cell,
                               const int* row_ptr_upper,
                               int* col_idx,
                               int* pair_counts,
                               cudaStream_t stream) {
    if (num_clusters == 0) return;
    int threads_needed = num_clusters * 32;
    int block = ATV::cuda::kDefaultBlock;
    int grid = (threads_needed + block - 1) / block;
    fill_pairs<work_t><<<grid, block, 0, stream>>>(
        cell_start, cell_clusters, cluster_cell_idx,
        dims_x, dims_y, dims_z,
        stencil_dx, stencil_dy, stencil_dz, num_stencil_ptr,
        aabb_min_x, aabb_min_y, aabb_min_z,
        aabb_max_x, aabb_max_y, aabb_max_z,
        num_clusters, periodic, cutoff_plus_skin_sq,
        cell_matrix, inv_cell, row_ptr_upper, col_idx, pair_counts);
    ATV_CUDA_CHECK(cudaGetLastError());
}

template void launch_fill_cluster_pairs<float>(
    const int*, const int*, const int*, int, int, int,
    const int*, const int*, const int*, const int*,
    const float*, const float*, const float*,
    const float*, const float*, const float*,
    int, bool, float, const float*, const float*,
    const int*, int*, int*, cudaStream_t);
template void launch_fill_cluster_pairs<double>(
    const int*, const int*, const int*, int, int, int,
    const int*, const int*, const int*, const int*,
    const double*, const double*, const double*,
    const double*, const double*, const double*,
    int, bool, double, const double*, const double*,
    const int*, int*, int*, cudaStream_t);

// ---------------------------------------------------------------------------
// Kernel: Write row_ptr sentinel on device (avoids D→H readback)
// ---------------------------------------------------------------------------

__global__ void write_row_ptr_sentinel(
    int* __restrict__ row_ptr, const int* __restrict__ pair_counts,
    int num_clusters) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        row_ptr[num_clusters] =
            row_ptr[num_clusters - 1] + pair_counts[num_clusters - 1];
    }
}

void launch_write_row_ptr_sentinel(int* row_ptr, const int* pair_counts,
                                    int num_clusters, cudaStream_t stream) {
    write_row_ptr_sentinel<<<1, 1, 0, stream>>>(
        row_ptr, pair_counts, num_clusters);
    ATV_CUDA_CHECK(cudaGetLastError());
}

// ---------------------------------------------------------------------------
// Kernel: In-place compaction of col_idx from upper-bound to exact CSR layout.
// Safe because row_ptr[ci] <= row_ptr_upper[ci] for all ci (exact count <=
// coarse count implies smaller prefix sums).
// ---------------------------------------------------------------------------

__global__ void compact_pairs(
    int* __restrict__ col_idx,
    const int* __restrict__ row_ptr_upper,
    const int* __restrict__ row_ptr,
    const int* __restrict__ pair_counts,
    int num_clusters) {
    // One warp per cluster for parallel copy
    int warp_id = (blockIdx.x * blockDim.x + threadIdx.x) / 32;
    int lane = threadIdx.x % 32;

    if (warp_id >= num_clusters) return;

    int ci = warp_id;
    int src = row_ptr_upper[ci];
    int dst = row_ptr[ci];
    int n = pair_counts[ci];

    // Warp-cooperative copy: each lane handles elements at stride 32
    for (int i = lane; i < n; i += 32) {
        col_idx[dst + i] = col_idx[src + i];
    }
}

void launch_compact_pairs(int* col_idx, const int* row_ptr_upper,
                           const int* row_ptr, const int* pair_counts,
                           int num_clusters, cudaStream_t stream) {
    if (num_clusters == 0) return;
    int block = ATV::cuda::kDefaultBlock;
    int threads_needed = num_clusters * 32;
    int grid = (threads_needed + block - 1) / block;
    compact_pairs<<<grid, block, 0, stream>>>(
        col_idx, row_ptr_upper, row_ptr, pair_counts, num_clusters);
    ATV_CUDA_CHECK(cudaGetLastError());
}

} // namespace ATV::md::nlist
