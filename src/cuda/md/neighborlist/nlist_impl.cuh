#pragma once

/// @file nlist_impl.cuh
/// @brief Internal implementation class bridging public API and CUDA kernels.

#include "nlist_data.cuh"
#include "nlist_memory.cuh"
#include "triclinic.cuh"

#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

#include <cstddef>
#include <cstdint>

namespace ATV::cuda::md::nlist {

/// @brief Manages all device memory and kernel orchestration.
template <typename work_t>
struct NListBuffers {
    // Morton sort temporaries
    DeviceVector<uint32_t> morton_codes;
    DeviceVector<int> atom_indices;
    DeviceVector<uint32_t> morton_sorted;
    DeviceVector<int> indices_sorted;
    DeviceVector<char> sort_temp;

    // Cluster SoA
    DeviceVector<work_t> cluster_x, cluster_y, cluster_z;
    DeviceVector<int> original_index;
    DeviceVector<int> cluster_count;

    // AABB
    DeviceVector<work_t> aabb_min_x, aabb_min_y, aabb_min_z;
    DeviceVector<work_t> aabb_max_x, aabb_max_y, aabb_max_z;

    // Grid
    DeviceVector<int> cluster_cell_idx; // [num_clusters] cell assignment
    DeviceVector<int> cell_start;       // [ncells + 1]
    DeviceVector<int> cell_clusters;    // [num_clusters]
    DeviceVector<int> histogram;        // temp for building CSR
    DeviceVector<char> scan_temp;

    // Pair list (CSR)
    DeviceVector<int> pair_counts;   // [num_clusters] per-cluster pair count
    DeviceVector<int> row_ptr;       // [num_clusters + 1]
    DeviceVector<int> row_ptr_upper; // [num_clusters + 1] coarse-count offsets
    DeviceVector<int> col_idx;       // [num_pairs]
    DeviceVector<uint64_t> interaction_mask; // [num_pairs]
    DeviceVector<char> prefix_temp;

    // Reference positions for the Verlet-skin rebuild check
    // (snapshot of frame.atoms.x/y/z taken at the end of every full
    // rebuild; consumed by launch_atom_displacement_max each step).
    DeviceVector<work_t> ref_x, ref_y, ref_z;

    // Stencil (device copy)
    DeviceVector<int> stencil_dx, stencil_dy, stencil_dz;

    // Flat pair list (fused pipeline)
    DeviceVector<int> flat_ci, flat_cj;
    DeviceVector<uint64_t> flat_mask;
    DeviceVector<int> flat_pair_counter; // single int, atomic counter

    // Dirty flags written by launch_cluster_displacement; currently
    // accepted-but-unread by smart_rolling_prune (reserved for a
    // future periodic-skip optimization, see cuda-md-review-cleanup-05*).
    DeviceVector<int> dirty_flags;

    void resize_atoms(int n, int n_clusters) {
        morton_codes.resize(n);
        atom_indices.resize(n);
        morton_sorted.resize(n);
        indices_sorted.resize(n);

        int padded = n_clusters * kClusterSize;
        cluster_x.resize(padded);
        cluster_y.resize(padded);
        cluster_z.resize(padded);
        original_index.resize(padded);
        cluster_count.resize(n_clusters);

        aabb_min_x.resize(n_clusters);
        aabb_min_y.resize(n_clusters);
        aabb_min_z.resize(n_clusters);
        aabb_max_x.resize(n_clusters);
        aabb_max_y.resize(n_clusters);
        aabb_max_z.resize(n_clusters);

        cluster_cell_idx.resize(n_clusters);

        pair_counts.resize(n_clusters);
        row_ptr.resize(n_clusters + 1);
        row_ptr_upper.resize(n_clusters + 1);

        ref_x.resize(n);
        ref_y.resize(n);
        ref_z.resize(n);
    }

    void resize_grid(int ncells) {
        cell_start.resize(ncells + 1);
        cell_clusters.resize(cluster_x.size() / kClusterSize);
        histogram.resize(ncells);
    }

    void resize_pairs(int num_pairs) {
        col_idx.resize(num_pairs);
        interaction_mask.resize(num_pairs);
    }

    void resize_flat_pairs(int max_pairs) {
        flat_ci.resize(max_pairs);
        flat_cj.resize(max_pairs);
        flat_mask.resize(max_pairs);
        flat_pair_counter.resize(1);
    }

    void resize_displacement(int n_clusters) { dirty_flags.resize(n_clusters); }

    ClusterData<work_t> cluster_data(int num_clusters, int num_atoms) {
        return {cluster_x.data(),
                cluster_y.data(),
                cluster_z.data(),
                original_index.data(),
                cluster_count.data(),
                aabb_min_x.data(),
                aabb_min_y.data(),
                aabb_min_z.data(),
                aabb_max_x.data(),
                aabb_max_y.data(),
                aabb_max_z.data(),
                num_clusters,
                num_atoms};
    }

    ClusterPairList<work_t> pair_list(int num_clusters, int num_pairs) {
        return {row_ptr.data(), col_idx.data(), interaction_mask.data(),
                num_pairs, num_clusters};
    }

    FlatPairList<work_t> flat_pair_list() {
        return {flat_ci.data(), flat_cj.data(), flat_mask.data(),
                flat_pair_counter.data(), static_cast<int>(flat_ci.size())};
    }
};

// Kernel launch declarations (defined in .cu files)
template <typename work_t>
void launch_compute_morton(const real_t* pos_x, const real_t* pos_y,
                           const real_t* pos_z, int num_atoms,
                           const real_t* inv_cell, bool periodic,
                           uint32_t* morton_codes, int* atom_indices,
                           cudaStream_t stream = 0);

template <typename work_t>
void launch_assign_clusters_to_cells(
    const real_t* aabb_min_x, const real_t* aabb_min_y,
    const real_t* aabb_min_z, const real_t* aabb_max_x,
    const real_t* aabb_max_y, const real_t* aabb_max_z, int num_clusters,
    const real_t* inv_cell, int dims_x, int dims_y, int dims_z, bool periodic,
    int* cluster_cell_idx, cudaStream_t stream = 0);

// Coarse count: upper-bound pair count per cluster (no AABB, O(stencil) per ci)
void launch_coarse_count_pairs(const int* cell_start,
                               const int* cluster_cell_idx, int dims_x,
                               int dims_y, int dims_z, const int* stencil_dx,
                               const int* stencil_dy, const int* stencil_dz,
                               const int* num_stencil_ptr, int num_clusters,
                               bool periodic, int* pair_counts,
                               cudaStream_t stream = 0);

// Single-pass AABB fill: writes pairs at row_ptr_upper offsets, outputs actual pair_counts
template <typename work_t>
void launch_fill_cluster_pairs(
    const int* cell_start, const int* cell_clusters,
    const int* cluster_cell_idx, int dims_x, int dims_y, int dims_z,
    const int* stencil_dx, const int* stencil_dy, const int* stencil_dz,
    const int* num_stencil_ptr, const real_t* aabb_min_x,
    const real_t* aabb_min_y, const real_t* aabb_min_z,
    const real_t* aabb_max_x, const real_t* aabb_max_y,
    const real_t* aabb_max_z, int num_clusters, bool periodic,
    work_t cutoff_plus_skin_sq, const real_t* cell_matrix,
    const real_t* inv_cell, const int* row_ptr_upper, int* col_idx,
    int* pair_counts, cudaStream_t stream = 0);

// Write row_ptr[N] = row_ptr[N-1] + pair_counts[N-1] on device (avoids D→H readback)
void launch_write_row_ptr_sentinel(int* row_ptr, const int* pair_counts,
                                   int num_clusters, cudaStream_t stream = 0);

// In-place compaction of col_idx from upper-bound to exact CSR layout
void launch_compact_pairs(int* col_idx, const int* row_ptr_upper,
                          const int* row_ptr, const int* pair_counts,
                          int num_clusters, cudaStream_t stream = 0);

// CUB sort wrapper (manages temp storage via DeviceVector)
void launch_radix_sort(uint32_t* keys_in, uint32_t* keys_out, int* values_in,
                       int* values_out, int num_items, DeviceVector<char>& temp,
                       cudaStream_t stream = 0);

// CUB prefix sum wrapper (manages temp storage via DeviceVector)
void launch_exclusive_sum(const int* input, int* output, int num_items,
                          DeviceVector<char>& temp, cudaStream_t stream = 0);

// Query the maximum temp-storage bytes that CUB radix sort and exclusive
// sum will require for `n_atoms` items / `n_clusters` items, so the
// caller can pre-size DeviceVector<char> temps in build() and avoid any
// cudaMalloc in the per-step launch() path.
//
// Returns {sort_bytes, scan_bytes_atoms, scan_bytes_clusters}. The
// extra `scan_bytes_clusters` covers `prefix_temp` used over
// `pair_counts[num_clusters]` in `launch_exclusive_sum` from
// `neighbor_list.cu::launch()`.
struct CubTempBytes {
    std::size_t sort_bytes;
    std::size_t scan_bytes_atoms;
    std::size_t scan_bytes_clusters;
};
CubTempBytes query_cub_temp_bytes(int num_atoms, int num_clusters);

// Grid CSR build
void launch_build_grid_csr(const int* cluster_cell_idx, int num_clusters,
                           int ncells, int* cell_start, int* cell_clusters,
                           int* histogram, DeviceVector<char>& scan_temp,
                           cudaStream_t stream = 0);

// --- Fused pipeline kernels (Phase 2-4) ---

// Fused AoS→SoA + Morton code computation (replaces separate K1+K2)
template <typename work_t>
void launch_aos_to_soa_morton(const double* positions_aos, int num_atoms,
                              real_t* x, real_t* y, real_t* z,
                              const real_t* inv_cell, bool periodic,
                              uint32_t* morton_codes, int* atom_indices,
                              cudaStream_t stream = 0);

// Fused scatter-to-cluster + AABB computation (replaces separate K4+K5)
template <typename work_t>
void launch_scatter_and_aabb(const int* sorted_indices, const real_t* pos_x,
                             const real_t* pos_y, const real_t* pos_z,
                             int num_atoms, int num_clusters, real_t* cluster_x,
                             real_t* cluster_y, real_t* cluster_z,
                             int* original_index, int* cluster_count,
                             real_t* aabb_min_x, real_t* aabb_min_y,
                             real_t* aabb_min_z, real_t* aabb_max_x,
                             real_t* aabb_max_y, real_t* aabb_max_z,
                             cudaStream_t stream = 0);

// CSR→Flat expansion (one warp per cluster, parallel copy)
void launch_csr_to_flat(const int* row_ptr, const int* col_idx,
                        const uint64_t* interaction_mask, int num_clusters,
                        int* flat_ci, int* flat_cj, uint64_t* flat_mask,
                        cudaStream_t stream = 0);

// Per-cluster displacement + dirty flags + rebuild flag
template <typename work_t>
void launch_cluster_displacement(
    const real_t* cluster_x, const real_t* cluster_y, const real_t* cluster_z,
    const real_t* ref_cluster_x, const real_t* ref_cluster_y,
    const real_t* ref_cluster_z, const int* cluster_count, int num_clusters,
    work_t dirty_threshold_sq, real_t* cluster_disp_sq, int* dirty_flags,
    int* rebuild_flag, work_t half_skin_sq, cudaStream_t stream = 0);

// Per-atom displacement vs ref_{x,y,z}; sets rebuild_flag=1 on first
// thread whose squared displacement exceeds half_skin_sq. Reads
// `frame.atoms.x/y/z` directly (always live), unlike
// `launch_cluster_displacement` which would compare stale sorted
// `cluster_x/y/z` to its ref snapshot between rebuilds. Used by
// `cuda::NeighborList::launch()` as the Verlet-skin gate.
template <typename work_t>
void launch_atom_displacement_max(const real_t* pos_x, const real_t* pos_y,
                                  const real_t* pos_z, const real_t* ref_x,
                                  const real_t* ref_y, const real_t* ref_z,
                                  int num_atoms, work_t half_skin_sq,
                                  int* rebuild_flag, cudaStream_t stream = 0);

// Snapshot current per-atom positions into ref_{x,y,z} buffers. Called
// after every full neighbor-list rebuild so the next step's
// `launch_atom_displacement_max` measures motion since the rebuild.
void launch_snapshot_ref_atoms(const real_t* pos_x, const real_t* pos_y,
                               const real_t* pos_z, real_t* ref_x,
                               real_t* ref_y, real_t* ref_z, int num_atoms,
                               cudaStream_t stream = 0);

// Smart rolling prune on flat pair list (skip clean pairs)
template <typename work_t>
void launch_smart_rolling_prune(
    const real_t* cluster_x, const real_t* cluster_y, const real_t* cluster_z,
    const int* flat_ci, const int* flat_cj, const int* num_flat_pairs_device,
    int pair_capacity, const int* dirty_flags, bool periodic, work_t cutoff_sq,
    const real_t* cell_matrix, const real_t* inv_cell, uint64_t* flat_mask,
    int prune_step, int prune_period, cudaStream_t stream = 0);

// --- GPU-side setup (setup.cu) ---

void launch_setup_cell_device(const prec_t* d_box, real_t* d_cell,
                              real_t* d_inv_cell, int* d_dims,
                              int* d_stencil_dx, int* d_stencil_dy,
                              int* d_stencil_dz, int* d_n_stencil,
                              int* d_ncells, work_t cutoff_plus_skin,
                              bool periodic, int max_stencil,
                              cudaStream_t stream = 0);

void launch_setup_cell(const real_t* box_data, real_t* d_cell,
                       real_t* d_inv_cell, int* d_dims, int* d_stencil_dx,
                       int* d_stencil_dy, int* d_stencil_dz, int* d_n_stencil,
                       int* d_ncells, work_t cutoff_plus_skin, bool periodic,
                       int max_stencil, cudaStream_t stream = 0);

void launch_write_int(int* arr, int idx, int val, cudaStream_t stream = 0);

// --- Export to atom pairs (export_pairs.cu) ---

void launch_export_atom_pairs(const int* flat_ci, const int* flat_cj,
                              const uint64_t* flat_mask,
                              const int* num_flat_device, int n_flat,
                              const int* original_index, const real_t* pos_x,
                              const real_t* pos_y, const real_t* pos_z,
                              bool periodic, const real_t* cell_matrix,
                              const real_t* inv_cell, int* out_i, int* out_j,
                              real_t* out_dx, real_t* out_dy, real_t* out_dz,
                              real_t* out_dist, int* out_count,
                              const int* capacity_ptr, cudaStream_t stream = 0);

} // namespace ATV::cuda::md::nlist
