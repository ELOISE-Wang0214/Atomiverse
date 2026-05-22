#pragma once

/// @file nlist_data.cuh
/// @brief GPU data structures for the cluster-based neighbor list.

#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/types.cuh"
#include "flat_pair_list.cuh"

#include <cstddef>
#include <cstdint>

namespace ATV::cuda::md::nlist {

/// @brief Number of atoms per cluster (fixed for warp-friendly sizing).
inline constexpr int kClusterSize = 8;

/// @brief Per-cluster atom data in SoA layout.
/// @note Non-owning view; lifetime managed by NListBuffers.
template <typename work_t>
struct ClusterData {
    work_t* x;           ///< [num_clusters * kClusterSize]
    work_t* y;           ///< [num_clusters * kClusterSize]
    work_t* z;           ///< [num_clusters * kClusterSize]
    int* original_index; ///< [num_clusters * kClusterSize] map to input order
    int* cluster_count;  ///< [num_clusters] actual atom count (<=kClusterSize)
    work_t* aabb_min_x;  ///< [num_clusters]
    work_t* aabb_min_y;  ///< [num_clusters]
    work_t* aabb_min_z;  ///< [num_clusters]
    work_t* aabb_max_x;  ///< [num_clusters]
    work_t* aabb_max_y;  ///< [num_clusters]
    work_t* aabb_max_z;  ///< [num_clusters]
    int num_clusters;
    int num_atoms;
};

/// @brief Cluster pair list in CSR format.
/// @note Non-owning view; lifetime managed by NListBuffers.
/// @note Half-Newton: pairs (ci, cj) with cj >= ci are stored.
///       For diagonal pairs (ci == cj), only upper-triangle interactions
///       (j > i) are set in the mask to avoid double-counting.
///       Force kernels must accumulate both i->j and j->i from each pair.
template <typename work_t>
struct ClusterPairList {
    int* row_ptr; ///< [num_clusters + 1] CSR row pointer
    int* col_idx; ///< [num_pairs] j-cluster index per pair
    uint64_t*
        interaction_mask; ///< [num_pairs] 64-bit: bit[i*8+j]=1 for interaction
    int num_pairs;
    int num_clusters;
};

/// @brief Spatial grid of clusters for stencil-based pair search.
/// @note Non-owning view; lifetime managed by NListBuffers.
template <typename work_t>
struct ClusterGrid {
    int* cell_start;    ///< [ncells + 1] CSR: cluster start per cell
    int* cell_clusters; ///< [num_clusters] cluster indices sorted by cell
    int dims_x, dims_y, dims_z;  ///< Grid dimensions
    work_t cell_matrix[9];       ///< Cell matrix (row-major)
    work_t inv_cell[9];          ///< Inverse cell matrix
    const real_t* d_cell_matrix; ///< Device pointer to row-major cell matrix
    const real_t* d_inv_cell;    ///< Device pointer to inverse cell matrix
};

} // namespace ATV::cuda::md::nlist
