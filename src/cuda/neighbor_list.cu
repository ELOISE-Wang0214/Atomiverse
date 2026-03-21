/// @file md/neighbor_list.cu
/// @brief NeighborList Component — build/launch/destroy lifecycle.
///
/// build():   allocates PairList + internal buffers (one-time, may D2H/sync)
/// launch():  rebuilds pair list every step (pure GPU, zero D2H/sync/malloc)
/// destroy(): frees all GPU resources

#include "atomiverse/cuda/neighbor_list.cuh"
#include "md/neighborlist/nlist_impl.cuh"

#include <algorithm>
#include <limits>

#include <cuda_runtime.h>

using Real = ATV::Real;
using Int = ATV::Int;
using namespace ATV::md::nlist;

static constexpr int kMaxStencil = 512;

namespace ATV {

// ============================================================================
// State
// ============================================================================

struct NeighborList::State {
    NListBuffers<Real> buf;

    DeviceVector<Real> d_cell{9};
    DeviceVector<Real> d_inv_cell{9};
    DeviceVector<int> d_dims{3};
    DeviceVector<int> d_ncells{1};
    DeviceVector<int> stencil_dx{kMaxStencil};
    DeviceVector<int> stencil_dy{kMaxStencil};
    DeviceVector<int> stencil_dz{kMaxStencil};
    DeviceVector<int> d_n_stencil{1};
    DeviceVector<Int> d_pair_count{1};

    int num_atoms = 0;
    int num_clusters = 0;
    int pair_capacity = 0;
    int dims_x = 0, dims_y = 0, dims_z = 0;
    bool periodic = true;
};

// ============================================================================
// build() — one-time resource preparation
// ============================================================================

void NeighborList::build(const BuildContext& ctx) {
    state_ = std::make_unique<State>();
    auto& s = *state_;
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    Int n_atoms = ctx.frame.n_atoms();

    s.num_atoms = static_cast<int>(n_atoms);
    int n = s.num_atoms;
    int nc = (n + kClusterSize - 1) / kClusterSize;
    s.num_clusters = nc;
    s.periodic = true;

    Real r_search = cutoff_ + skin_;

    // GPU setup: compute cell/inv_cell/dims/stencil from device box
    launch_setup_cell_device(
        frame.box, s.d_cell.data(), s.d_inv_cell.data(),
        s.d_dims.data(), s.stencil_dx.data(), s.stencil_dy.data(),
        s.stencil_dz.data(), s.d_n_stencil.data(), s.d_ncells.data(),
        r_search, s.periodic, kMaxStencil, stream);

    // Read dims from device (build-time, allowed)
    cudaStreamSynchronize(stream);
    cudaMemcpy(&s.dims_x, s.d_dims.data(), sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s.dims_y, s.d_dims.data() + 1, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s.dims_z, s.d_dims.data() + 2, sizeof(int), cudaMemcpyDeviceToHost);
    int ncells = s.dims_x * s.dims_y * s.dims_z;

    // Pair capacity (upper bound)
    long long pair_cap_ll = static_cast<long long>(nc) * (nc + 1) / 2;
    s.pair_capacity = static_cast<int>(
        std::min(pair_cap_ll, static_cast<long long>(std::numeric_limits<int>::max())));

    // Allocate ALL internal buffers
    s.buf.resize_atoms(n, nc);
    s.buf.resize_displacement(nc);
    s.buf.resize_grid(ncells);
    int pc = s.pair_capacity;
    s.buf.row_ptr_upper.resize(nc + 1);
    s.buf.row_ptr.resize(nc + 1);
    s.buf.col_idx.resize(pc > 0 ? pc : 1);
    int ms = pc > 0 ? pc : 1;
    s.buf.interaction_mask.resize(ms);
    s.buf.resize_flat_pairs(ms);

    // Allocate PairList in cuda::Frame
    int atom_pair_cap = std::max(1, n * 50);
    cudaMalloc(&frame.pairs.i, atom_pair_cap * sizeof(Int));
    cudaMalloc(&frame.pairs.j, atom_pair_cap * sizeof(Int));
    cudaMalloc(&frame.pairs.dx, atom_pair_cap * sizeof(Real));
    cudaMalloc(&frame.pairs.dy, atom_pair_cap * sizeof(Real));
    cudaMalloc(&frame.pairs.dz, atom_pair_cap * sizeof(Real));
    cudaMalloc(&frame.pairs.dist, atom_pair_cap * sizeof(Real));
    cudaMalloc(&frame.pairs.n_pairs, sizeof(Int));
    cudaMemset(frame.pairs.n_pairs, 0, sizeof(Int));
    cudaMalloc(&frame.pairs.capacity, sizeof(Int));
    Int cap_val = static_cast<Int>(atom_pair_cap);
    cudaMemcpy(frame.pairs.capacity, &cap_val, sizeof(Int), cudaMemcpyHostToDevice);
}

// ============================================================================
// launch() — per-step GPU work (pure GPU, zero D2H/sync/malloc)
// ============================================================================

void NeighborList::launch(const LaunchContext& ctx) {
    auto& s = *state_;
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;

    int n = s.num_atoms;
    int nc = s.num_clusters;
    Real r_search = cutoff_ + skin_;

    // Recompute cell setup on GPU (box might change for NPT)
    launch_setup_cell_device(
        frame.box, s.d_cell.data(), s.d_inv_cell.data(),
        s.d_dims.data(), s.stencil_dx.data(), s.stencil_dy.data(),
        s.stencil_dz.data(), s.d_n_stencil.data(), s.d_ncells.data(),
        r_search, s.periodic, kMaxStencil, stream);

    // K1: Morton codes
    launch_compute_morton<Real>(
        frame.x, frame.y, frame.z, n,
        s.d_inv_cell.data(), s.periodic,
        s.buf.morton_codes.data(), s.buf.atom_indices.data(), stream);

    // K2: Radix sort
    launch_radix_sort(
        s.buf.morton_codes.data(), s.buf.morton_sorted.data(),
        s.buf.atom_indices.data(), s.buf.indices_sorted.data(),
        n, s.buf.sort_temp, stream);

    // K3: Fused scatter + AABB
    launch_scatter_and_aabb<Real>(
        s.buf.indices_sorted.data(),
        frame.x, frame.y, frame.z, n, nc,
        s.buf.cluster_x.data(), s.buf.cluster_y.data(), s.buf.cluster_z.data(),
        s.buf.original_index.data(), s.buf.cluster_count.data(),
        s.buf.aabb_min_x.data(), s.buf.aabb_min_y.data(), s.buf.aabb_min_z.data(),
        s.buf.aabb_max_x.data(), s.buf.aabb_max_y.data(), s.buf.aabb_max_z.data(),
        stream);

    // K4: Assign clusters to grid cells
    launch_assign_clusters_to_cells<Real>(
        s.buf.aabb_min_x.data(), s.buf.aabb_min_y.data(), s.buf.aabb_min_z.data(),
        s.buf.aabb_max_x.data(), s.buf.aabb_max_y.data(), s.buf.aabb_max_z.data(),
        nc, s.d_inv_cell.data(), s.dims_x, s.dims_y, s.dims_z,
        s.periodic, s.buf.cluster_cell_idx.data(), stream);

    // K5: Build grid CSR
    int ncells = s.dims_x * s.dims_y * s.dims_z;
    launch_build_grid_csr(
        s.buf.cluster_cell_idx.data(), nc, ncells,
        s.buf.cell_start.data(), s.buf.cell_clusters.data(),
        s.buf.histogram.data(), s.buf.scan_temp, stream);

    Real cut_skin_sq = r_search * r_search;

    // K6: Pair search pipeline
    launch_coarse_count_pairs(
        s.buf.cell_start.data(), s.buf.cluster_cell_idx.data(),
        s.dims_x, s.dims_y, s.dims_z,
        s.stencil_dx.data(), s.stencil_dy.data(), s.stencil_dz.data(),
        s.d_n_stencil.data(), nc, s.periodic,
        s.buf.pair_counts.data(), stream);

    launch_exclusive_sum(
        s.buf.pair_counts.data(), s.buf.row_ptr_upper.data(),
        nc, s.buf.prefix_temp, stream);
    launch_write_row_ptr_sentinel(
        s.buf.row_ptr_upper.data(), s.buf.pair_counts.data(), nc, stream);

    launch_fill_cluster_pairs<Real>(
        s.buf.cell_start.data(), s.buf.cell_clusters.data(),
        s.buf.cluster_cell_idx.data(),
        s.dims_x, s.dims_y, s.dims_z,
        s.stencil_dx.data(), s.stencil_dy.data(), s.stencil_dz.data(),
        s.d_n_stencil.data(),
        s.buf.aabb_min_x.data(), s.buf.aabb_min_y.data(), s.buf.aabb_min_z.data(),
        s.buf.aabb_max_x.data(), s.buf.aabb_max_y.data(), s.buf.aabb_max_z.data(),
        nc, s.periodic, cut_skin_sq,
        s.d_cell.data(), s.d_inv_cell.data(),
        s.buf.row_ptr_upper.data(), s.buf.col_idx.data(),
        s.buf.pair_counts.data(), stream);

    launch_exclusive_sum(
        s.buf.pair_counts.data(), s.buf.row_ptr.data(),
        nc, s.buf.prefix_temp, stream);
    launch_write_row_ptr_sentinel(
        s.buf.row_ptr.data(), s.buf.pair_counts.data(), nc, stream);

    launch_compact_pairs(
        s.buf.col_idx.data(), s.buf.row_ptr_upper.data(),
        s.buf.row_ptr.data(), s.buf.pair_counts.data(), nc, stream);

    // Masks + flatten
    int pc = s.pair_capacity;
    int ms2 = pc > 0 ? pc : 1;
    cudaMemsetAsync(s.buf.interaction_mask.data(), 0xFF,
                    ms2 * sizeof(uint64_t), stream);

    cudaMemcpyAsync(s.buf.flat_pair_counter.data(),
                    s.buf.row_ptr.data() + nc,
                    sizeof(int), cudaMemcpyDeviceToDevice, stream);
    launch_csr_to_flat(
        s.buf.row_ptr.data(), s.buf.col_idx.data(),
        s.buf.interaction_mask.data(), nc,
        s.buf.flat_ci.data(), s.buf.flat_cj.data(),
        s.buf.flat_mask.data(), stream);

    // Export to atom-level PairList
    int cap = s.pair_capacity;
    if (cap <= 0) {
        cudaMemsetAsync(frame.pairs.n_pairs, 0, sizeof(Int), stream);
        return;
    }

    cudaMemsetAsync(s.d_pair_count.data(), 0, sizeof(Int), stream);
    launch_export_atom_pairs(
        s.buf.flat_ci.data(), s.buf.flat_cj.data(), s.buf.flat_mask.data(),
        s.buf.flat_pair_counter.data(), cap, s.buf.original_index.data(),
        frame.x, frame.y, frame.z,
        s.periodic, s.d_cell.data(), s.d_inv_cell.data(),
        frame.pairs.i, frame.pairs.j,
        frame.pairs.dx, frame.pairs.dy, frame.pairs.dz, frame.pairs.dist,
        s.d_pair_count.data(), frame.pairs.capacity, stream);

    cudaMemcpyAsync(frame.pairs.n_pairs, s.d_pair_count.data(),
                    sizeof(Int), cudaMemcpyDeviceToDevice, stream);
}

// ============================================================================
// destroy() — release build-time resources
// ============================================================================

void NeighborList::destroy(const DestroyContext& ctx) {
    auto& frame = ctx.frame;

    // Free PairList (allocated in build)
    cudaFree(frame.pairs.i); cudaFree(frame.pairs.j);
    cudaFree(frame.pairs.dx); cudaFree(frame.pairs.dy);
    cudaFree(frame.pairs.dz); cudaFree(frame.pairs.dist);
    cudaFree(frame.pairs.n_pairs); cudaFree(frame.pairs.capacity);
    frame.pairs = {};

    // Free internal state (DeviceVectors auto-free in State destructor)
    state_.reset();
}

// ============================================================================
// Constructor / move
// ============================================================================

NeighborList::NeighborList(Real cutoff, Real skin)
    : cutoff_(cutoff), skin_(skin < Real(0) ? Real(0.3) * cutoff : skin) {}

NeighborList::~NeighborList() = default;
NeighborList::NeighborList(NeighborList&&) noexcept = default;
NeighborList& NeighborList::operator=(NeighborList&&) noexcept = default;

} // namespace ATV
