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

#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>

using namespace ATV::cuda::md::nlist;

static constexpr int kMaxStencil = 512;

namespace ATV {
namespace cuda {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

// ============================================================================
// State
// ============================================================================

struct NeighborList::State {
    NListBuffers<real_t> buf;

    DeviceVector<real_t> d_cell{9};
    DeviceVector<real_t> d_inv_cell{9};
    DeviceVector<int> d_dims{3};
    DeviceVector<int> d_ncells{1};
    DeviceVector<int> stencil_dx{kMaxStencil};
    DeviceVector<int> stencil_dy{kMaxStencil};
    DeviceVector<int> stencil_dz{kMaxStencil};
    DeviceVector<int> d_n_stencil{1};
    DeviceVector<int> d_pair_count{1};

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
    int n_atoms = ctx.frame.n_atoms();

    s.num_atoms = static_cast<int>(n_atoms);
    int n = s.num_atoms;
    int nc = (n + kClusterSize - 1) / kClusterSize;
    s.num_clusters = nc;
    s.periodic = true;

    work_t r_search = cutoff_ + skin_;

    // GPU setup: compute cell/inv_cell/dims/stencil from device box
    launch_setup_cell_device(frame.box, s.d_cell.data(), s.d_inv_cell.data(),
                             s.d_dims.data(), s.stencil_dx.data(),
                             s.stencil_dy.data(), s.stencil_dz.data(),
                             s.d_n_stencil.data(), s.d_ncells.data(), r_search,
                             s.periodic, kMaxStencil, stream);

    // Read dims from device (build-time, allowed)
    cudaStreamSynchronize(stream);
    cudaMemcpy(&s.dims_x, s.d_dims.data(), sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&s.dims_y, s.d_dims.data() + 1, sizeof(int),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(&s.dims_z, s.d_dims.data() + 2, sizeof(int),
               cudaMemcpyDeviceToHost);
    int ncells = s.dims_x * s.dims_y * s.dims_z;

    // Pair capacity (upper bound)
    long long pair_cap_ll = static_cast<long long>(nc) * (nc + 1) / 2;
    s.pair_capacity = static_cast<int>(std::min(
        pair_cap_ll, static_cast<long long>(std::numeric_limits<int>::max())));

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

    // Allocate the owned PairList.
    //
    // Capacity sized from physical density × neighborhood volume rather
    // than a flat per-atom guess. For dense molecular systems (e.g. UA
    // polyethylene at ρ ≈ 0.0195 atoms/Å³, cutoff 10.5 Å + skin 2 Å) a
    // pair can have ~160 neighbors within the search shell; a single
    // (i,j) entry is emitted per pair from the cluster export, so the
    // tight upper bound is n_pairs ≈ n · N_neighbors / 2. We multiply
    // by a 2× safety factor for inhomogeneous regions (the typical hot
    // spot during stage-1 NVELimit relaxation of an "as-built" PE
    // configuration is local overlap that briefly inflates the neighbor
    // count). Falls back to `max(n·50, 64k)` when the density estimate
    // is degenerate (atomless or zero-volume frame).
    //
    // The previous flat `n · 50` cap silently TRUNCATED overflow pairs:
    // a 10000-atom PE system with ~160 neighbors per atom hit ~800k
    // pairs while capacity was 500k → ~37% of pairs missing per step.
    // The dropped pairs had no LJ contribution, so atoms drifted toward
    // each other under bond/angle forces alone until LJ "rediscovered"
    // them at near-zero r → catastrophic E_pot spikes (10¹⁴ eV
    // observed) and irrecoverable energy injection that no Langevin
    // damp could compensate.
    double box_vol = 0.0;
    double cut_search = static_cast<double>(cutoff_ + skin_);
    {
        // Read the orthorhombic box volume from the device-side box[0],4,8.
        prec_t box_diag[3] = {prec_t(0), prec_t(0), prec_t(0)};
        cudaMemcpy(&box_diag[0], &ctx.frame.box[0], sizeof(prec_t),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(&box_diag[1], &ctx.frame.box[4], sizeof(prec_t),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(&box_diag[2], &ctx.frame.box[8], sizeof(prec_t),
                   cudaMemcpyDeviceToHost);
        box_vol = static_cast<double>(box_diag[0]) * box_diag[1] * box_diag[2];
    }
    int atom_pair_cap = 0;
    if (box_vol > 0.0 && n > 0) {
        double density = double(n) / box_vol;
        double sphere_vol = (4.0 / 3.0) * 3.141592653589793 * cut_search
                            * cut_search * cut_search;
        double n_neighbors_est = density * sphere_vol;
        // Cluster-based export emits each pair TWICE (Ci→Cj and Cj→Ci);
        // empirical entry count is ~2.4× density·sphere_vol per atom
        // (cluster padding rounds neighbor lists up to the next
        // 8-atom group, so the effective sphere is slightly larger
        // than the geometric r_search). Use 3× as a comfortable cap.
        double per_atom = std::max(100.0, 3.0 * n_neighbors_est);
        atom_pair_cap = static_cast<int>(std::min(
            double(n) * per_atom,  // each atom appears in `per_atom` entries
            double(std::numeric_limits<int>::max())));
    }
    atom_pair_cap = std::max(atom_pair_cap, std::max(1, n * 50));
    cudaMalloc(&pairs_.i, atom_pair_cap * sizeof(int));
    cudaMalloc(&pairs_.j, atom_pair_cap * sizeof(int));
    cudaMalloc(&pairs_.dx, atom_pair_cap * sizeof(real_t));
    cudaMalloc(&pairs_.dy, atom_pair_cap * sizeof(real_t));
    cudaMalloc(&pairs_.dz, atom_pair_cap * sizeof(real_t));
    cudaMalloc(&pairs_.dist, atom_pair_cap * sizeof(real_t));
    cudaMalloc(&pairs_.n_pairs, sizeof(int));
    cudaMemset(pairs_.n_pairs, 0, sizeof(int));
    cudaMalloc(&pairs_.capacity, sizeof(int));
    int cap_val = static_cast<int>(atom_pair_cap);
    cudaMemcpy(pairs_.capacity, &cap_val, sizeof(int), cudaMemcpyHostToDevice);

    // Publish the pair list to the execution context for downstream
    // consumers (LJ, RDF), mirroring cpu::md::NeighborList.
    ctx.runtime.pair_list = &pairs_;
}

// ============================================================================
// launch() — per-step GPU work (pure GPU, zero D2H/sync/malloc)
// ============================================================================

void NeighborList::launch(const LaunchContext& ctx) {
    auto& s = *state_;
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;

    // Re-publish the pair list (address is stable across the run).
    ctx.runtime.pair_list = &pairs_;

    int n = s.num_atoms;
    int nc = s.num_clusters;
    work_t r_search = cutoff_ + skin_;

    // Recompute cell setup on GPU (box might change for NPT)
    launch_setup_cell_device(frame.box, s.d_cell.data(), s.d_inv_cell.data(),
                             s.d_dims.data(), s.stencil_dx.data(),
                             s.stencil_dy.data(), s.stencil_dz.data(),
                             s.d_n_stencil.data(), s.d_ncells.data(), r_search,
                             s.periodic, kMaxStencil, stream);

    // K1: Morton codes
    launch_compute_morton<work_t>(frame.atoms.x, frame.atoms.y, frame.atoms.z,
                                  n, s.d_inv_cell.data(), s.periodic,
                                  s.buf.morton_codes.data(),
                                  s.buf.atom_indices.data(), stream);

    // K2: Radix sort
    launch_radix_sort(s.buf.morton_codes.data(), s.buf.morton_sorted.data(),
                      s.buf.atom_indices.data(), s.buf.indices_sorted.data(), n,
                      s.buf.sort_temp, stream);

    // K3: Fused scatter + AABB
    launch_scatter_and_aabb<work_t>(
        s.buf.indices_sorted.data(), frame.atoms.x, frame.atoms.y,
        frame.atoms.z, n, nc, s.buf.cluster_x.data(), s.buf.cluster_y.data(),
        s.buf.cluster_z.data(), s.buf.original_index.data(),
        s.buf.cluster_count.data(), s.buf.aabb_min_x.data(),
        s.buf.aabb_min_y.data(), s.buf.aabb_min_z.data(),
        s.buf.aabb_max_x.data(), s.buf.aabb_max_y.data(),
        s.buf.aabb_max_z.data(), stream);

    // K4: Assign clusters to grid cells
    launch_assign_clusters_to_cells<work_t>(
        s.buf.aabb_min_x.data(), s.buf.aabb_min_y.data(),
        s.buf.aabb_min_z.data(), s.buf.aabb_max_x.data(),
        s.buf.aabb_max_y.data(), s.buf.aabb_max_z.data(), nc,
        s.d_inv_cell.data(), s.dims_x, s.dims_y, s.dims_z, s.periodic,
        s.buf.cluster_cell_idx.data(), stream);

    // K5: Build grid CSR
    int ncells = s.dims_x * s.dims_y * s.dims_z;
    launch_build_grid_csr(s.buf.cluster_cell_idx.data(), nc, ncells,
                          s.buf.cell_start.data(), s.buf.cell_clusters.data(),
                          s.buf.histogram.data(), s.buf.scan_temp, stream);

    work_t cut_skin_sq = r_search * r_search;

    // K6: Pair search pipeline
    launch_coarse_count_pairs(
        s.buf.cell_start.data(), s.buf.cluster_cell_idx.data(), s.dims_x,
        s.dims_y, s.dims_z, s.stencil_dx.data(), s.stencil_dy.data(),
        s.stencil_dz.data(), s.d_n_stencil.data(), nc, s.periodic,
        s.buf.pair_counts.data(), stream);

    launch_exclusive_sum(s.buf.pair_counts.data(), s.buf.row_ptr_upper.data(),
                         nc, s.buf.prefix_temp, stream);
    launch_write_row_ptr_sentinel(s.buf.row_ptr_upper.data(),
                                  s.buf.pair_counts.data(), nc, stream);

    launch_fill_cluster_pairs<work_t>(
        s.buf.cell_start.data(), s.buf.cell_clusters.data(),
        s.buf.cluster_cell_idx.data(), s.dims_x, s.dims_y, s.dims_z,
        s.stencil_dx.data(), s.stencil_dy.data(), s.stencil_dz.data(),
        s.d_n_stencil.data(), s.buf.aabb_min_x.data(), s.buf.aabb_min_y.data(),
        s.buf.aabb_min_z.data(), s.buf.aabb_max_x.data(),
        s.buf.aabb_max_y.data(), s.buf.aabb_max_z.data(), nc, s.periodic,
        cut_skin_sq, s.d_cell.data(), s.d_inv_cell.data(),
        s.buf.row_ptr_upper.data(), s.buf.col_idx.data(),
        s.buf.pair_counts.data(), stream);

    launch_exclusive_sum(s.buf.pair_counts.data(), s.buf.row_ptr.data(), nc,
                         s.buf.prefix_temp, stream);
    launch_write_row_ptr_sentinel(s.buf.row_ptr.data(),
                                  s.buf.pair_counts.data(), nc, stream);

    launch_compact_pairs(s.buf.col_idx.data(), s.buf.row_ptr_upper.data(),
                         s.buf.row_ptr.data(), s.buf.pair_counts.data(), nc,
                         stream);

    // Masks + flatten
    int pc = s.pair_capacity;
    int ms2 = pc > 0 ? pc : 1;
    cudaMemsetAsync(s.buf.interaction_mask.data(), 0xFF, ms2 * sizeof(uint64_t),
                    stream);

    cudaMemcpyAsync(s.buf.flat_pair_counter.data(), s.buf.row_ptr.data() + nc,
                    sizeof(int), cudaMemcpyDeviceToDevice, stream);
    launch_csr_to_flat(s.buf.row_ptr.data(), s.buf.col_idx.data(),
                       s.buf.interaction_mask.data(), nc, s.buf.flat_ci.data(),
                       s.buf.flat_cj.data(), s.buf.flat_mask.data(), stream);

    // Mask refresh: replaces the all-ones initial mask above with the actual
    // per-atom-pair interaction mask computed from current cluster positions.
    //
    // CRITICAL — without this call:
    //   * flat_mask stays at 0xFF for every cluster pair (memset above).
    //   * export_atom_pairs emits ALL 64 atom-pair combos per cluster pair.
    //   * Intra-cluster (ci == cj) emits both (a,b) and (b,a) directions PLUS
    //     diagonal self-pairs → forces on intra-cluster pairs are DOUBLED
    //     by the LJ kernel's Newton-3 accumulation; self-pairs at r=0 are
    //     filtered by the `r >= kMinPairDist` guard but waste pair slots.
    //   * Inter-cluster pairs at r > cutoff but within cutoff+skin AABB get
    //     exported and only filtered (silently) by the LJ kernel's `r <
    //     cutoff` test — wasted compute and pair slots, but no force.
    // smart_rolling_prune both filters by interaction cutoff and applies the
    // half-Newton exclusion for ci == cj (only j_lane > i_lane bits kept),
    // which gives the correct per-atom-pair entry count downstream.
    work_t cutoff_sq = cutoff_ * cutoff_;
    launch_smart_rolling_prune<work_t>(
        s.buf.cluster_x.data(), s.buf.cluster_y.data(), s.buf.cluster_z.data(),
        s.buf.flat_ci.data(), s.buf.flat_cj.data(),
        s.buf.flat_pair_counter.data(), s.pair_capacity,
        s.buf.dirty_flags.data(), s.periodic, cutoff_sq, s.d_cell.data(),
        s.d_inv_cell.data(), s.buf.flat_mask.data(), /*prune_step*/ 0,
        /*prune_period*/ 1, stream);

    // Export to atom-level PairList
    int cap = s.pair_capacity;
    if (cap <= 0) {
        cudaMemsetAsync(pairs_.n_pairs, 0, sizeof(int), stream);
        return;
    }

    cudaMemsetAsync(s.d_pair_count.data(), 0, sizeof(int), stream);
    launch_export_atom_pairs(
        s.buf.flat_ci.data(), s.buf.flat_cj.data(), s.buf.flat_mask.data(),
        s.buf.flat_pair_counter.data(), cap, s.buf.original_index.data(),
        frame.atoms.x, frame.atoms.y, frame.atoms.z, s.periodic,
        s.d_cell.data(), s.d_inv_cell.data(), pairs_.i, pairs_.j, pairs_.dx,
        pairs_.dy, pairs_.dz, pairs_.dist, s.d_pair_count.data(),
        pairs_.capacity, stream);

    cudaMemcpyAsync(pairs_.n_pairs, s.d_pair_count.data(), sizeof(int),
                    cudaMemcpyDeviceToDevice, stream);
}

// ============================================================================
// destroy() — release build-time resources
// ============================================================================

void NeighborList::destroy(const DestroyContext& ctx) {
    // Free the owned PairList (allocated in build)
    cudaFree(pairs_.i);
    cudaFree(pairs_.j);
    cudaFree(pairs_.dx);
    cudaFree(pairs_.dy);
    cudaFree(pairs_.dz);
    cudaFree(pairs_.dist);
    cudaFree(pairs_.n_pairs);
    cudaFree(pairs_.capacity);
    pairs_ = {};
    ctx.runtime.pair_list = nullptr;

    // Free internal state (DeviceVectors auto-free in State destructor)
    state_.reset();
}

// ============================================================================
// Constructor / move
// ============================================================================

NeighborList::NeighborList(double cutoff, double skin)
    : cutoff_(cutoff), skin_(skin < 0.0 ? 0.3 * cutoff : skin) {}

NeighborList::~NeighborList() = default;
NeighborList::NeighborList(NeighborList&&) noexcept = default;
NeighborList& NeighborList::operator=(NeighborList&&) noexcept = default;

} // namespace cuda
} // namespace ATV
