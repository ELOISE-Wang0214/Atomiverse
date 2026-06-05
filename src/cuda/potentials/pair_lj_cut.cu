/// @file cuda/potentials/pair_lj_cut.cu
/// @brief CUDA PairLJCut — Lennard-Jones 12-6 shifted (kernel + class).
///        Includes LAMMPS-style `special_bonds` 1-2/1-3/1-4 exclusions.

#include "atomiverse/cuda/box.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/lj_kernel.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "atomiverse/cuda/types.cuh"
#include <cuda_runtime.h>
#include <set>
#include <stdexcept>
#include <vector>

namespace ATV {

using cuda::prec_t;
using cuda::real_t;
using cuda::work_t;

namespace cuda {
namespace kernels {

/// @brief Device-side binary search of a sorted int range.
///        Returns true if `target` is in `[beg, end)`.
__device__ inline bool excl_contains(const int* __restrict__ partners, int beg,
                                     int end, int target) {
    int lo = beg, hi = end;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        int v = partners[mid];
        if (v == target)
            return true;
        if (v < target)
            lo = mid + 1;
        else
            hi = mid;
    }
    return false;
}

/// @brief Atom-centric, **atomic-free** LJ force. One thread owns one atom
///        (cluster slot). It walks its cluster's FULL neighbor-cluster list
///        (cl_col_idx over [row_ptr[ci], row_ptr[ci+1])), computes each
///        pair distance from CURRENT atom positions (min-image), accumulates
///        the owned atom's force in registers, and writes it ONCE — no
///        atomicAdd, no flat pair list, no per-step export.
///
/// Replaces the old one-thread-per-pair kernel that did 6 atomicAdd/pair
/// (~160-way contention per force slot in a dense melt → the 76%-of-runtime
/// bottleneck). With a FULL neighbor list each interaction is computed twice
/// (once by each endpoint's owning thread), so the force on each atom is
/// complete without writing its neighbors; energy and virial are halved to
/// undo the double count. kClusterSize is the neighbor list's fixed cluster
/// width; pad slots carry a negative orig_index sentinel.
// Gather current atom positions into the cluster-sorted slot layout so the
// force kernel reads neighbour positions COALESCED (spx[cj*8+sj] is contiguous
// per cluster) instead of scattering x[orig_index[..]] ~80×/atom. One scattered
// read per slot here amortises the scatter that previously happened per pair —
// the dominant cost of lj_force_cluster (latency-bound uncoalesced gathers).
__global__ void gather_cluster_pos(const int* __restrict__ orig_index,
                                   const real_t* __restrict__ x,
                                   const real_t* __restrict__ y,
                                   const real_t* __restrict__ z, int n_slots,
                                   real_t* __restrict__ spx,
                                   real_t* __restrict__ spy,
                                   real_t* __restrict__ spz) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n_slots)
        return;
    int a = orig_index[k];
    if (a >= 0) {
        spx[k] = x[a];
        spy[k] = y[a];
        spz[k] = z[a];
    }
}

template <bool Ortho>
__global__ __launch_bounds__(ATV::cuda::kDefaultBlock, 4) void
lj_force_cluster(real_t* __restrict__ fx, real_t* __restrict__ fy,
                 real_t* __restrict__ fz, prec_t* __restrict__ energy,
                 prec_t* __restrict__ stress, const real_t* __restrict__ spx,
                 const real_t* __restrict__ spy, const real_t* __restrict__ spz,
                 const int* __restrict__ row_ptr,
                 const int* __restrict__ col_idx,
                 const int* __restrict__ orig_index, int num_clusters,
                 cuda::Box box, bool periodic, work_t epsilon, work_t sigma_sq,
                 work_t cutoff_sq, work_t min_pair_dist_sq, work_t v_shift,
                 const int* __restrict__ excl_offsets,
                 const int* __restrict__ excl_partners, work_t s_excl) {
    constexpr int CS = 8; // kClusterSize (neighbor-list cluster width)
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int ci = tid / CS;
    int slot = tid % CS;

    work_t my_energy = work_t(0);
    work_t my_vxx = work_t(0), my_vyy = work_t(0), my_vzz = work_t(0);
    work_t my_vxy = work_t(0), my_vxz = work_t(0), my_vyz = work_t(0);

    // Mixed precision (matches LAMMPS KOKKOS SINGLE_DOUBLE): geometry + force
    // in real_t (float — FP32 is 2x FP64 throughput on GH200 and the divide is
    // far cheaper), energy/virial accumulated in work_t (double). Force is
    // stored as real_t anyway, so float accumulation loses nothing.
    const real_t eps_r = real_t(epsilon);
    const real_t sigsq_r = real_t(sigma_sq);
    const real_t cut_r = real_t(cutoff_sq);
    const real_t minr_r = real_t(min_pair_dist_sq);
    const real_t vshift_r = real_t(v_shift);
    const real_t sexcl_r = real_t(s_excl);

    if (ci < num_clusters) {
        int ai = orig_index[ci * CS + slot];
        if (ai >= 0) {
            real_t xi = spx[ci * CS + slot], yi = spy[ci * CS + slot],
                   zi = spz[ci * CS + slot];
            real_t fxi = real_t(0), fyi = real_t(0), fzi = real_t(0);
            const bool has_excl = (excl_offsets != nullptr);
            int beg_excl = has_excl ? excl_offsets[ai] : 0;
            int end_excl = has_excl ? excl_offsets[ai + 1] : 0;

            int rstart = row_ptr[ci], rend = row_ptr[ci + 1];
            for (int p = rstart; p < rend; ++p) {
                int cj = col_idx[p];
#pragma unroll
                for (int sj = 0; sj < CS; ++sj) {
                    int aj = orig_index[cj * CS + sj];
                    if (aj < 0 || aj == ai)
                        continue;
                    real_t dx = spx[cj * CS + sj] - xi;
                    real_t dy = spy[cj * CS + sj] - yi;
                    real_t dz = spz[cj * CS + sj] - zi;
                    if (periodic) {
                        if constexpr (Ortho) {
                            dx -= box.cell[0] * nearbyintf(dx * box.inv_l[0]);
                            dy -= box.cell[4] * nearbyintf(dy * box.inv_l[1]);
                            dz -= box.cell[8] * nearbyintf(dz * box.inv_l[2]);
                        } else {
                            box.min_image(dx, dy, dz);
                        }
                    }
                    real_t r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 >= cut_r || r2 < minr_r)
                        continue;
                    real_t scale = real_t(1);
                    if (has_excl && end_excl > beg_excl
                        && excl_contains(excl_partners, beg_excl, end_excl,
                                         aj)) {
                        if (sexcl_r == real_t(0))
                            continue;
                        scale = sexcl_r;
                    }
                    real_t inv_r2 = real_t(1) / r2;
                    real_t sr2 = sigsq_r * inv_r2;
                    real_t sr6 = sr2 * sr2 * sr2;
                    real_t sr12 = sr6 * sr6;
                    real_t f_over_r = scale * real_t(24) * eps_r
                                      * (real_t(2) * sr12 - sr6) * inv_r2;
                    // dx = x[aj] - x[ai], so the force ON ai is -f_over_r·dx.
                    fxi += -f_over_r * dx;
                    fyi += -f_over_r * dy;
                    fzi += -f_over_r * dz;
                    // Each pair is visited twice (ai's and aj's threads): halve
                    // the energy/virial to recover the unique-pair sum. Promote
                    // the per-pair term to work_t for the double accumulation.
                    real_t e_pair =
                        scale * (real_t(4) * eps_r * (sr12 - sr6) - vshift_r);
                    my_energy += work_t(0.5) * work_t(e_pair);
                    my_vxx += work_t(0.5) * work_t(dx * (f_over_r * dx));
                    my_vyy += work_t(0.5) * work_t(dy * (f_over_r * dy));
                    my_vzz += work_t(0.5) * work_t(dz * (f_over_r * dz));
                    my_vxy += work_t(0.5) * work_t(dx * (f_over_r * dy));
                    my_vxz += work_t(0.5) * work_t(dx * (f_over_r * dz));
                    my_vyz += work_t(0.5) * work_t(dy * (f_over_r * dz));
                }
            }
            // Single non-atomic accumulate: this thread uniquely owns atom ai
            // (orig_index is a permutation), and force-phase kernels serialize
            // on the stream, so += safely adds to bonded contributions.
            fx[ai] += fxi;
            fy[ai] += fyi;
            fz[ai] += fzi;
        }
    }

    // Packed warp-shuffle reduction: 7 scalars (energy + 6 stress
    // components) reduced via cg::reduce. Replaces the original
    // sequential 7-pass shared-mem reductions with 7 × `__syncthreads`
    // (one per scalar inside `block_reduce_sum`) instead of 7 × 8
    // shared-mem barrier-steps. Drops shared-mem footprint from
    // `blockDim.x × sizeof(work_t)` to (blockDim.x / 32) × sizeof(work_t)
    // (cuda-md-review-cleanup-05-reductions-and-atomics).
    extern __shared__ work_t warp_scratch[]; // per-warp staging cells

    work_t my_vals[7] = {my_energy, my_vxx, my_vyy, my_vzz,
                         my_vxy,    my_vxz, my_vyz};
    work_t out[7];
    block_reduce_sum_n<7>(my_vals, warp_scratch, out);
    if (threadIdx.x == 0) {
        atomicAdd(energy, static_cast<prec_t>(out[0]));
        if (stress != nullptr) {
#pragma unroll
            for (int c = 0; c < 6; ++c)
                atomicAdd(&stress[c], static_cast<prec_t>(out[c + 1]));
        }
    }
}

} // namespace kernels
} // namespace cuda

namespace cuda {

PairLJCut::PairLJCut(double epsilon, double sigma, double cutoff)
    : Pair(cutoff), epsilon_(epsilon), sigma_(sigma) {
    v_shift_ = kernels::lj_shift(epsilon, sigma, cutoff);
}

PairLJCut& PairLJCut::set_special_bonds(double s12, double s13,
                                        double s14) & {
    // Supported combinations: all equal (incl. (0,0,0)) or all 1.0.
    if (s12 == 1.0 && s13 == 1.0 && s14 == 1.0) {
        s12_ = s13_ = s14_ = 1.0;
        s_excl_ = 1.0;
        has_exclusions_ = false;
        return *this;
    }
    if (s12 != s13 || s12 != s14) {
        throw std::invalid_argument(
            "cuda::PairLJCut::set_special_bonds: mixed scaling combinations "
            "(s12 != s13 or s12 != s14) are not supported in this revision. "
            "Pass three equal values (typically (0,0,0) for Dreiding-UA).");
    }
    s12_ = s13_ = s14_ = s12;
    s_excl_ = s12;
    has_exclusions_ = true; // exclusions enumerated at build()
    return *this;
}

PairLJCut PairLJCut::set_special_bonds(double s12, double s13,
                                       double s14) && {
    // Rvalue overload: delegate to the lvalue body, then move-return so
    // the chained `cuda::PairLJCut{…}.set_special_bonds(…)` expression
    // yields a movable rvalue (not a dangling lvalue reference to a
    // temporary). cf. `cpu::PairLJCut::set_special_bonds(…) &&`.
    this->set_special_bonds(s12, s13, s14);
    return std::move(*this);
}

void PairLJCut::build(const BuildContext& ctx) {
    if (!has_exclusions_) {
        n_excl_total_ = 0;
        return;
    }
    auto& f = ctx.frame;
    int n = f.n_atoms();
    int n_b = f.n_bonds();
    int n_a = f.n_angles();
    int n_d = f.n_dihedrals();

    // D2H read of topology (build-time cost; not in hot path).
    std::vector<int> bi(n_b), bj(n_b);
    std::vector<int> ai(n_a), ak(n_a);
    std::vector<int> di(n_d), dl(n_d);
    if (n_b > 0) {
        cudaMemcpy(bi.data(), f.bonds.atomi, n_b * sizeof(int),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(bj.data(), f.bonds.atomj, n_b * sizeof(int),
                   cudaMemcpyDeviceToHost);
    }
    if (n_a > 0) {
        cudaMemcpy(ai.data(), f.angles.atomi, n_a * sizeof(int),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(ak.data(), f.angles.atomk, n_a * sizeof(int),
                   cudaMemcpyDeviceToHost);
    }
    if (n_d > 0) {
        cudaMemcpy(di.data(), f.dihedrals.atomi, n_d * sizeof(int),
                   cudaMemcpyDeviceToHost);
        cudaMemcpy(dl.data(), f.dihedrals.atoml, n_d * sizeof(int),
                   cudaMemcpyDeviceToHost);
    }

    // Per-atom sorted-unique partner sets.
    std::vector<std::set<int>> excl(n);
    auto add_pair = [&](int a, int b) {
        if (a >= 0 && a < n && b >= 0 && b < n && a != b) {
            excl[a].insert(b);
            excl[b].insert(a);
        }
    };
    for (int k = 0; k < n_b; ++k) add_pair(bi[k], bj[k]);
    for (int k = 0; k < n_a; ++k) add_pair(ai[k], ak[k]);
    for (int k = 0; k < n_d; ++k) add_pair(di[k], dl[k]);

    // Pack into CSR.
    std::vector<int> offsets(n + 1, 0);
    for (int i = 0; i < n; ++i)
        offsets[i + 1] = offsets[i] + static_cast<int>(excl[i].size());
    n_excl_total_ = offsets[n];
    std::vector<int> partners(n_excl_total_);
    for (int i = 0; i < n; ++i) {
        int cursor = offsets[i];
        for (int p : excl[i])
            partners[cursor++] = p;
        // std::set already sorted ascending; nothing to do.
    }

    // Upload to device.
    cudaMalloc(&d_excl_offsets_, (n + 1) * sizeof(int));
    cudaMemcpy(d_excl_offsets_, offsets.data(), (n + 1) * sizeof(int),
               cudaMemcpyHostToDevice);
    if (n_excl_total_ > 0) {
        cudaMalloc(&d_excl_partners_, n_excl_total_ * sizeof(int));
        cudaMemcpy(d_excl_partners_, partners.data(),
                   n_excl_total_ * sizeof(int), cudaMemcpyHostToDevice);
    }
}

void PairLJCut::destroy(const DestroyContext& /*ctx*/) {
    if (d_excl_offsets_) {
        cudaFree(d_excl_offsets_);
        d_excl_offsets_ = nullptr;
    }
    if (d_excl_partners_) {
        cudaFree(d_excl_partners_);
        d_excl_partners_ = nullptr;
    }
    n_excl_total_ = 0;
    if (spx_) {
        cudaFree(spx_);
        cudaFree(spy_);
        cudaFree(spz_);
        spx_ = spy_ = spz_ = nullptr;
        spos_cap_ = 0;
    }
}

void PairLJCut::compute(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    const PairList* pairs = ctx.runtime.pair_list;
    if (!pairs)
        return;
    // Grid size = the NeighborList's pair_capacity upper bound (read once
    // from device). The kernel early-exits via `if (idx >= *n_pairs)` at
    // line 63, so over-launching by a few × is cheap (one global load +
    // early return per wasted thread) and far cheaper than the per-step
    // D2H + cudaStreamSynchronize that the previous grid-by-n_pairs
    // sizing required. The pair_capacity was previously read at build
    // time and is stable; we read it once here at first call and cache
    // the host int statically (per-instance via PIMPL would be cleaner
    // but adds State plumbing — the read is build-time-equivalent and
    // pays only on the first launch).
    //
    // The earlier bug — n·50 host estimate silently truncating overflow
    // pairs at ~37% loss on a 10k-atom PE system — was solved by the
    // density-based capacity rewrite in NeighborList::build()
    // (cuda-md-review-cleanup-02). With the new capacity sized to
    // ~3× n_neighbors_est, this kernel over-launches by ~3× actual
    // count; the early-exit absorbs it.
    // Atom-centric, atomic-free path: one thread per atom (cluster slot)
    // traverses the FULL cluster neighbor CSR published by NeighborList and
    // accumulates its own force from current positions. No flat pair list,
    // no per-step export, no atomicAdd on forces.
    if (pairs->cl_row_ptr == nullptr || pairs->num_clusters <= 0)
        return;
    constexpr int CS = 8; // kClusterSize
    int block = kDefaultBlock;
    int nslots = pairs->num_clusters * CS;
    int grid = grid_size(nslots);
    int shmem = warps_per_block(block) * sizeof(work_t);

    // Coalesced-position path: gather current positions into the cluster-sorted
    // slot layout once, then the force kernel reads neighbours contiguously.
    // (Re)allocate the sorted-position buffers when the slot count grows
    // (num_clusters changes across multi-GPU reneighbor intervals).
    if (nslots > spos_cap_) {
        if (spx_) {
            cudaFree(spx_);
            cudaFree(spy_);
            cudaFree(spz_);
        }
        cudaMalloc(&spx_, nslots * sizeof(real_t));
        cudaMalloc(&spy_, nslots * sizeof(real_t));
        cudaMalloc(&spz_, nslots * sizeof(real_t));
        spos_cap_ = nslots;
    }
    kernels::gather_cluster_pos<<<grid, block, 0, stream>>>(
        pairs->cl_orig_index, frame.atoms.x, frame.atoms.y, frame.atoms.z,
        nslots, spx_, spy_, spz_);
    // Pass nullptr for exclusion CSR when none configured — kernel takes
    // fast path and avoids the partner lookup overhead entirely.
    const int* excl_off = (has_exclusions_) ? d_excl_offsets_ : nullptr;
    const int* excl_par = (has_exclusions_) ? d_excl_partners_ : nullptr;
    const work_t cutoff_sq = static_cast<work_t>(cutoff_)
                             * static_cast<work_t>(cutoff_);
    const work_t sigma_sq = static_cast<work_t>(sigma_)
                            * static_cast<work_t>(sigma_);
    const work_t min_pair_dist_sq =
        static_cast<work_t>(kernels::kMinPairDist)
        * static_cast<work_t>(kernels::kMinPairDist);
    cuda::Box box = ctx.runtime.host_box;
    if (box.ortho) {
        kernels::lj_force_cluster<true><<<grid, block, shmem, stream>>>(
            frame.atoms.fx, frame.atoms.fy, frame.atoms.fz,
            &frame.scalars->energy, &frame.scalars->stress[0], spx_, spy_, spz_,
            pairs->cl_row_ptr, pairs->cl_col_idx, pairs->cl_orig_index,
            pairs->num_clusters, box, /*periodic=*/true, epsilon_, sigma_sq,
            cutoff_sq, min_pair_dist_sq, v_shift_, excl_off, excl_par,
            static_cast<work_t>(s_excl_));
    } else {
        kernels::lj_force_cluster<false><<<grid, block, shmem, stream>>>(
            frame.atoms.fx, frame.atoms.fy, frame.atoms.fz,
            &frame.scalars->energy, &frame.scalars->stress[0], spx_, spy_, spz_,
            pairs->cl_row_ptr, pairs->cl_col_idx, pairs->cl_orig_index,
            pairs->num_clusters, box, /*periodic=*/true, epsilon_, sigma_sq,
            cutoff_sq, min_pair_dist_sq, v_shift_, excl_off, excl_par,
            static_cast<work_t>(s_excl_));
    }
}

} // namespace cuda
} // namespace ATV
