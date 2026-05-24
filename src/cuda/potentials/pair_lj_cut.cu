/// @file cuda/potentials/pair_lj_cut.cu
/// @brief CUDA PairLJCut — Lennard-Jones 12-6 shifted (kernel + class).
///        Includes LAMMPS-style `special_bonds` 1-2/1-3/1-4 exclusions.

#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/neighbor_list.cuh"
#include "atomiverse/cuda/potentials/lj_kernel.cuh"
#include "atomiverse/cuda/potentials/pair_lj_cut.cuh"
#include "atomiverse/cuda/types.cuh"
#include <algorithm>
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

__global__ void lj_force(real_t* __restrict__ fx, real_t* __restrict__ fy,
                         real_t* __restrict__ fz, prec_t* __restrict__ energy,
                         prec_t* __restrict__ stress,
                         const int* __restrict__ pi, const int* __restrict__ pj,
                         const real_t* __restrict__ pdx,
                         const real_t* __restrict__ pdy,
                         const real_t* __restrict__ pdz,
                         const real_t* __restrict__ pdist, work_t epsilon,
                         work_t sigma, work_t cutoff, work_t v_shift,
                         const int* __restrict__ n_pairs_ptr,
                         const int* __restrict__ excl_offsets,
                         const int* __restrict__ excl_partners,
                         work_t s_excl) {

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int n_pairs = *n_pairs_ptr;

    work_t my_energy = work_t(0);
    work_t my_vxx = work_t(0), my_vyy = work_t(0), my_vzz = work_t(0);
    work_t my_vxy = work_t(0), my_vxz = work_t(0), my_vyz = work_t(0);
    if (idx < n_pairs) {
        int i = pi[idx], j = pj[idx];
        // Reject ghost / padding pairs emitted by the cluster-based neighbor
        // list export. The cluster representation pads each cluster to
        // kClusterSize atoms; padded slots leak into the flat PairList with
        // j = -1 (sentinel) and garbage `dx/dy/dz` that can coincidentally
        // produce a finite `pdist[idx]` not caught by the `r >= kMinPairDist`
        // floor.
        //
        const bool valid_pair = (i >= 0) && (j >= 0);
        work_t r = valid_pair ? pdist[idx] : work_t(0);
        if (valid_pair && r < cutoff && r >= kMinPairDist) {
            // Special-bonds exclusion check.
            // `excl_offsets == nullptr` ⇒ no exclusions configured ⇒ fast path.
            work_t scale = work_t(1);
            bool skip_pair = false;
            if (excl_offsets != nullptr) {
                int beg = excl_offsets[i];
                int end = excl_offsets[i + 1];
                if (end > beg
                    && excl_contains(excl_partners, beg, end, j)) {
                    if (s_excl == work_t(0))
                        skip_pair = true;
                    else
                        scale = s_excl;
                }
            }

            if (!skip_pair) {
                work_t inv_r = work_t(1) / r;
                work_t sr = sigma * inv_r;
                work_t sr2 = sr * sr;
                work_t sr6 = sr2 * sr2 * sr2;
                work_t sr12 = sr6 * sr6;

                work_t f_over_r = scale * work_t(24) * epsilon
                                  * (work_t(2) * sr12 - sr6) * inv_r * inv_r;
                work_t dx = pdx[idx], dy = pdy[idx], dz = pdz[idx];

                atomicAdd(&fx[i], static_cast<real_t>(-f_over_r * dx));
                atomicAdd(&fy[i], static_cast<real_t>(-f_over_r * dy));
                atomicAdd(&fz[i], static_cast<real_t>(-f_over_r * dz));
                atomicAdd(&fx[j], static_cast<real_t>(f_over_r * dx));
                atomicAdd(&fy[j], static_cast<real_t>(f_over_r * dy));
                atomicAdd(&fz[j], static_cast<real_t>(f_over_r * dz));

                my_energy =
                    scale * (work_t(4) * epsilon * (sr12 - sr6) - v_shift);

                // Pair virial (LAMMPS sign: w_ab = Σ_pairs r_a · F_b, where
                // F on atom j = +f_over_r · d_ab). The diagonal of the
                // stress tensor (used by NoseHooverNPT) is the sum of pair
                // virials plus the kinetic component, divided by 3·V.
                my_vxx = dx * (f_over_r * dx);
                my_vyy = dy * (f_over_r * dy);
                my_vzz = dz * (f_over_r * dz);
                my_vxy = dx * (f_over_r * dy);
                my_vxz = dx * (f_over_r * dz);
                my_vyz = dy * (f_over_r * dz);
            }
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
}

void PairLJCut::compute(const LaunchContext& ctx) {
    auto stream = static_cast<cudaStream_t>(ctx.runtime.stream);
    auto& frame = ctx.frame;
    const PairList* pairs = ctx.runtime.pair_list;
    if (!pairs)
        return;
    int n = ctx.frame.n_atoms();
    (void)n; // kept for symmetry with sibling pair styles
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
    if (pair_capacity_host_ == 0) {
        cudaMemcpyAsync(&pair_capacity_host_, pairs->capacity, sizeof(int),
                        cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream); // build-time-equivalent: first launch only
    }
    int block = kDefaultBlock;
    int grid = grid_size(pair_capacity_host_ > 0 ? pair_capacity_host_ : 1);
    int shmem = warps_per_block(block) * sizeof(work_t);
    // Pass nullptr for exclusion CSR when none configured — kernel takes
    // fast path and avoids the partner lookup overhead entirely.
    const int* excl_off = (has_exclusions_) ? d_excl_offsets_ : nullptr;
    const int* excl_par = (has_exclusions_) ? d_excl_partners_ : nullptr;
    kernels::lj_force<<<grid, block, shmem, stream>>>(
        frame.atoms.fx, frame.atoms.fy, frame.atoms.fz, &frame.scalars->energy,
        &frame.scalars->stress[0], pairs->i, pairs->j, pairs->dx, pairs->dy,
        pairs->dz, pairs->dist, epsilon_, sigma_, cutoff_, v_shift_,
        pairs->n_pairs, excl_off, excl_par, static_cast<work_t>(s_excl_));
}

} // namespace cuda
} // namespace ATV
