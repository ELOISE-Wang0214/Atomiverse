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
    if (idx < n_pairs) {
        int i = pi[idx], j = pj[idx];
        // Reject ghost / padding pairs emitted by the cluster-based neighbor
        // list export. The cluster representation pads each cluster to
        // kClusterSize atoms; padded slots leak into the flat PairList with
        // j = -1 (sentinel) and garbage `dx/dy/dz` that can coincidentally
        // produce a finite `pdist[idx]` not caught by the `r >= kMinPairDist`
        // floor. Without this guard those ghost pairs bypass the
        // exclusion-CSR mask (-1 is in no partner list) and erroneously
        // contribute LJ energy. Both i and j are checked for symmetry.
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
            }
        }
    }

    extern __shared__ work_t s_energy[];
    s_energy[threadIdx.x] = my_energy;
    __syncthreads();
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s)
            s_energy[threadIdx.x] += s_energy[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0 && s_energy[0] != work_t(0))
        atomicAdd(energy, static_cast<prec_t>(s_energy[0]));
}

} // namespace kernels
} // namespace cuda

namespace {
constexpr int kEstMaxNeighbors = 100;
} // anonymous namespace

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
    int max_pairs = static_cast<int>(n) * kEstMaxNeighbors / 2;
    int block = kDefaultBlock;
    int grid = grid_size(max_pairs);
    // Pass nullptr for exclusion CSR when none configured — kernel takes
    // fast path and avoids the partner lookup overhead entirely.
    const int* excl_off = (has_exclusions_) ? d_excl_offsets_ : nullptr;
    const int* excl_par = (has_exclusions_) ? d_excl_partners_ : nullptr;
    kernels::lj_force<<<grid, block, block * sizeof(work_t), stream>>>(
        frame.atoms.fx, frame.atoms.fy, frame.atoms.fz, &frame.scalars->energy,
        pairs->i, pairs->j, pairs->dx, pairs->dy, pairs->dz, pairs->dist,
        epsilon_, sigma_, cutoff_, v_shift_, pairs->n_pairs, excl_off, excl_par,
        static_cast<work_t>(s_excl_));
}

} // namespace cuda
} // namespace ATV
