/// @file cuda/potentials/pair_compute.cuh
/// @brief Generic per-atom Verlet-list force launcher (force-infra-03).
///
/// The shared CUDA pair traversal + reduction, parameterized on a
/// `ATV::potentials::PairKernel` (LJ / Coulomb-real / combined). One thread per
/// atom iterates its FULL per-atom Verlet list, computes the per-pair force via
/// the kernel, accumulates the owning atom's force in registers (atomic-free —
/// thread i uniquely owns fx[i]), and block-reduces energy + the 6 virial
/// components. The full list double-counts, so energy/virial are halved —
/// identical ownership/structure to the prior hand-written `lj_force_verlet`.
///
/// The kernel is instantiated at `real_t` (the prior float-math path), so the
/// numerics match bit-for-bit. `Kernel::needs_charge` gates the per-pair charge
/// read so charge-free styles (LJ) pay nothing for it.
#pragma once

#include <cstddef>

#include "atomiverse/cuda/box.cuh"
#include "atomiverse/cuda/launch.cuh"
#include "atomiverse/cuda/types.cuh"
#include "atomiverse/potentials/pair_kernel.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief Binary search for `target` in the sorted exclusion partner slice.
__device__ inline bool pc_excl_contains(const int* __restrict__ partners,
                                        int beg, int end, int target) {
    int lo = beg, hi = end;
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        const int v = partners[mid];
        if (v == target)
            return true;
        if (v < target)
            lo = mid + 1;
        else
            hi = mid;
    }
    return false;
}

// Virial path needs ~7 extra accumulators; capping at 6 resident blocks/SM
// starved it of registers (ncu: ~80M local-memory spill requests at 1M atoms).
// Give the virial instantiation more registers (min 4 blocks/SM) while keeping
// the lean no-virial path dense (6 blocks/SM).
template <bool Ortho, bool Virial, class Kernel>
__global__ __launch_bounds__(ATV::cuda::kDefaultBlock, Virial ? 4 : 6) void pair_compute_verlet(
    real_t* __restrict__ fx, real_t* __restrict__ fy, real_t* __restrict__ fz,
    prec_t* __restrict__ energy, prec_t* __restrict__ stress,
    const real_t* __restrict__ x, const real_t* __restrict__ y,
    const real_t* __restrict__ z, const real_t* __restrict__ q,
    const int* __restrict__ vl_nbr, const int* __restrict__ vl_cnt,
    int vl_stride, int n_atoms, cuda::Box box, bool periodic, Kernel kernel,
    work_t cutoff_sq, work_t min_pair_dist_sq,
    const int* __restrict__ excl_offsets, const int* __restrict__ excl_partners,
    work_t s_excl) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    work_t my_energy = work_t(0);
    work_t my_vxx = work_t(0), my_vyy = work_t(0), my_vzz = work_t(0);
    work_t my_vxy = work_t(0), my_vxz = work_t(0), my_vyz = work_t(0);

    const real_t cut_r = real_t(cutoff_sq);
    const real_t minr_r = real_t(min_pair_dist_sq);
    const real_t sexcl_r = real_t(s_excl);

    if (i < n_atoms) {
        const real_t xi = x[i], yi = y[i], zi = z[i];
        real_t fxi = real_t(0), fyi = real_t(0), fzi = real_t(0);
        real_t qi = real_t(0);
        if constexpr (Kernel::needs_charge)
            qi = q[i];
        const bool has_excl = (excl_offsets != nullptr);
        const int beg_excl = has_excl ? excl_offsets[i] : 0;
        const int end_excl = has_excl ? excl_offsets[i + 1] : 0;

        const int cnt = vl_cnt[i];
        for (int k = 0; k < cnt; ++k) {
            const int aj = vl_nbr[static_cast<std::size_t>(k) * vl_stride + i];
            real_t dx = x[aj] - xi;
            real_t dy = y[aj] - yi;
            real_t dz = z[aj] - zi;
            if (periodic) {
                if constexpr (Ortho) {
                    dx -= box.cell[0] * nearbyintf(dx * box.inv_l[0]);
                    dy -= box.cell[4] * nearbyintf(dy * box.inv_l[1]);
                    dz -= box.cell[8] * nearbyintf(dz * box.inv_l[2]);
                } else {
                    box.min_image(dx, dy, dz);
                }
            }
            const real_t r2 = dx * dx + dy * dy + dz * dz;
            if (r2 >= cut_r || r2 < minr_r)
                continue;
            real_t scale = real_t(1);
            if (has_excl && end_excl > beg_excl
                && pc_excl_contains(excl_partners, beg_excl, end_excl, aj)) {
                if (sexcl_r == real_t(0))
                    continue;
                scale = sexcl_r;
            }
            real_t qj = real_t(0);
            if constexpr (Kernel::needs_charge)
                qj = q[aj];

            // raw(): the launcher already pre-skipped r2 >= cutoff_sq, so the
            // leaf kernel need not re-mask (perf-parity with the hand kernel).
            const ATV::potentials::PairContribution<real_t> c =
                kernel.raw(r2, qi, qj);
            const real_t f_over_r = scale * c.fpair;
            fxi += -f_over_r * dx;
            fyi += -f_over_r * dy;
            fzi += -f_over_r * dz;
            my_energy += work_t(0.5) * work_t(scale * c.energy);
            if constexpr (Virial) {
                my_vxx += work_t(0.5) * work_t(dx * (f_over_r * dx));
                my_vyy += work_t(0.5) * work_t(dy * (f_over_r * dy));
                my_vzz += work_t(0.5) * work_t(dz * (f_over_r * dz));
                my_vxy += work_t(0.5) * work_t(dx * (f_over_r * dy));
                my_vxz += work_t(0.5) * work_t(dx * (f_over_r * dz));
                my_vyz += work_t(0.5) * work_t(dy * (f_over_r * dz));
            }
        }
        fx[i] += fxi;
        fy[i] += fyi;
        fz[i] += fzi;
    }

    extern __shared__ work_t warp_scratch[];
    if constexpr (Virial) {
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
    } else {
        work_t my_vals[1] = {my_energy};
        work_t out[1];
        block_reduce_sum_n<1>(my_vals, warp_scratch, out);
        if (threadIdx.x == 0)
            atomicAdd(energy, static_cast<prec_t>(out[0]));
    }
}

/// @brief Host launcher: Ortho/Virial dispatch + geometry, mirroring the prior
///        PairLJCut::compute. `stress != nullptr` selects the virial template.
template <class Kernel>
void launch_pair_compute(const Kernel& kernel, real_t* fx, real_t* fy,
                         real_t* fz, prec_t* energy, prec_t* stress,
                         const real_t* x, const real_t* y, const real_t* z,
                         const real_t* q, const int* vl_nbr, const int* vl_cnt,
                         int vl_stride, int n_atoms, cuda::Box box,
                         bool periodic, work_t cutoff_sq, work_t min_pair_dist_sq,
                         const int* excl_offsets, const int* excl_partners,
                         work_t s_excl, cudaStream_t stream) {
    if (n_atoms <= 0)
        return;
    const int block = kDefaultBlock;
    const int grid = grid_size(n_atoms);
    const int shmem = warps_per_block(block) * sizeof(work_t);
    const bool virial = (stress != nullptr);
    auto go = [&](auto ortho_tag, auto virial_tag) {
        constexpr bool O = decltype(ortho_tag)::value;
        constexpr bool V = decltype(virial_tag)::value;
        pair_compute_verlet<O, V, Kernel><<<grid, block, shmem, stream>>>(
            fx, fy, fz, energy, stress, x, y, z, q, vl_nbr, vl_cnt, vl_stride,
            n_atoms, box, periodic, kernel, cutoff_sq, min_pair_dist_sq,
            excl_offsets, excl_partners, s_excl);
    };
    using T = std::true_type;
    using F = std::false_type;
    if (box.ortho) {
        if (virial)
            go(T{}, T{});
        else
            go(T{}, F{});
    } else {
        if (virial)
            go(F{}, T{});
        else
            go(F{}, F{});
    }
}

} // namespace md
} // namespace cuda
} // namespace ATV
