/// @file cuda/potentials/lj_kernel.cuh
/// @brief Pointer-level LJ kernel -- shared between cuda::LJ component and torch backend.
///
/// This kernel is the single source of truth for LJ force+energy computation.
/// It accepts raw pointers so it can be called from any backend.
#pragma once


#include "atomiverse/cuda/types.cuh"

namespace ATV {
namespace cuda {
namespace kernels {

/// @brief Minimum pair distance to avoid LJ singularity [Angstrom].
inline constexpr cuda::work_t kMinPairDist = cuda::work_t(1e-10);

/// @brief LJ 12-6 force + energy kernel. One thread per pair.
///
/// Forces are accumulated via atomicAdd (Newton's 3rd law).
/// Energy is reduced per-block using shared memory, then atomicAdd to global scalar.
///
/// @param fx,fy,fz        Force arrays [n_atoms, device]. Accumulated via atomicAdd.
/// @param energy          Potential energy scalar [1, device]. Accumulated via atomicAdd.
/// @param pi,pj           Pair atom indices [n_pairs, device].
/// @param pdx,pdy,pdz     Pair displacements [n_pairs, device].
/// @param pdist           Pair distances [n_pairs, device].
/// @param epsilon         LJ well depth [eV].
/// @param sigma           LJ collision diameter [Angstrom].
/// @param cutoff          Interaction cutoff [Angstrom].
/// @param v_shift         Potential shift at cutoff [eV].
/// @param n_pairs_ptr     Pair count [1, device].
__global__ void lj_force(
    cuda::real_t* __restrict__ fx, cuda::real_t* __restrict__ fy, cuda::real_t* __restrict__ fz,
    cuda::prec_t* __restrict__ energy,
    const int* __restrict__ pi, const int* __restrict__ pj,
    const cuda::real_t* __restrict__ pdx, const cuda::real_t* __restrict__ pdy,
    const cuda::real_t* __restrict__ pdz, const cuda::real_t* __restrict__ pdist,
    cuda::work_t epsilon, cuda::work_t sigma, cuda::work_t cutoff, cuda::work_t v_shift,
    const int* __restrict__ n_pairs_ptr);

/// @brief Compute LJ potential shift at cutoff.
/// @param epsilon Well depth [eV].
/// @param sigma Collision diameter [Angstrom].
/// @param cutoff Interaction cutoff [Angstrom].
/// @return Shift value V(cutoff).
inline double lj_shift(double epsilon, double sigma, double cutoff) {
    double sr = sigma / cutoff;
    double sr2 = sr * sr;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;
    return 4.0 * epsilon * (sr12 - sr6);
}

} // namespace kernels
} // namespace cuda
} // namespace ATV
