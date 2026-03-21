/// @file cuda/potentials/lj_kernel.cuh
/// @brief Pointer-level LJ kernel -- shared between cuda::LJ component and torch backend.
///
/// This kernel is the single source of truth for LJ force+energy computation.
/// It accepts raw pointers so it can be called from any backend.
#pragma once

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cuda {
namespace kernels {

/// @brief Minimum pair distance to avoid LJ singularity [Angstrom].
inline constexpr Real kMinPairDist = Real(1e-10);

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
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ pi, const Int* __restrict__ pj,
    const Real* __restrict__ pdx, const Real* __restrict__ pdy,
    const Real* __restrict__ pdz, const Real* __restrict__ pdist,
    Real epsilon, Real sigma, Real cutoff, Real v_shift,
    const Int* __restrict__ n_pairs_ptr);

/// @brief Compute LJ potential shift at cutoff.
/// @param epsilon Well depth [eV].
/// @param sigma Collision diameter [Angstrom].
/// @param cutoff Interaction cutoff [Angstrom].
/// @return Shift value V(cutoff).
inline Real lj_shift(Real epsilon, Real sigma, Real cutoff) {
    Real sr = sigma / cutoff;
    Real sr2 = sr * sr;
    Real sr6 = sr2 * sr2 * sr2;
    Real sr12 = sr6 * sr6;
    return Real(4) * epsilon * (sr12 - sr6);
}

} // namespace kernels
} // namespace cuda
} // namespace ATV
