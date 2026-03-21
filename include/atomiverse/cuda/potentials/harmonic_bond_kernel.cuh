/// @file cuda/potentials/harmonic_bond_kernel.cuh
/// @brief Pointer-level harmonic bond kernel -- shared between cuda component and torch backend.
///
/// V(r) = K (r - r0)^2,  F = -2K(r - r0) r_hat
///
/// Single source of truth for harmonic bond force+energy computation.
#pragma once

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace cuda {
namespace kernels {

/// @brief Harmonic bond stretching kernel. One thread per bond.
///
/// Forces accumulated via atomicAdd (Newton's 3rd law).
/// Energy reduced per-block using shared memory, then atomicAdd to global scalar.
/// Uses minimum image convention (orthorhombic box).
///
/// @param x,y,z           Position arrays [n_atoms, device].
/// @param fx,fy,fz        Force arrays [n_atoms, device]. Accumulated via atomicAdd.
/// @param energy          Potential energy scalar [1, device]. Accumulated via atomicAdd.
/// @param bi,bj,bt        Bond atom indices and type [n_bonds, device].
/// @param K               Spring constants per bond type [n_types, device].
/// @param r0              Equilibrium distances per bond type [n_types, device].
/// @param box             Simulation cell 3x3 row-major [9, device].
/// @param n_bonds         Number of bonds.
__global__ void harmonic_bond_k(
    const Real* __restrict__ x, const Real* __restrict__ y, const Real* __restrict__ z,
    Real* __restrict__ fx, Real* __restrict__ fy, Real* __restrict__ fz,
    Real* __restrict__ energy,
    const Int* __restrict__ bi, const Int* __restrict__ bj, const Int* __restrict__ bt,
    const Real* __restrict__ K, const Real* __restrict__ r0,
    const Real* __restrict__ box, Int n_bonds);

} // namespace kernels
} // namespace cuda
} // namespace ATV
