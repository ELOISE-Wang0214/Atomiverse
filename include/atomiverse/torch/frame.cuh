/// @file torch/frame.cuh
/// @brief torch::Frame -- differentiable system state backed by at::Tensor.
///
/// Parallel to cpu::Frame (xtensor) and cuda::Frame (raw device pointers).
/// Tensors may have requires_grad=true for autograd-enabled workflows.
///
/// @see cpu::Frame for the host-side counterpart
/// @see cuda::Frame for the raw-pointer CUDA counterpart
/// @see torch::to_torch() / torch::to_cuda() for transfer
#pragma once

#ifdef ATV_HAS_TORCH

#include <torch/torch.h>

#include "atomiverse/core/int.h"
#include "atomiverse/core/real.h"

namespace ATV {
namespace torch {

/// @brief LibTorch type aliases to resolve ATV::torch vs ::torch collision.
/// @{
using Tensor = ::torch::Tensor;
using ScalarType = ::at::ScalarType;
using Device = ::torch::Device;
/// @}

/// @brief Compile-time mapping: ATV Real -> Torch ScalarType.
inline constexpr ScalarType kRealDtype =
    std::is_same_v<Real, float> ? ScalarType::Float : ScalarType::Double;

/// @brief Compile-time mapping: ATV Int -> Torch ScalarType.
inline constexpr ScalarType kIntDtype =
    std::is_same_v<Int, int> ? ScalarType::Int : ScalarType::Long;

/// @brief Differentiable system state backed by at::Tensor.
///
/// Each field is a CUDA Tensor in SoA layout. Positions and parameters
/// may carry requires_grad=true for autograd; topology and indices do not.
///
/// Usage:
/// @code{.cpp}
/// auto tframe = ATV::torch::to_torch(cuda_frame);
/// auto energy = ATV::torch::lj(tframe, epsilon, sigma, cutoff);
/// energy.backward();
/// @endcode
struct Frame {
    // ---- Per-atom tensors [N], dtype=kRealDtype, device=CUDA ----

    Tensor x, y, z;              ///< Positions [Angstrom].
    Tensor vx, vy, vz;           ///< Velocities [Angstrom/fs].
    Tensor fx, fy, fz;           ///< Forces [eV/Angstrom].
    Tensor mass;                 ///< Masses [amu].
    Tensor type_id;              ///< Atom type indices (kIntDtype, no grad).

    // ---- Topology (kIntDtype, no grad) ----

    Tensor bond_i, bond_j, bond_type;              ///< Bond topology [n_bonds].
    Tensor angle_i, angle_j, angle_k, angle_type;  ///< Angle topology [n_angles], j is center.
    Tensor dihedral_i, dihedral_j, dihedral_k, dihedral_l, dihedral_type; ///< Dihedral topology [n_dihedrals].

    // ---- Pair list (populated by neighbor_list, no grad on indices) ----

    Tensor pair_i, pair_j;                   ///< Pair indices [n_pairs], kIntDtype.
    Tensor pair_dx, pair_dy, pair_dz;        ///< Pair displacements [n_pairs], kRealDtype.
    Tensor pair_dist;                        ///< Pair distances [n_pairs], kRealDtype.

    // ---- System scalars ----

    Tensor energy;               ///< Potential energy [1], [eV].
    Tensor kinetic;              ///< Kinetic energy [1], [eV].
    Tensor stress;               ///< Stress tensor [6], Voigt [eV/Angstrom^3].
    Tensor box;                  ///< Simulation cell [9], row-major [Angstrom].

    // ---- Metadata (host-cached) ----

    /// @brief Number of atoms in the system.
    Int n_atoms() const { return n_atoms_; }

    /// @brief Number of bonds.
    Int n_bonds() const { return n_bonds_; }

    /// @brief Number of pairs in the pair list.
    Int n_pairs() const { return n_pairs_; }

    Int n_atoms_ = 0;      ///< Host-cached atom count.
    Int n_bonds_ = 0;      ///< Host-cached bond count.
    Int n_angles_ = 0;     ///< Host-cached angle count.
    Int n_dihedrals_ = 0;  ///< Host-cached dihedral count.
    Int n_pairs_ = 0;      ///< Host-cached pair count.
};

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
