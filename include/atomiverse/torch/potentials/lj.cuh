/// @file torch/potentials/lj.cuh
/// @brief Differentiable LJ 12-6 operator for the torch backend.
///
/// Forward reuses cuda::kernels::lj_force. Backward provides analytical
/// gradients w.r.t. positions, epsilon, and sigma via autograd::Function.
///
/// @see cuda::kernels::lj_force for the shared kernel
/// @see torch::Frame for the Tensor-based system state
#pragma once

#ifdef ATV_HAS_TORCH

#include "atomiverse/torch/frame.cuh"

namespace ATV {
namespace torch {

/// @brief Compute differentiable LJ energy and forces.
///
/// Accumulates energy into frame.energy and forces into frame.fx/fy/fz.
/// Supports autograd w.r.t. positions (frame.x/y/z), epsilon, and sigma.
///
/// @param frame     torch::Frame with positions and populated pair list.
/// @param epsilon   LJ well depth [eV]. Scalar Tensor, may require grad.
/// @param sigma     LJ collision diameter [Angstrom]. Scalar Tensor, may require grad.
/// @param cutoff    Interaction cutoff [Angstrom]. Not differentiable.
/// @return Energy Tensor [1] (for .backward() calls).
Tensor lj(Frame& frame, const Tensor& epsilon, const Tensor& sigma, double cutoff);

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
