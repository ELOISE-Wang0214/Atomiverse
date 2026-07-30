/// @file torch/frame/molrs_exchange.cuh
/// @brief torch::from_molrs / torch::to_molrs -- molrs.Frame <-> torch::Frame.
///
/// `molrs.Frame` is the canonical exchange layer; `torch::Frame` is the
/// torch backend's differentiable execution frame. The materializer is
/// shared infrastructure (sibling of `torch/transfer.cuh` at the backend
/// root, per the CLAUDE.md backend-layering rule) and is meant to run once
/// per `run()` -- never inside a step loop.
///
/// Implemented in terms of the already-tested CPU materializer: a transient
/// `cpu::Frame` is the staging point, and the host<->device hop is the
/// explicit `torch::from_cpu` / `torch::to_cpu` transfer.
///
/// @see torch::Frame for the torch execution frame
/// @see cpu::from_molrs / cpu::to_molrs for the CPU staging materializer
/// @see torch::from_cpu / torch::to_cpu for the host<->device transfer
/// @see molrs.hpp for the molrs CXX FFI bridge
#pragma once

#ifdef ATV_HAS_TORCH

#include "atomiverse/molrs.hpp"
#include "atomiverse/torch/frame/frame.cuh"

namespace ATV {
namespace torch {

/// @brief Build a .Frame into a new torch::Frame (on CUDA).
///
/// Stages through `cpu::from_molrs` (which reads the `atoms` block, the
/// `bonds` / `angles` / `dihedrals` topology blocks, and the box) and then
/// performs an explicit H2D transfer via `torch::from_cpu`. The resulting
/// torch::Frame holds all available fields as CUDA Tensors.
///
/// Required fields are the atom positions (`x`, `y`, `z`) and `type`; their
/// absence throws `std::runtime_error` naming the missing field. Every other
/// field is optional and defaults to zero / empty when absent.
///
/// @param src molrs.Frame to read (introspected by block/column name).
/// @return Newly constructed torch::Frame with data on the CUDA device.
/// @throws std::runtime_error if a required field is missing.
Frame from_molrs(const molrs::FrameRef& src);

/// @brief Write a torch::Frame's results back into an existing molrs.Frame.
///
/// Stages through an explicit D2H transfer (`torch::to_cpu`) into a
/// transient `cpu::Frame`, then delegates to `cpu::to_molrs`, which
/// overwrites the `atoms` block result columns (`x, y, z, vx, vy, vz, fx,
/// fy, fz, charge`) and the system scalars. Topology and box are left
/// untouched.
///
/// @param src torch::Frame whose state is exported.
/// @param dst molrs.Frame to update in place (must be the same frame that
///            was materialized, or any frame with a compatible atom count).
void to_molrs(const Frame& src, molrs::FrameRef& dst);

} // namespace torch
} // namespace ATV

#endif // ATV_HAS_TORCH
