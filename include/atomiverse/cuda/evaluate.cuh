/// @file cuda/evaluate.cuh
/// @brief cuda::Evaluate — concrete one-pass CUDA evaluator driver.
///
/// `cuda::Evaluate` is the user-facing CUDA outer-driver under the
/// `EvaluateDriver` category. It runs the registered evaluator components
/// exactly once on the current `cpu::Frame` (transferred H2D, evaluated on
/// the device, retrieved D2H) and returns — no time advance, no
/// minimization loop.
///
/// Single-point CUDA SCF energy / density workflows go here:
/// `cuda::Evaluate + cuda::scf::SCFEvaluator`.
///
/// `cuda::Evaluate` is a thin concrete driver — it adds no runtime
/// behaviour beyond the `cuda::EvaluateDriver` base and contributes no
/// infrastructure components in its constructor. It mirrors `cpu::Evaluate`
/// minus the CPU-only categorical setters (`set_pair` / `set_neighbor` …),
/// which the CUDA driver layer does not provide. Components are registered
/// through the inherited public `add()` / `on(phase).add()`.
///
/// Inheritance: Evaluate → cuda::EvaluateDriver →
/// cuda::CUDADriver<evaluate::Phase> → Driver<evaluate::Phase,
/// cuda::device_tag>.
///
/// @code{.cpp}
///   cuda::Evaluate eval;
///   eval.add(cuda::scf::SCFEvaluator{integrals, tol, max_iter});
///   eval.build();
///   eval.run(frame);            // frame.energy populated
/// @endcode
///
/// @see cuda::EvaluateDriver for the inherited one-pass loop + validation
/// @see cuda::single_point for the free-function single-evaluator helper
#pragma once

#include "atomiverse/cuda/evaluate_driver.cuh"

namespace ATV {
namespace cuda {

/// @brief One-pass CUDA evaluator driver.
///
/// Validation rule (inherited from `cuda::EvaluateDriver`): at least one
/// component with `role()==Evaluator` and a non-empty `capabilities()`
/// bitset.
class Evaluate : public EvaluateDriver {};

} // namespace cuda
} // namespace ATV
