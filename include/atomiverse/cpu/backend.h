/// @file cpu/backend.h
/// @brief cpu::Backend -- CPU execution engine and component lifecycle manager.
///
/// Operates directly on the user's cpu::Frame (zero-copy, zero transfer).
/// Manages OpenMP thread configuration and MPI awareness.
///
/// @see Driver for the user-facing facade that delegates to Backend
#pragma once

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/facade/plan.h"

namespace ATV {
namespace cpu {

/// @brief CPU execution engine -- lifecycle manager + dispatcher.
///
/// Key difference from CUDA Backend: no materialize/retrieve.
/// CPU Backend operates directly on the user's cpu::Frame (zero-copy).
///
/// Lifecycle:
/// 1. init() -- configure OpenMP, detect MPI
/// 2. set_frame() -- store pointer (zero-copy)
/// 3. build_components() -- one-time resource allocation
/// 4. run() / step() -- per-step dispatch
/// 5. destroy_components() -- release resources (reverse order)
/// 6. shutdown() -- cleanup
class Backend {
public:
    /// @brief Configure OpenMP threads and detect MPI environment.
    void init();

    /// @brief Bind a compiled Plan for execution.
    /// @param plan Plan to execute (non-owning reference, must outlive Backend usage).
    void load(const Plan<cpu::ComponentBase>& plan) { plan_ = &plan; }

    /// @brief Set the Frame to operate on (zero-copy, non-owning).
    /// @param frame Host-side system state.
    void set_frame(cpu::Frame& frame) { frame_ = &frame; }

    /// @brief Build all components in Plan order.
    void build_components();

    /// @brief Execute one simulation step.
    /// @param step_idx Zero-based step index.
    void step(int step_idx);

    /// @brief Execute multiple simulation steps.
    /// @param n_steps Number of steps to execute.
    void run(int n_steps);

    /// @brief Destroy all components in reverse Plan order.
    void destroy_components();

    /// @brief Release all resources.
    void shutdown();

    /// @brief Access the CPU context.
    Context& context() { return context_; }

    /// @brief Access the Frame.
    Frame& frame() { return *frame_; }

private:
    Context context_{};
    cpu::Frame* frame_ = nullptr;
    const Plan<cpu::ComponentBase>* plan_ = nullptr;
};

} // namespace cpu
} // namespace ATV
