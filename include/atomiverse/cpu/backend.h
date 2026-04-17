/// @file cpu/backend.h
/// @brief cpu::Backend -- CPU execution engine and component lifecycle manager.
///
/// Operates directly on the user's cpu::Frame (zero-copy, zero transfer).
/// Manages OpenMP thread configuration, MPI awareness, and the MD spatial
/// `DomainDecomposition` published through `Context::dd`.
///
/// @see Driver for the user-facing facade that delegates to Backend
#pragma once

#include <memory>

#include "atomiverse/cpu/component.h"
#include "atomiverse/cpu/frame.h"
#include "atomiverse/cpu/context.h"
#include "atomiverse/facade/plan.h"

namespace ATV {
namespace cpu {

class AnalysisExecutor;

namespace md {
namespace mpi {
struct DomainDecomposition;
} // namespace mpi
} // namespace md

/// @brief CPU execution engine -- lifecycle manager + dispatcher.
///
/// Key difference from CUDA Backend: no materialize/retrieve.
/// CPU Backend operates directly on the user's cpu::Frame (zero-copy).
///
/// Lifecycle:
/// 1. init() -- configure OpenMP, detect MPI, allocate single-rank DomainDecomposition
/// 2. set_frame() -- store pointer (zero-copy), bind DomainDecomposition to the frame
/// 3. build_components() -- one-time resource allocation
/// 4. run() / step() -- per-step dispatch
/// 5. destroy_components() -- release resources (reverse order)
/// 6. shutdown() -- cleanup
class Backend {
public:
    /// @brief Default-construct an uninitialized Backend.
    Backend();

    /// @brief Release owned resources (DomainDecomposition).
    ~Backend();

    Backend(const Backend&) = delete;
    Backend& operator=(const Backend&) = delete;
    Backend(Backend&&) noexcept;
    Backend& operator=(Backend&&) noexcept;

    /// @brief Configure OpenMP threads, detect MPI, allocate single-rank
    /// DomainDecomposition, and publish it through `context().dd`.
    void init();

    /// @brief Bind a compiled Plan for execution.
    /// @param plan Plan to execute (non-owning reference, must outlive Backend usage).
    void load(const Plan<cpu::ComponentBase>& plan) { plan_ = &plan; }

    /// @brief Set the Frame to operate on (zero-copy, non-owning).
    ///
    /// Also binds the owned `DomainDecomposition` to the frame: populates
    /// `dd->n_local`, `dd->capacity`, and `dd->tag_to_local`. Safe to call
    /// multiple times (rebinds to the new frame).
    /// @param frame Host-side system state.
    void set_frame(cpu::Frame& frame);

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
    std::unique_ptr<md::mpi::DomainDecomposition> dd_;
    std::unique_ptr<AnalysisExecutor> executor_;
};

} // namespace cpu
} // namespace ATV
