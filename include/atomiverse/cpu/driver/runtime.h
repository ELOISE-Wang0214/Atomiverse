/// @file cpu/driver/runtime.h
/// @brief cpu::Runtime -- the CPU loop runtime: Plan dispatch + component
///        lifecycle, composed from communication / device / asynchronous
///        sub-runtimes.
///
/// backend-mirror-05-dissolve replaced the former god-class: domain
/// decomposition now lives in `communication::DomainRuntime`, thread + timing
/// setup in `device::DeviceRuntime`, the async pipeline in
/// `asynchronous::AnalysisRuntime`. This type owns those three, binds the Frame,
/// and dispatches Plans (one per driver node in the tree). It operates directly
/// on the user's cpu::Frame (zero-copy, zero transfer).
///
/// One Runtime is shared by a whole driver tree: a sub-driver runs its own Plan
/// on its parent's Runtime (published on `Context::exec.runtime`), so nesting
/// does not double the runtime, the decomposition, or the frame binding.
#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atomiverse/core/plan.h"
#include "atomiverse/core/stage.h"
#include "atomiverse/core/timer.h"
#include "atomiverse/cpu/asynchronous/analysis_runtime.h"
#include "atomiverse/cpu/communication/domain_runtime.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/device/device_runtime.h"
#include "atomiverse/cpu/frame/frame.h"

namespace ATV {
namespace cpu {

namespace md {
namespace mpi {
struct DomainDecomposition;
} // namespace mpi
} // namespace md

/// @brief CPU loop runtime -- lifecycle manager + Plan dispatcher.
class Runtime {
  public:
    Runtime();
    ~Runtime();
    Runtime(const Runtime&) = delete;
    Runtime& operator=(const Runtime&) = delete;
    Runtime(Runtime&&) noexcept;
    Runtime& operator=(Runtime&&) noexcept;

    /// @brief Configure threads (device), detect MPI, allocate the DD
    ///        (communication), start the analysis pipeline, publish services.
    void init();

    /// @brief Bind a compiled Plan for execution (non-owning).
    void load(const Plan<cpu::Component>& plan) { plan_ = &plan; }

    /// @brief Set the Frame to operate on (zero-copy) and bind the DD to it.
    void set_frame(cpu::Frame& frame);

    /// @brief Build all components in Plan order.
    void build_components();

    /// @brief Execute one simulation step.
    void step(int step_idx);

    // ── Nested plans: a sub-driver runs its own Plan on THIS runtime ──────────
    void build_plan(const Plan<cpu::Component>& plan);
    void step_plan(const Plan<cpu::Component>& plan, int step_idx);
    void destroy_plan(const Plan<cpu::Component>& plan);

    /// @brief Re-run only the bindings of a single phase (line-search trials).
    /// @tparam PhaseEnum The workflow Phase enum this runtime's Plan was compiled from.
    template <typename PhaseEnum>
    void run_stage(PhaseEnum phase, int step_idx) {
        run_slot_(ATV::slot_of(phase), step_idx);
    }

    /// @brief Execute multiple simulation steps.
    void run(int n_steps);

    /// @brief Destroy all components in reverse Plan order.
    void destroy_components();

    /// @brief Block until the async-analysis executor has drained.
    void drain_analysis() { analysis_.drain(); }

    /// @brief Release all resources.
    void shutdown();

    /// @brief Access the CPU context.
    Context& context() { return context_; }

    /// @brief Cap the OpenMP thread count. Must be called before `init()`.
    void set_threads(int n) { device_.set_threads(context_, n); }

    /// @brief Access the Frame.
    Frame& frame() { return *frame_; }

    /// @brief The per-component timing registry (opt-in), owned by the device runtime.
    TimerRegistry& timing() { return device_.timing(); }
    const TimerRegistry& timing() const { return device_.timing(); }

    /// @brief Per-(unique-)component owned heap bytes for the driver memory summary.
    std::vector<std::pair<std::string, std::size_t>>
    collect_component_memory() const;

    /// @brief Enable/disable per-launch timing (device runtime; preference sticks).
    void enable_timing(bool on) { device_.enable_timing(on); }

    /// @brief Clear accumulated step timing.
    void reset_timing() { device_.reset_timing(); }

    /// @brief Adopt an externally-owned DomainDecomposition (before `init()`).
    /// Out-of-line: the unique_ptr is destroyed where DomainDecomposition is complete.
    void adopt_dd(std::unique_ptr<md::mpi::DomainDecomposition> dd);

    /// @brief Release the owned DomainDecomposition for reuse.
    std::unique_ptr<md::mpi::DomainDecomposition> release_dd();

    /// @brief Run serially on a full-replica frame even under multi-rank MPI.
    void set_serial_on_replica() noexcept { domain_.set_serial_on_replica(); }

  private:
    /// @brief Replay one compiled Plan slot. The erased half of `run_stage`.
    void run_slot_(int slot, int step_idx);

    Context context_{};
    cpu::Frame* frame_ = nullptr;
    const Plan<cpu::Component>* plan_ = nullptr;

    communication::DomainRuntime domain_{};
    device::DeviceRuntime device_{};
    asynchronous::AnalysisRuntime analysis_{};

    /// Per-slot, per-op dense timing handles, precomputed in `build_components()`
    /// so the hot `step()` loop never calls `component->name()`.
    std::vector<std::vector<int>> op_handle_;
};

} // namespace cpu
} // namespace ATV
