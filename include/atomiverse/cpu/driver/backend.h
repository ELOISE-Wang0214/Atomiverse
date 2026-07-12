/// @file cpu/driver/backend.h
/// @brief cpu::Backend -- CPU execution engine and component lifecycle manager.
///
/// Operates directly on the user's cpu::Frame (zero-copy, zero transfer).
/// Manages OpenMP thread configuration, MPI awareness, and the MD spatial
/// `DomainDecomposition` published through `Context::dd`.
///
/// @see Driver for the user-facing facade that delegates to Backend
#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atomiverse/core/timer.h"
#include "atomiverse/cpu/component/component.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
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
    /// DomainDecomposition, and publish it through `context().md.dd`.
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

    /// @brief Free `atoms.map` and `dd->tag_to_local` when the frame has
    /// no bonded topology (nothing resolves global tags at run time —
    /// they cost ~50 B/atom each when populated). Called from
    /// `set_frame()` after the DD bind/partition.
    void release_tag_maps_if_unused(cpu::Frame& frame);

    /// @brief Build all components in Plan order.
    void build_components();

    /// @brief Execute one simulation step.
    /// @param step_idx Zero-based step index.
    void step(int step_idx);

    /// @brief Re-run only the bindings of a single phase on the current frame.
    ///
    /// Used for line-search trial evaluations: an Optimizer component, mid
    /// `launch()`, writes a trial geometry into the frame and calls this to
    /// recompute energy/forces by replaying the force phase (ZeroForces →
    /// evaluator → MPI reduce). Every op in matching slots is launched
    /// unconditionally (frequency gating is bypassed — a trial is not a step)
    /// and timing is not recorded (this runs re-entrantly inside `step()`).
    /// @param phase_value Phase index to replay (e.g. minimizer force phase).
    /// @param step_idx    Step index forwarded to the LaunchContext.
    void run_phase(int phase_value, int step_idx);

    /// @brief Execute multiple simulation steps.
    /// @param n_steps Number of steps to execute.
    void run(int n_steps);

    /// @brief Destroy all components in reverse Plan order.
    void destroy_components();

    /// @brief Block until the async-analysis executor has drained every
    /// queued and in-flight job, leaving its workers alive.
    ///
    /// A continuable driver calls this at the end of each `run()` so all
    /// snapshots from that segment reach their subscribers before
    /// control returns. No-op when no executor exists.
    void drain_analysis();

    /// @brief Release all resources.
    void shutdown();

    /// @brief Access the CPU context.
    Context& context() { return context_; }

    /// @brief Cap the OpenMP thread count. Must be called before `init()`.
    ///
    /// `init()` already honours `context_.n_threads` when it is positive (it calls
    /// `omp_set_num_threads`) — the plumbing was there all along; there was simply no
    /// public way in, so an in-process caller with a small system had no answer to
    /// OpenMP fork/join overhead except the OMP_NUM_THREADS env var.
    /// @param n Thread count (>= 1).
    void set_threads(int n) { context_.n_threads = (n < 1) ? 0 : n; }

    /// @brief Access the Frame.
    Frame& frame() { return *frame_; }

    /// @brief The per-component timing registry (opt-in). Records each
    /// component launch when enabled; reached by components via
    /// `Context::timer`. @see core/timer.h, cpu/md/timing_report.h
    TimerRegistry& timing() { return timing_; }
    const TimerRegistry& timing() const { return timing_; }

    /// @brief Per-(unique-)component owned heap bytes for the built-in driver
    /// memory summary: each registered component's `memory_bytes()` keyed by
    /// `name()`, plus synthetic entries `[atoms]` (per-atom SoA + tag map) and
    /// `[domain]` (decomposition maps). The MD layer classifies the names into
    /// categories and reduces across ranks. @see cpu/md/timing_report.h
    std::vector<std::pair<std::string, std::size_t>>
    collect_component_memory() const;

    /// @brief Enable/disable per-launch timing. Default ON (so every driver
    /// emits the built-in end-of-run summary). The preference is remembered and
    /// re-applied by `init()`, so calling this BEFORE the first `run()` sticks;
    /// the `ATV_TIMER` env var (`0`/`off`) still overrides in `init()`.
    void enable_timing(bool on) {
        timing_pref_ = on;
        timing_.set_enabled(on);
    }

    /// @brief Clear accumulated step timing (e.g. between protocol stages).
    void reset_timing() { timing_.reset(); }

    /// @brief Adopt an externally-owned `DomainDecomposition`.
    ///
    /// Takes ownership of `dd`. Must be called BEFORE `init()`; once the
    /// dd is set, `init()` skips the default `make_unique<DD>()` step and
    /// rebuilds topology on the adopted instance instead (preserving its
    /// `n_local`, `tag_to_local`, `sublo/subhi`, and `cuts`, while
    /// recreating `cart_comm` on the current `MPI_COMM_WORLD`).
    ///
    /// Intended for stage-orchestrators that create one Backend per
    /// stage but want the spatial partition state to survive the
    /// teardown — without this, each fresh Backend's
    /// `partition_atoms(frame)` re-iterates the whole atom array
    /// (`first_partition = (n_local == 0)`) and misinterprets the
    /// previous stage's stale tail as owned atoms.
    void adopt_dd(std::unique_ptr<md::mpi::DomainDecomposition> dd);

    /// @brief Release the owned `DomainDecomposition` for reuse.
    ///
    /// Call between stages to extract the DD before tearing this Backend
    /// down; pass it back to the next stage's Backend via `adopt_dd()`.
    /// Returns nullptr if no DD is owned.
    std::unique_ptr<md::mpi::DomainDecomposition> release_dd();

  private:
    Context context_{};
    cpu::Frame* frame_ = nullptr;
    const Plan<cpu::ComponentBase>* plan_ = nullptr;
    std::unique_ptr<md::mpi::DomainDecomposition> dd_;
    std::unique_ptr<AnalysisExecutor> executor_;

    TimerRegistry timing_{};
    /// @brief Desired timing-enabled state, re-applied by `init()`. Default
    /// true (built-in summary always-on); `enable_timing()` updates it.
    bool timing_pref_ = true;
    /// Per-slot, per-op dense timing handles into `timing_`, precomputed in
    /// `build_components()` so the hot `step()` loop never calls
    /// `component->name()`. Mirrors `plan_->slots[s].ops` order.
    std::vector<std::vector<int>> op_handle_;
};

} // namespace cpu
} // namespace ATV
