/// @file cpu/task/task_materializer.h
/// @brief TaskMaterializer -- rank-agnostic TaskSpec -> runnable driver
///        (mpi-launcher-embedded-python-02).
///
/// The constructor materializes the driver a TaskSpec describes and
/// establishes the run-ready invariant (registration replay + `build()`):
///
///  - `TaskKind::MD` -- a `cpu::md::MD`, rebuilt by replaying each entry
///    exactly the way the Python bindings registered it (`md.add(comp)` for
///    default-binding entries, `md.on(phase).add(comp, order)` for explicit
///    ones -- never the categorical setters), then `build()`.
///  - `TaskKind::SinglePoint` -- a `cpu::Evaluate` carrying an
///    `AM1Evaluator` bound at `evaluate::Phase::evaluate`, then `build()`.
///
/// `run(frame)` is the ONE task-kind-agnostic consumption entry (chain 03's
/// WorkerLoop needs no per-kind branch): MD dispatches its step budget,
/// SinglePoint dispatches one Evaluate pass.
///
/// Single-point results contract: the AM1 pass writes the total energy [eV]
/// onto `frame.scalars.energy` and per-atom Mulliken charges [e] onto
/// `frame.atoms.charge` -- both fields are part of the
/// `comm::pack_exchange` schema, so they ride the existing Frame gather
/// back to rank 0 with no extra channel.
///
/// @see cpu/task/task_spec.h
#pragma once

#include <cstdint>
#include <optional>

#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/md/md.h"
#include "atomiverse/cpu/task/task_spec.h"

namespace ATV {
namespace cpu {
namespace task {

/// @brief Materializes a TaskSpec into a run-ready driver pipeline.
class TaskMaterializer {
  public:
    /// @brief Replay the spec's registrations and build the driver.
    /// @param spec Validated task description (see TaskSpec's ctors).
    explicit TaskMaterializer(const TaskSpec& spec);

    /// @brief Execute the task on `frame` (task-kind-agnostic).
    /// @param frame Host-resident system state, updated in place
    ///              (positions [Angstrom], energies [eV], charges [e]).
    void run(cpu::Frame& frame);

    /// @brief The kind of task this materializer runs.
    /// @return TaskKind::MD or TaskKind::SinglePoint.
    TaskKind kind() const noexcept { return kind_; }

  private:
    TaskKind kind_;
    std::int64_t n_steps_;
    std::optional<md::MD> md_;
    std::optional<Evaluate> eval_;
};

} // namespace task
} // namespace cpu
} // namespace ATV
