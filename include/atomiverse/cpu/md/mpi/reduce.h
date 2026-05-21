/// @file cpu/md/mpi/reduce.h
/// @brief End-of-step MPI reductions for MD: energy, kinetic, virial.
///
/// Phase 1: no-op stub with default binding at `end_of_step @ -50`.
/// Phase 3: batched `MPI_Allreduce` of `{energy, kinetic, virial[6], n_local}`
/// (9 doubles + 1 int) into a single collective so per-rank bookkeeping is
/// consistent before observers or the snapshot publisher read the frame's
/// scalar state.
///
/// Order is tight: must run after `KineticEnergy @ -100` so the reduced KE
/// is the locally-summed one, and before `SnapshotPublisher @ +10` so async
/// observers see globally-reduced scalars.
///
/// @see .claude/specs/parallel-md-engine.md  § A.4
#pragma once

#include <string>
#include <vector>

#include "atomiverse/core/binding.h"
#include "atomiverse/cpu/component.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

/// @brief Batched allreduce of `{energy, kinetic, virial[6], n_local}`.
///
/// Phase 1: structural no-op; single-rank state is already globally consistent.
class ReduceEnergy : public cpu::ComponentBase {
  public:
    std::string name() const override { return "reduce_energy"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::end_of_step), -50}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
