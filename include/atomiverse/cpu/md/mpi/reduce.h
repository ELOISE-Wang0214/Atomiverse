/// @file cpu/md/mpi/reduce.h
/// @brief End-of-step MPI reductions for MD: energy, kinetic, stress.
///
/// `ReduceEnergy` does one batched `MPI_Allreduce(SUM)` of the 18 per-rank
/// scalar partials `{energy, energy_pair, energy_bond, energy_angle,
/// energy_dihedral, kinetic, stress[6], stress_kinetic[6]}` into a single
/// collective, so per-rank bookkeeping is globally consistent before
/// observers or the snapshot publisher read the frame's scalar state.
/// `n_local` is NOT reduced — it is authoritative per rank from
/// `dd->n_local`. Single-rank: no-op early-return.
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

/// @brief Batched 18-double `MPI_Allreduce(SUM)` of the per-rank energy,
/// kinetic, and stress partials. Single-rank: no-op (state already global).
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
