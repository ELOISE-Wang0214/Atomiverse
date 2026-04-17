/// @file cpu/md/mpi/atom_migration.h
/// @brief Atom migration component: moves owned atoms across rank boundaries.
///
/// Phase 1: structural no-op with a default binding at `pre_force @ -30`.
/// Phase 3: `MPI_Alltoallv` of atoms + bonded tuples to new owners whenever
/// this rank's ownership set has drifted beyond the halo skin; bumps
/// `DomainDecomposition::neighbor_generation` on each migration so
/// downstream consumers (neighbor list, halo send lists, bonded index
/// remap) rebuild from tag→local.
///
/// `launch()` never allocates. Capacity growth is handled out of band via
/// the driver-level `rebuild_frame()` hook — see spec § A.4.
///
/// @see .claude/specs/parallel-md-engine.md  § A.4, § A.5
#pragma once

#include <string>
#include <vector>

#include "atomiverse/cpu/component.h"
#include "atomiverse/core/binding.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {
namespace mpi {

/// @brief Atom migration across rank boundaries.
///
/// Phase 1: structural no-op. Registers `pre_force @ -30` so the component
/// graph has the correct slot reserved from the start.
class AtomMigration : public cpu::ComponentBase {
public:
    std::string name() const override { return "atom_migration"; }
    std::vector<DefaultBinding> default_bindings() const override {
        return {{static_cast<int>(ATV::md::Phase::pre_force), -30}};
    }
    void launch(const cpu::LaunchContext& ctx) override;
};

} // namespace mpi
} // namespace md
} // namespace cpu
} // namespace ATV
