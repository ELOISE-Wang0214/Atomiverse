/// @file cpu/md/md.h
/// @brief cpu::md::MD — concrete time-advance driver under DynamicsDriver.
#pragma once

#include "atomiverse/core/device.h"
#include "atomiverse/cpu/backend_traits.h"
#include "atomiverse/cpu/dynamics_driver.h"
#include "atomiverse/cpu/snapshot_publisher.h"
#include "atomiverse/cpu/md/init.h"
#include "atomiverse/cpu/md/mpi/atom_migration.h"
#include "atomiverse/cpu/md/mpi/halo_exchange.h"
#include "atomiverse/cpu/md/mpi/reduce.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cpu {
namespace md {

/// @brief CPU MD simulation driver.
///
/// Convenience class that auto-adds infrastructure components.
///
/// End-of-step ordering: `ZeroKinetic(-1000) → KineticEnergy(-100) →
/// ReduceEnergy(-50) → SnapshotPublisher(+10) → user observers (default 0+)`.
/// pre_force ordering includes `AtomMigration(-30)`, which follows the
/// LAMMPS reneighbor cadence by attempting exchange every 10 steps.
/// Multi-rank pair simulations must explicitly add `HaloExchangeForward{cutoff}`
/// at `pre_force @ -20` and `HaloExchangeReverse{}` at `post_force @ -1000`.
/// (Negative-order infrastructure runs before user-positive components in each phase.)
/// In single-rank mode every MPI stub is a no-op and `SnapshotPublisher`
/// takes its 5 ns zero-subscriber early return path.
///
/// Usage:
/// @code{.cpp}
/// cpu::md::MD sim;
/// sim.add(cpu::md::VelocityVerlet{0.001});
/// sim.on(ATV::md::Phase::force).add(cpu::md::NeighborList{6.8}, -10);
/// sim.build();
/// sim.run(frame, 1000);
/// @endcode
class MD : public DynamicsDriver {
public:
    MD() {
        // Pre-existing zero/kinetic infrastructure.
        add(ZeroForces{});
        add(ZeroEnergy{});
        add(ZeroStress{});   // pair-virial accumulator must start at 0
        add(ZeroKinetic{});
        add(KineticEnergy{});

        // MPI infrastructure (single-rank no-ops in Phase 1).
        add(mpi::AtomMigration{10});
        add(mpi::ReduceEnergy{});

        // Async-analysis publisher -- shared infrastructure, explicitly
        // bound to *this driver's* end_of_step @ +10. Subscribers
        // inheriting `cpu::md::AsyncObserver` default to +20 on the same
        // phase and see this publisher's slot.
        this->on(ATV::md::Phase::end_of_step)
            .add(cpu::SnapshotPublisher{}, 10);
    }
};

} // namespace md
} // namespace cpu
} // namespace ATV
