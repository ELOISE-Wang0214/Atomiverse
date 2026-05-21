/// @file cuda/md/md.cuh
/// @brief md::MD -- molecular dynamics driver with LAMMPS-style phases.
///
/// Automatically adds infrastructure components:
/// - ZeroForces at pre_force (order -1000)
/// - ZeroEnergy at pre_force (order -999)
///
/// Usage:
/// @code{.cpp}
/// md::MD sim;
/// sim.add(md::VelocityVerlet{0.001});
/// sim.on(md::Phase::force)
///    .add(md::NeighborList{3.0}, -10)
///    .add(md::LJ{1.0, 1.0, 3.0});
/// sim.build();
/// sim.run(frame, 10000);
/// @endcode
///
/// @see md::Phase for the LAMMPS-style phase sequence
/// @see Driver for the inherited fluent API
#pragma once

#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/md/kinetic_energy.cuh"
#include "atomiverse/facade/driver.h"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace md {

/// @brief MD driver with LAMMPS-style phase sequence.
///
/// Automatically adds infrastructure components:
/// - ZeroForces at pre_force (order -1000) -- zeros fx, fy, fz.
/// - ZeroEnergy at pre_force (order -999) -- zeros potential energy.
///
/// @see minimizer::CG for the energy minimizer driver
class MD : public Driver<Phase> {
  public:
    /// @brief Construct MD driver with auto-added zeroing infrastructure.
    MD() {
        on(Phase::pre_force).add(ZeroForces{}, -1000);
        on(Phase::pre_force).add(ZeroEnergy{}, -999);
        on(Phase::end_of_step).add(ZeroKinetic{}, -1000);
        on(Phase::end_of_step).add(KineticEnergy{}, -100);
    }
};

} // namespace md
} // namespace ATV
