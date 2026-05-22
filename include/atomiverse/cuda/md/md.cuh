/// @file cuda/md/md.cuh
/// @brief cuda::md::MD -- molecular dynamics driver with LAMMPS-style phases.
///
/// `MD` is the concrete CUDA molecular-dynamics driver. It is a
/// `cuda::DynamicsDriver` — the time-advance outer category — and so
/// inherits that base's continuable `run(frame, n_steps)` (full step
/// budget, honours `outer_status.stop_requested`, ignores
/// `outer_status.converged`) and its build-time `validate_components()`
/// contract (≥ 1 `Integrator`-role component and ≥ 1 component carrying
/// `Capability::Force`).
///
/// The constructor adds only infrastructure components — it contributes
/// no runtime behaviour beyond the base. The Zero* initializers are
/// phase-agnostic (`cuda/md/init.cuh`); `MD` binds them explicitly to
/// `md::Phase` via `on(phase).add()`:
/// - ZeroForces at pre_force (order -1000) -- zeros fx, fy, fz.
/// - ZeroEnergy at pre_force (order -999) -- zeros potential energy.
/// - ZeroKinetic at end_of_step (order -1000) -- zeros kinetic energy.
/// - KineticEnergy at end_of_step (order -100) -- accumulates KE.
///
/// Usage:
/// @code{.cpp}
/// cuda::md::MD sim;
/// sim.add(cuda::md::VelocityVerlet{0.001});
/// sim.on(md::Phase::force)
///    .add(cuda::NeighborList{3.0}, -10)
///    .add(cuda::LJ{1.0, 1.0, 3.0});
/// sim.build();          // validate_components() then compile the Plan
/// sim.run(frame, 10000);
/// @endcode
///
/// @see cuda::DynamicsDriver for the inherited loop + validation
/// @see md::Phase for the LAMMPS-style phase sequence
#pragma once

#include "atomiverse/cuda/dynamics_driver.cuh"
#include "atomiverse/cuda/md/init.cuh"
#include "atomiverse/cuda/md/kinetic_energy.cuh"
#include "atomiverse/facade/stage.h"

namespace ATV {
namespace cuda {
namespace md {

/// @brief CUDA MD driver — a `cuda::DynamicsDriver` with LAMMPS-style
///        phase sequence.
///
/// Inherits the time-advance loop and `validate_components()` from
/// `cuda::DynamicsDriver`; the constructor only binds the zeroing and
/// kinetic-energy infrastructure components.
///
/// @see cuda::minimizer::CG for the energy minimizer driver
class MD : public cuda::DynamicsDriver {
  public:
    /// @brief Construct MD driver, binding the phase-agnostic Zero* and
    ///        KineticEnergy infrastructure components to `md::Phase`.
    MD() {
        on(ATV::md::Phase::pre_force).add(ZeroForces{}, -1000);
        on(ATV::md::Phase::pre_force).add(ZeroEnergy{}, -999);
        on(ATV::md::Phase::end_of_step).add(ZeroKinetic{}, -1000);
        on(ATV::md::Phase::end_of_step).add(KineticEnergy{}, -100);
    }
};

} // namespace md
} // namespace cuda
} // namespace ATV
