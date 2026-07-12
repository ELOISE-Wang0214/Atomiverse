/// @file facade/component_role.h
/// @brief ComponentRole + Capability — orthogonal tagging axes for components.
///
/// Two independent axes describe a component:
///
///   1. **Role** (structural): what slot does this component fill in the
///      outer driver loop? (Integrator, Optimizer, Thermostat, Evaluator, ...)
///   2. **Capability** (data-producing): what fields can this component write
///      to the frame? (Energy, Force, Stress, DensityMatrix, Position, Velocity, ...)
///
/// **Capability is the write-set, and it means it.** The bitset originally carried
/// only the fields an *Evaluator* produces, so every component that moved atoms —
/// VelocityVerlet, NoseHoover, Deform, ZeroLinearMomentum, the whole integrator and
/// thermostat population — declared `Capability::None` while writing the two most
/// fundamental fields in the frame. The doc said write-set; the code meant
/// evaluator-output. `Position` and `Velocity` close that hole: a component now
/// declares what it actually touches, which is the only basis on which a future
/// engine could detect two components fighting over `v` in the same phase.
///
/// Splitting them is what lets a single `SCFEvaluator` class be valid as
/// energy-only (in `Evaluate`), energy+density (single-point), and
/// energy+force (in `MD` / `Minimize`) — by adjusting its
/// `capabilities()` at construction time.
///
/// Each **preset** validates the registered component set at `build()` time using these
/// tags — a run that cannot possibly do what was asked of it fails at build, not after
/// an hour of wall clock:
///
///   - `cpu::Evaluate` requires ≥ 1 component with role==Evaluator and
///     capabilities() != None.
///   - `cpu::md::MD` requires ≥ 1 Integrator AND ≥ 1 component advertising
///     `Capability::Force` (capability, not role: a Thermostat-role Langevin that writes
///     forces satisfies "the integrator has something to integrate").
///   - `cpu::minimizer::Minimize` requires ≥ 1 Optimizer AND ≥ 1 Evaluator that
///     advertises `Capability::Energy` or `Capability::Force`.
///   - `cpu::vib::Hessian` requires `Capability::Force` AND `Capability::Hessian`.
///   - `cpu::mc::MC` requires ≥ 1 component with role==Mover AND ≥ 1 component
///     advertising `Capability::Energy` (the trial-energy pipeline).
///   - `cpu::PathDriver` requires an Optimizer AND a Force evaluator.
///
/// These are *requirements*, not a taxonomy. There are no driver "category" base classes
/// any more: a driver is a staged loop, and termination is a `LoopHandler`.
///
/// @see cpu/driver/cpu_driver.h — the loop node the presets derive
/// @see cpu/driver/loop_handler.h — termination as a policy
#pragma once

#include <cstdint>

namespace ATV {

/// @brief Structural role of a component within an outer driver loop.
///
/// Roles are orthogonal to capabilities. A component at the
/// "Evaluator" position (Role::Evaluator) describes *what* it computes
/// through `capabilities()`, not through the role tag.
enum class ComponentRole : uint8_t {
    None =
        0, ///< Untagged / structural infrastructure (Zero*, KineticEnergy, recorders, observers — set explicitly via Recorder/Observer if needed).
    Evaluator, ///< Computes one or more fields on the frame; capabilities() says which.
    Integrator, ///< Time-step integrator (VelocityVerlet, ...).
    Optimizer,  ///< Optimization step (LBFGS, CG, FIRE, ...).
    Thermostat, ///< Temperature controller (NoseHoover, Langevin, ...).
    Barostat,   ///< Pressure controller (Berendsen, MTK, ...).
    Constraint, ///< Holds a geometric constraint (SETTLE, RATTLE): projects positions
                ///< and velocities back onto the constraint manifold after the
                ///< integrator has moved them. It is NOT an Integrator — tagging it
                ///< as one would let a Plan carrying a constraint but no actual
                ///< integrator pass DynamicsDriver's "at least one Integrator"
                ///< validation, i.e. run with nothing propagating the equations of
                ///< motion. It also REDUCES the system's degrees of freedom, which
                ///< the driver must publish before any thermostat caches Ndof.
    Mover, ///< Proposes a trial state change for a sampling driver (TranslateMove, ...).
    ConvergenceCheck, ///< Reports loop convergence to the outer driver via outer_status.
    Recorder, ///< Writes I/O (XYZ, Zarr, thermo).
    Observer, ///< Computes derived diagnostics (RDF, density, orbital info).
};

/// @brief Bitset of computable quantities a component can write.
///
/// A component declares the union of fields it may write each launch.
/// New bits are added as the project grows; reserved up to 32 bits.
/// When the bitspace is exhausted, switch to `std::bitset<N>`.
enum class Capability : uint32_t {
    None = 0,
    Energy = 1u << 0,
    Force = 1u << 1,
    Stress = 1u << 2,
    Charge = 1u << 3,
    Dipole = 1u << 4,
    DensityMatrix = 1u << 5,
    Orbitals = 1u << 6,
    Hessian =
        1u << 7, ///< Second-derivative matrix d2E/dRdR (Response category).
    Position = 1u << 8, ///< Atom coordinates x/y/z (integrators, constraints, Deform).
    Velocity =
        1u << 9, ///< Atom velocities vx/vy/vz (integrators, thermostats, constraints).
};

constexpr Capability operator|(Capability a, Capability b) noexcept {
    return static_cast<Capability>(static_cast<uint32_t>(a)
                                   | static_cast<uint32_t>(b));
}

constexpr Capability operator&(Capability a, Capability b) noexcept {
    return static_cast<Capability>(static_cast<uint32_t>(a)
                                   & static_cast<uint32_t>(b));
}

constexpr Capability& operator|=(Capability& a, Capability b) noexcept {
    a = a | b;
    return a;
}

/// @brief Test whether `set` carries every bit in `bits`.
constexpr bool has_capability(Capability set, Capability bits) noexcept {
    return (set & bits) == bits;
}

/// @brief True iff `c` carries no capability bits.
constexpr bool is_none(Capability c) noexcept { return c == Capability::None; }

} // namespace ATV
