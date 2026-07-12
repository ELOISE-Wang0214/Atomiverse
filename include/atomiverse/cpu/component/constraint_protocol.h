/// @file cpu/component/constraint_protocol.h
/// @brief ConstraintProvider — opt-in protocol for components that remove degrees
///        of freedom. CUDA-free header.
///
/// Follows the `NeighborConsumer` / `HaloParticipant` precedent: `ComponentBase`
/// stays minimal, and a component that participates inherits the protocol and
/// overrides the matching `as_*` accessor. The driver probes the typed accessor —
/// no `dynamic_cast`, no `name()` matching.
#pragma once

namespace ATV {
namespace cpu {

class Frame;

/// @brief A component that constrains the system, and therefore removes degrees
///        of freedom from it.
///
/// **Why the engine has to know.** Temperature is `T = 2 KE / (Ndof kB)`. A rigid
/// water molecule has 6 degrees of freedom, not 9 — SETTLE removes 3 (two O-H bond
/// lengths and the H-O-H angle). If the engine keeps assuming `3N - 3`, a rigid-water
/// run reports a temperature ~33 % too low, AND every thermostat drives the system to
/// the wrong setpoint, because the thermostat's target kinetic energy is
/// `Ndof kB T / 2`. The bug is silent: nothing crashes, the trajectory just
/// equilibrates to the wrong temperature.
///
/// The count is per-frame (it scales with the number of constrained molecules), so
/// the driver asks for it after the frame is bound and before any component builds —
/// thermostats cache `Ndof` in their own `build()`.
class ConstraintProvider {
  public:
    virtual ~ConstraintProvider() = default;

    /// @brief Number of scalar constraints this component imposes.
    /// @param frame The bound system state (the count scales with it).
    /// @return Count of independent scalar constraints (e.g. 3 per rigid water).
    virtual int n_constraints(const cpu::Frame& frame) const = 0;
};

} // namespace cpu
} // namespace ATV
