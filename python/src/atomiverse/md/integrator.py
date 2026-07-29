"""Time integrators (mirrors ``cpu/dynamics/`` integrator components).

Canonical deep path: ``atomiverse.md.integrator.VelocityVerlet``; the
top-level convenience name ``atomiverse.VelocityVerlet`` is the same class.

Install via ``MD.set_integrator(VelocityVerlet(dt))``.
"""

from atomiverse._binding import _cpp

#: Symplectic velocity-Verlet integrator — ``cpu::md::VelocityVerlet(dt [fs])``.
VelocityVerlet = _cpp.VelocityVerlet

__all__ = ["VelocityVerlet"]
