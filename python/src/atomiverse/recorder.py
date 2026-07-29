"""I/O recorders (mirrors ``cpu/record/``).

Canonical deep path: ``atomiverse.recorder.Thermo``; top-level
``atomiverse.Thermo`` is the same class.

Prefer ``MD.set_thermo(N)`` for the common LAMMPS-style cadence; construct
``Thermo(freq)`` only when passing through ``set_record`` or the phase
escape hatch.
"""

from atomiverse._binding import _cpp

#: Thermodynamic output at ``end_of_step`` — ``cpu::Thermo(freq [steps])``.
Thermo = _cpp.Thermo

__all__ = ["Thermo"]
