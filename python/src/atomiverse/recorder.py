"""I/O recorders (mirrors ``cpu/record/``).

Canonical deep paths: ``atomiverse.recorder.Thermo`` and
``atomiverse.recorder.XYZ``; the top-level ``atomiverse.Thermo`` /
``atomiverse.XYZ`` are the same classes.

Prefer ``MD.set_thermo(N)`` for the common LAMMPS-style cadence; construct
``Thermo(freq)`` only when passing through ``set_record`` or the phase
escape hatch. ``XYZ(filename, freq)`` always goes through ``set_record``.

Both recorders are bound with the MD driver, so they exist on an extension
built with the ``md`` module (``ATV_BUILD_MD=ON``); on one built without it,
importing this module raises ``MissingDomainError`` naming the flag.
"""

from atomiverse._binding import _cpp, require_domain

require_domain(__name__, "md")

#: Thermodynamic output at ``end_of_step`` — ``cpu::Thermo(freq [steps])``.
Thermo = _cpp.Thermo

#: XYZ trajectory output at ``end_of_step`` —
#: ``cpu::XyzRecorder(filename, freq [steps])``.
XYZ = _cpp.XYZ

__all__ = ["Thermo", "XYZ"]
