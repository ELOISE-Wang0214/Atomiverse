"""Plane-wave DFT method API for periodic systems.

The PW-enabled wheel exposes the same object lifecycle as C++:

``PW`` (members) → ``build(frame)`` → ``run(frame)`` → ``PWResult``.

There is no bag named ``PWConfig`` / ``PWSolver`` and no transient functional
solve entry point.

Example::

    import atomiverse as av

    pw = av.PW()
    pw.upf_path = "Si.upf"
    pw.ecut_ry = 30.0
    pw.k_mesh = 4
    pw.build(frame)
    result = pw.run(frame)
    print(result.energy)  # Ry

When Atomiverse is built without ``ATV_BUILD_PW=ON``, this module remains
importable but exports no PW types. The package surface therefore reflects the
compiled capability instead of pretending an unavailable solver exists.
"""

from __future__ import annotations

from atomiverse._binding import _cpp

if hasattr(_cpp, "PW"):
    PW = _cpp.PW
    PWResult = _cpp.PWResult
    __all__ = ["PW", "PWResult"]
else:
    __all__: list[str] = []
