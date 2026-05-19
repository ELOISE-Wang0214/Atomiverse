"""Atomiverse: molecular simulation compute backend for molpy.

Example — Plane-Wave DFT::

    import numpy as np
    import atomiverse as av

    frame = av.Frame(2)
    frame.atoms.type[:] = [14, 14]
    frame.atoms.x[:] = [0.0, 1.3575]
    frame.atoms.y[:] = [0.0, 1.3575]
    frame.atoms.z[:] = [0.0, 1.3575]
    frame.box = ...  # FCC cell
    result = av.pw_scf(frame, "Si.upf", ecut_ry=30.0, k_mesh=4)
    print(result.energy)
"""

from atomiverse._binding import _cpp as _cpp
from atomiverse.pw import pw_scf, PWResult

Frame = _cpp.CpuFrame

__all__ = ["Frame", "pw_scf", "PWResult"]
