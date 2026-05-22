"""Atomiverse: molecular simulation compute backend for molrs.

Atomiverse binds drivers and components only — it does not bind a Frame
data type. Simulation state crosses the boundary as a ``molrs.Frame``;
each driver ``run()`` borrows it and writes results back in place.

Example — Molecular Dynamics::

    import numpy as np
    import molrs
    import atomiverse as av

    frame = molrs.Frame()
    atoms = molrs.Block()
    atoms.insert("x", np.zeros(8)); atoms.insert("y", np.zeros(8))
    atoms.insert("z", np.zeros(8))
    atoms.insert("mass", np.full(8, 39.948))
    atoms.insert("type", np.zeros(8, dtype=np.int32))
    frame["atoms"] = atoms
    frame.box = molrs.Box.cube(20.0)

    sim = av.MD(dt=0.001)
    sim.add(av.InitVelocity(10.0, 42))
    sim.run(frame, 100)

Example — Plane-Wave DFT::

    result = av.pw_scf(frame, "Si.upf", ecut_ry=30.0, k_mesh=4)
    print(result.energy)
"""

from atomiverse.md import (
    MD,
    HarmonicAngle,
    HarmonicBond,
    InitVelocity,
    LJ,
    MultiHarmonicDihedral,
    NeighborList,
    NoseHoover,
    Thermo,
    VelocityVerlet,
)
from atomiverse.pw import PWResult, pw_scf

__all__ = [
    "MD",
    "VelocityVerlet",
    "InitVelocity",
    "NoseHoover",
    "LJ",
    "NeighborList",
    "HarmonicBond",
    "HarmonicAngle",
    "MultiHarmonicDihedral",
    "Thermo",
    "pw_scf",
    "PWResult",
]
