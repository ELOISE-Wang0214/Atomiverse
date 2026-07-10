"""Atomiverse: molecular simulation compute backend for molrs.

Atomiverse binds drivers and components only — it does not bind a Frame
data type. Simulation state crosses the boundary as a ``molrs.Frame``;
each driver ``run()`` borrows it and writes results back in place.

Current surface (honest capability boundary)
--------------------------------------------
- **MD** — the ``MD`` driver plus its components, exported under the
  canonical C++ names (LAMMPS-style): ``VelocityVerlet``, ``InitVelocity``,
  ``NoseHoover``, ``PairLJCut``, ``NeighborList``, ``BondHarmonic``,
  ``AngleHarmonic``, ``DihedralMultiHarmonic``, ``Thermo``.
- **Plane-wave DFT** — ``pw_scf`` / ``PWResult`` (requires a wheel built
  with ``ATV_BUILD_PW=ON``; calling it on a PW-less wheel raises a
  ``RuntimeError`` with rebuild guidance).

Roadmap: Gaussian-SCF and geometry-optimization bindings (the C++
``SCFEvaluator`` family and the ``Minimize`` driver with ``LBFGS``/``CG``
components) are future work with their own spec — they are deliberately
NOT exposed here rather than exposed broken.

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
    sim.on("force").add(av.NeighborList(8.5), order=-10) \\
                   .add(av.PairLJCut(0.0104, 3.4, 8.5))
    sim.run(frame, 100)

Example — Plane-Wave DFT::

    result = av.pw_scf(frame, "Si.upf", ecut_ry=30.0, k_mesh=4)
    print(result.energy)
"""

from atomiverse.md import (
    MD,
    AngleHarmonic,
    BondHarmonic,
    DihedralMultiHarmonic,
    InitVelocity,
    NeighborList,
    NoseHoover,
    PairLJCut,
    Thermo,
    VelocityVerlet,
)
from atomiverse.pw import PWResult, pw_scf

__all__ = [
    "MD",
    "VelocityVerlet",
    "InitVelocity",
    "NoseHoover",
    "PairLJCut",
    "NeighborList",
    "BondHarmonic",
    "AngleHarmonic",
    "DihedralMultiHarmonic",
    "Thermo",
    "pw_scf",
    "PWResult",
]
