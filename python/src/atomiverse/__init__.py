"""Atomiverse: molecular simulation compute backend for molpy.

Users work with ``mp.Frame`` — Atomiverse handles conversion internally.

Example — MD::

    import molpy as mp
    import atomiverse as av

    frame = mp.Frame({"atoms": mp.Block({
        "x": [0, 3.4], "y": [0, 0], "z": [0, 0],
        "type": [0, 0], "mass": [39.948, 39.948],
    })})
    frame.box = mp.Box(np.eye(3) * 30.0)

    sim = av.MD(dt=0.001)
    sim.add(av.InitVelocity(temperature=10.0, seed=42))
    sim.on("force").add(av.NeighborList(cutoff=8.5), order=-10)
    sim.on("force").add(av.LJ(epsilon=0.0104, sigma=3.4, cutoff=8.5))
    sim.on("end_of_step").add(av.Thermo(freq=100))
    sim.run(frame, 1000)

Example — RHF::

    frame = mp.Frame({"atoms": mp.Block({
        "x": [0, 0.74], "y": [0, 0], "z": [0, 0],
        "type": [1, 1], "mass": [1.008, 1.008],
    })})
    frame.box = mp.Box(np.eye(3) * 100.0)

    rhf = av.RHF(frame, basis="sto-3g")
    result = rhf.solve(tol=1e-10)
    print(result.energy)

Example — Unified SCF::

    scf = av.SCF(frame, method="uks", xc="PBE")
    result = scf.solve()
    print(result.xc_energy)

Example — Geometry Optimization::

    opt = av.GeometryOptimizer(frame, method="rhf")
    result = opt.run()
    print(result.final_energy)
"""

# ── Drivers (OOP wrappers) ──
from atomiverse.md import MD
from atomiverse.scf import (
    SCF, SCFResult, GradientResult, Integrals,
    RHF, RHFResult,
    UHF, UHFResult,
    UKS, UKSResult,
)
from atomiverse.optimizer import GeometryOptimizer, OptResult

# ── Components (value objects — direct from C++) ──
from atomiverse._binding import _cpp as _cpp

VelocityVerlet = _cpp.VelocityVerlet
InitVelocity = _cpp.InitVelocity
NoseHoover = _cpp.NoseHoover
LJ = _cpp.LJ
NeighborList = _cpp.NeighborList
HarmonicBond = _cpp.HarmonicBond
HarmonicAngle = _cpp.HarmonicAngle
MultiHarmonicDihedral = _cpp.MultiHarmonicDihedral
Thermo = _cpp.Thermo


def register_molpy() -> None:
    """Register Atomiverse device handlers into molpy's device registry."""
    from atomiverse._device import register_devices
    register_devices()


__all__ = [
    # Drivers
    "MD",
    # SCF (unified)
    "SCF", "SCFResult", "GradientResult", "Integrals",
    # SCF (convenience wrappers, backward compatible)
    "RHF", "RHFResult",
    "UHF", "UHFResult",
    "UKS", "UKSResult",
    # Optimizer
    "GeometryOptimizer", "OptResult",
    # Components
    "VelocityVerlet", "InitVelocity", "NoseHoover",
    "LJ", "NeighborList", "HarmonicBond", "HarmonicAngle",
    "MultiHarmonicDihedral",
    "Thermo",
    # molpy integration
    "register_molpy",
]
