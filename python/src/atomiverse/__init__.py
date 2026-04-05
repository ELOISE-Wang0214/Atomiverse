"""Atomiverse: molecular simulation compute backend.

Provides CPU (and optionally CUDA) simulation drivers that operate on
molpy Frame objects via zero-copy numpy transfer.
"""

try:
    from atomiverse._atomiverse import (  # installed wheel
        CpuFrame, MD, Phase, PhaseBuilder,
        VelocityVerlet, InitVelocity, NoseHoover,
        LJ, NeighborList, HarmonicBond, HarmonicAngle, MultiHarmonicDihedral,
        Thermo,
        # SCF
        IntegralPack, RHFResult, UHFResult, UKSResult, XCType,
        build_integrals_sto3g, rhf_solve, uhf_solve, uks_solve,
    )
except ImportError:
    from _atomiverse import (  # development: .so in PYTHONPATH
        CpuFrame, MD, Phase, PhaseBuilder,
        VelocityVerlet, InitVelocity, NoseHoover,
        LJ, NeighborList, HarmonicBond, HarmonicAngle, MultiHarmonicDihedral,
        Thermo,
        # SCF
        IntegralPack, RHFResult, UHFResult, UKSResult, XCType,
        build_integrals_sto3g, rhf_solve, uhf_solve, uks_solve,
    )

__all__ = [
    # Frame
    "CpuFrame",
    # MD
    "MD", "Phase", "PhaseBuilder",
    "VelocityVerlet", "InitVelocity", "NoseHoover",
    "LJ", "NeighborList", "HarmonicBond", "HarmonicAngle", "MultiHarmonicDihedral",
    "Thermo",
    # SCF
    "IntegralPack", "RHFResult", "UHFResult", "UKSResult", "XCType",
    "build_integrals_sto3g", "rhf_solve", "uhf_solve", "uks_solve",
    # Helpers
    "run_md",
]


def run_md(frame, sim, n_steps):
    """Run MD simulation on a molpy Frame.

    Transfers data from molpy Frame to cpu::Frame, runs the simulation,
    and writes results back.

    Parameters
    ----------
    frame : molpy.Frame
        Input frame with "atoms" block (x, y, z, mass, type required).
    sim : MD
        Configured and built MD driver.
    n_steps : int
        Number of simulation steps.

    Returns
    -------
    molpy.Frame
        The same frame object, updated with simulation results.
    """
    import numpy as np

    atoms = frame["atoms"]
    n_atoms = atoms.nrows

    cf = CpuFrame(n_atoms)

    # Populate atoms
    cf.set_positions(
        np.ascontiguousarray(atoms["x"], dtype=np.float32),
        np.ascontiguousarray(atoms["y"], dtype=np.float32),
        np.ascontiguousarray(atoms["z"], dtype=np.float32),
    )
    cf.set_mass(np.ascontiguousarray(atoms["mass"], dtype=np.float32))
    cf.set_type(np.ascontiguousarray(atoms["type"], dtype=np.float32))

    if "vx" in atoms:
        cf.set_velocities(
            np.ascontiguousarray(atoms["vx"], dtype=np.float32),
            np.ascontiguousarray(atoms["vy"], dtype=np.float32),
            np.ascontiguousarray(atoms["vz"], dtype=np.float32),
        )

    # Box
    if frame.box is not None:
        cf.set_box(np.ascontiguousarray(
            np.asarray(frame.box.matrix), dtype=np.float32
        ))

    # Topology
    if "bonds" in frame:
        bonds = frame["bonds"]
        cf.set_bonds(
            np.ascontiguousarray(bonds["atomi"], dtype=np.int32),
            np.ascontiguousarray(bonds["atomj"], dtype=np.int32),
            np.ascontiguousarray(bonds["type"], dtype=np.int32),
        )

    if "angles" in frame:
        angles = frame["angles"]
        cf.set_angles(
            np.ascontiguousarray(angles["atomi"], dtype=np.int32),
            np.ascontiguousarray(angles["atomj"], dtype=np.int32),
            np.ascontiguousarray(angles["atomk"], dtype=np.int32),
            np.ascontiguousarray(angles["type"], dtype=np.int32),
        )

    if "dihedrals" in frame:
        dihedrals = frame["dihedrals"]
        cf.set_dihedrals(
            np.ascontiguousarray(dihedrals["atomi"], dtype=np.int32),
            np.ascontiguousarray(dihedrals["atomj"], dtype=np.int32),
            np.ascontiguousarray(dihedrals["atomk"], dtype=np.int32),
            np.ascontiguousarray(dihedrals["atoml"], dtype=np.int32),
            np.ascontiguousarray(dihedrals["type"], dtype=np.int32),
        )

    # Run
    sim.run(cf, n_steps)

    # Write results back
    atoms["x"] = cf.get_x()
    atoms["y"] = cf.get_y()
    atoms["z"] = cf.get_z()
    atoms["vx"] = cf.get_vx()
    atoms["vy"] = cf.get_vy()
    atoms["vz"] = cf.get_vz()
    atoms["fx"] = cf.get_fx()
    atoms["fy"] = cf.get_fy()
    atoms["fz"] = cf.get_fz()

    frame.metadata["energy"] = float(cf.energy)
    frame.metadata["kinetic"] = float(cf.kinetic)

    return frame
