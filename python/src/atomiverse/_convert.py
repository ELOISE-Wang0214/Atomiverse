"""mp.Frame <-> CpuFrame conversion.

Internal module. Users never call these directly — they are used inside
``MD.run()``, ``RHF()``, etc.

All failures are loud:
  - Missing blocks/fields → ValueError with exact field name
  - Wrong shapes → ValueError with expected vs actual
  - dtype mismatch → silently cast (this is fine, not a logic error)
"""

from __future__ import annotations

import numpy as np
import molpy as mp

from atomiverse._binding import _cpp

# ---------------------------------------------------------------------------
# Detect C++ Real / Int dtypes (runs once at import time)
# ---------------------------------------------------------------------------
_tmp = _cpp.CpuFrame(1)
REAL_DTYPE: np.dtype = np.asarray(_tmp.get_x()).dtype
INT_DTYPE: np.dtype = np.asarray(_tmp.get_type()).dtype
del _tmp


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _real(arr: np.ndarray) -> np.ndarray:
    """Contiguous Real array."""
    return np.ascontiguousarray(np.asarray(arr, dtype=REAL_DTYPE))


def _int(arr: np.ndarray) -> np.ndarray:
    """Contiguous Int array."""
    return np.ascontiguousarray(np.asarray(arr, dtype=INT_DTYPE))


def _require_block(frame: mp.Frame, name: str) -> mp.Block:
    if name not in frame:
        raise ValueError(f"mp.Frame must contain a '{name}' block")
    return frame[name]


def _require_fields(block: mp.Block, block_name: str, *fields: str) -> None:
    missing = [f for f in fields if f not in block]
    if missing:
        raise ValueError(
            f"'{block_name}' block missing required fields: {missing}"
        )


# ---------------------------------------------------------------------------
# mp.Frame → CpuFrame
# ---------------------------------------------------------------------------

def mp_to_cpu(frame: mp.Frame) -> _cpp.CpuFrame:
    """Convert an ``mp.Frame`` to an internal ``CpuFrame``.

    Required: ``frame["atoms"]`` with fields ``x, y, z, type``.

    Raises
    ------
    ValueError
        If required blocks or fields are missing.
    """
    atoms = _require_block(frame, "atoms")
    _require_fields(atoms, "atoms", "x", "y", "z", "type")

    n = atoms.nrows
    cf = _cpp.CpuFrame(n)

    cf.set_positions(_real(atoms["x"]), _real(atoms["y"]), _real(atoms["z"]))
    cf.set_type(_int(atoms["type"]))

    if "mass" in atoms:
        cf.set_mass(_real(atoms["mass"]))

    if "vx" in atoms and "vy" in atoms and "vz" in atoms:
        cf.set_velocities(
            _real(atoms["vx"]), _real(atoms["vy"]), _real(atoms["vz"])
        )

    if frame.box is not None:
        mat = np.asarray(frame.box.matrix, dtype=REAL_DTYPE)
        cf.set_box(np.ascontiguousarray(mat))

    # Topology
    if "bonds" in frame:
        bonds = frame["bonds"]
        _require_fields(bonds, "bonds", "atomi", "atomj", "type")
        cf.set_bonds(
            _int(bonds["atomi"]), _int(bonds["atomj"]), _int(bonds["type"])
        )

    if "angles" in frame:
        angles = frame["angles"]
        _require_fields(angles, "angles", "atomi", "atomj", "atomk", "type")
        cf.set_angles(
            _int(angles["atomi"]),
            _int(angles["atomj"]),
            _int(angles["atomk"]),
            _int(angles["type"]),
        )

    if "dihedrals" in frame:
        dihs = frame["dihedrals"]
        _require_fields(
            dihs, "dihedrals", "atomi", "atomj", "atomk", "atoml", "type"
        )
        cf.set_dihedrals(
            _int(dihs["atomi"]),
            _int(dihs["atomj"]),
            _int(dihs["atomk"]),
            _int(dihs["atoml"]),
            _int(dihs["type"]),
        )

    return cf


# ---------------------------------------------------------------------------
# CpuFrame → mp.Frame (write results back, in-place)
# ---------------------------------------------------------------------------

def cpu_to_mp(cf: _cpp.CpuFrame, frame: mp.Frame) -> None:
    """Write simulation results from ``CpuFrame`` back into ``mp.Frame``.

    Updates positions, velocities, forces in the atoms block and
    energy/kinetic in metadata.  Operates **in-place** on *frame*.
    """
    atoms = frame["atoms"]

    atoms["x"] = np.asarray(cf.get_x())
    atoms["y"] = np.asarray(cf.get_y())
    atoms["z"] = np.asarray(cf.get_z())
    atoms["vx"] = np.asarray(cf.get_vx())
    atoms["vy"] = np.asarray(cf.get_vy())
    atoms["vz"] = np.asarray(cf.get_vz())
    atoms["fx"] = np.asarray(cf.get_fx())
    atoms["fy"] = np.asarray(cf.get_fy())
    atoms["fz"] = np.asarray(cf.get_fz())

    frame.metadata["energy"] = float(cf.energy)
    frame.metadata["kinetic"] = float(cf.kinetic)
