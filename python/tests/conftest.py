"""Shared pytest fixtures for Atomiverse Python binding tests.

Atomiverse drivers operate directly on a ``molrs.Frame`` — there is no
separate Atomiverse Frame type. Fixtures build ``molrs.Frame`` instances.
"""

import numpy as np
import pytest
import molrs

import atomiverse as av


def _argon_columns(n_atoms: int):
    """Simple-cubic Ar lattice columns for ``n_atoms`` atoms."""
    side = int(round(n_atoms ** (1 / 3)))
    scale = 3.4  # ~sigma_Ar [A]
    xs, ys, zs = [], [], []
    for ix in range(side):
        for iy in range(side):
            for iz in range(side):
                xs.append(ix * scale)
                ys.append(iy * scale)
                zs.append(iz * scale)
    box_len = side * scale
    return (
        np.array(xs, dtype=np.float64),
        np.array(ys, dtype=np.float64),
        np.array(zs, dtype=np.float64),
        box_len,
    )


@pytest.fixture
def n_atoms() -> int:
    return 8


@pytest.fixture
def argon_frame(n_atoms: int) -> "molrs.Frame":
    """Single-type Ar molrs.Frame on a simple cubic lattice."""
    xs, ys, zs, box_len = _argon_columns(n_atoms)
    frame = molrs.Frame()
    atoms = molrs.Block()
    atoms.insert("x", xs)
    atoms.insert("y", ys)
    atoms.insert("z", zs)
    atoms.insert("mass", np.full(n_atoms, 39.948, dtype=np.float64))
    atoms.insert("type", np.zeros(n_atoms, dtype=np.int32))
    frame["atoms"] = atoms
    frame.box = molrs.Box.cube(box_len)
    return frame


@pytest.fixture
def full_lj_md() -> av.MD:
    """Full LJ MD: VelocityVerlet + InitVelocity + NeighborList + LJ."""
    cutoff = 8.5
    sim = av.MD(dt=0.001)
    sim.add(av.InitVelocity(10.0, 42))
    sim.on("force").add(av.NeighborList(cutoff), order=-10) \
                   .add(av.LJ(0.0104, 3.4, cutoff))
    return sim
