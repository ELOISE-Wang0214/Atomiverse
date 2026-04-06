"""Shared pytest fixtures for Atomiverse Python binding tests."""

import numpy as np
import pytest
import molpy as mp

import atomiverse as av


@pytest.fixture
def n_atoms() -> int:
    return 8


@pytest.fixture
def argon_frame(n_atoms: int) -> mp.Frame:
    """Single-type Ar frame on a simple cubic lattice, box ~6.8 A."""
    side = int(round(n_atoms ** (1 / 3)))
    scale = 3.4  # ~sigma_Ar [A]

    xs, ys, zs = [], [], []
    for ix in range(side):
        for iy in range(side):
            for iz in range(side):
                xs.append(ix * scale)
                ys.append(iy * scale)
                zs.append(iz * scale)

    frame = mp.Frame({"atoms": mp.Block({
        "x": np.array(xs, dtype=np.float32),
        "y": np.array(ys, dtype=np.float32),
        "z": np.array(zs, dtype=np.float32),
        "mass": np.full(n_atoms, 39.948, dtype=np.float32),
        "type": np.zeros(n_atoms, dtype=np.int32),
    })})
    frame.box = mp.Box(np.eye(3) * (side * scale))
    return frame


@pytest.fixture
def lj_md() -> av.MD:
    """Minimal LJ MD: VelocityVerlet + InitVelocity; no force components."""
    sim = av.MD(dt=0.001)
    sim.add(av.InitVelocity(10.0, 42))
    return sim


@pytest.fixture
def full_lj_md() -> av.MD:
    """Full LJ MD with NeighborList + LJ potential."""
    cutoff = 8.5
    sim = av.MD(dt=0.001)
    sim.add(av.InitVelocity(10.0, 42))
    sim.on("force").add(av.NeighborList(cutoff), order=-10) \
                   .add(av.LJ(0.0104, 3.4, cutoff))
    return sim
