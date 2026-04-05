"""Shared pytest fixtures for Atomiverse Python binding tests."""

import numpy as np
import pytest

import atomiverse as av


@pytest.fixture
def n_atoms() -> int:
    return 8


@pytest.fixture
def argon_frame(n_atoms: int) -> av.CpuFrame:
    """Single-type Ar frame on a simple cubic lattice, box 10 Å."""
    cf = av.CpuFrame(n_atoms)
    side = int(round(n_atoms ** (1 / 3)))
    positions = [
        (ix, iy, iz)
        for ix in range(side)
        for iy in range(side)
        for iz in range(side)
    ]
    scale = 3.4  # ~σ_Ar [Å]
    xs = np.array([p[0] * scale for p in positions], dtype=np.float32)
    ys = np.array([p[1] * scale for p in positions], dtype=np.float32)
    zs = np.array([p[2] * scale for p in positions], dtype=np.float32)
    cf.set_positions(xs, ys, zs)
    cf.set_mass(np.full(n_atoms, 39.948, dtype=np.float32))
    cf.set_type(np.zeros(n_atoms, dtype=np.int32))
    box = np.eye(3, dtype=np.float32) * (side * scale)
    cf.set_box(box)
    return cf


@pytest.fixture
def lj_md() -> av.MD:
    """Minimal LJ MD: VelocityVerlet + InitVelocity; no force components."""
    sim = av.MD()
    sim.add_velocity_verlet(av.VelocityVerlet(0.001))
    sim.add_init_velocity(av.InitVelocity(10.0, 42))
    return sim


@pytest.fixture
def full_lj_md() -> av.MD:
    """Full LJ MD with NeighborList + LJ potential."""
    cutoff = 8.5
    sim = av.MD()
    sim.add_velocity_verlet(av.VelocityVerlet(0.001))
    sim.add_init_velocity(av.InitVelocity(10.0, 42))
    sim.on(av.Phase.force).add_neighbor_list(av.NeighborList(cutoff), -10).add_lj(
        av.LJ(0.0104, 3.4, cutoff)
    )
    return sim
