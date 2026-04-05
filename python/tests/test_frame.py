"""Tests for CpuFrame numpy <-> cpu::Frame bridge."""

import numpy as np
import pytest

import atomiverse as av


class TestCpuFrameConstruction:
    def test_creates_with_n_atoms(self):
        cf = av.CpuFrame(100)
        assert cf.n_atoms == 100

    def test_zero_atoms_allowed(self):
        cf = av.CpuFrame(0)
        assert cf.n_atoms == 0

    def test_initial_energy_zero(self):
        cf = av.CpuFrame(4)
        assert cf.energy == pytest.approx(0.0)
        assert cf.kinetic == pytest.approx(0.0)


class TestPositions:
    def test_set_and_get_roundtrip(self):
        cf = av.CpuFrame(4)
        x = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
        y = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32)
        z = np.zeros(4, dtype=np.float32)
        cf.set_positions(x, y, z)
        np.testing.assert_allclose(cf.get_x(), x, rtol=1e-6)
        np.testing.assert_allclose(cf.get_y(), y, rtol=1e-6)
        np.testing.assert_allclose(cf.get_z(), z, rtol=1e-6)

    def test_positions_are_copies(self):
        cf = av.CpuFrame(2)
        x = np.array([1.0, 2.0], dtype=np.float32)
        cf.set_positions(x, np.zeros(2, dtype=np.float32), np.zeros(2, dtype=np.float32))
        x[0] = 99.0  # mutate original — should NOT affect frame
        assert cf.get_x()[0] == pytest.approx(1.0)

    def test_get_returns_numpy_array(self):
        cf = av.CpuFrame(3)
        cf.set_positions(
            np.zeros(3, dtype=np.float32),
            np.zeros(3, dtype=np.float32),
            np.zeros(3, dtype=np.float32),
        )
        assert isinstance(cf.get_x(), np.ndarray)
        assert cf.get_x().dtype == np.float32


class TestMassAndType:
    def test_set_and_get_mass(self):
        cf = av.CpuFrame(3)
        mass = np.array([39.948, 1.008, 12.011], dtype=np.float32)
        cf.set_mass(mass)
        np.testing.assert_allclose(cf.get_mass(), mass, rtol=1e-5)

    def test_set_and_get_type(self):
        cf = av.CpuFrame(3)
        types = np.array([0, 1, 2], dtype=np.int32)
        cf.set_type(types)
        np.testing.assert_array_equal(cf.get_type(), types)


class TestBox:
    def test_set_and_get_cubic_box(self):
        cf = av.CpuFrame(1)
        box = np.eye(3, dtype=np.float32) * 10.0
        cf.set_box(box)
        result = cf.get_box()
        np.testing.assert_allclose(result, box, rtol=1e-6)
        assert result.shape == (3, 3)

    def test_set_and_get_triclinic_box(self):
        cf = av.CpuFrame(1)
        box = np.array(
            [[10.0, 0.0, 0.0], [1.0, 9.8, 0.0], [0.5, 0.3, 9.6]], dtype=np.float32
        )
        cf.set_box(box)
        np.testing.assert_allclose(cf.get_box(), box, rtol=1e-5)


class TestVelocitiesAndForces:
    def test_set_and_get_velocities(self):
        cf = av.CpuFrame(2)
        vx = np.array([0.5, -0.5], dtype=np.float32)
        vy = np.zeros(2, dtype=np.float32)
        vz = np.zeros(2, dtype=np.float32)
        cf.set_velocities(vx, vy, vz)
        np.testing.assert_allclose(cf.get_vx(), vx, rtol=1e-6)

    def test_initial_forces_zero(self):
        cf = av.CpuFrame(4)
        # Forces are uninitialized until an MD step runs.
        # Just check shape/type.
        fx = cf.get_fx()
        assert fx.shape == (4,)
        assert fx.dtype == np.float32


class TestTopology:
    def test_set_bonds(self):
        cf = av.CpuFrame(4)
        atomi = np.array([0, 1, 2], dtype=np.int32)
        atomj = np.array([1, 2, 3], dtype=np.int32)
        btype = np.zeros(3, dtype=np.int32)
        cf.set_bonds(atomi, atomj, btype)  # should not raise

    def test_set_angles(self):
        cf = av.CpuFrame(4)
        cf.set_angles(
            np.array([0], dtype=np.int32),
            np.array([1], dtype=np.int32),
            np.array([2], dtype=np.int32),
            np.zeros(1, dtype=np.int32),
        )

    def test_set_dihedrals(self):
        cf = av.CpuFrame(4)
        cf.set_dihedrals(
            np.array([0], dtype=np.int32),
            np.array([1], dtype=np.int32),
            np.array([2], dtype=np.int32),
            np.array([3], dtype=np.int32),
            np.zeros(1, dtype=np.int32),
        )
