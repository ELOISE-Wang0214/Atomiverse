"""Tests for mp.Frame <-> CpuFrame conversion layer."""

import numpy as np
import pytest
import molpy as mp

from atomiverse._convert import mp_to_cpu, cpu_to_mp


class TestMpToCpu:
    def test_basic_roundtrip(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": np.array([1.0, 2.0], dtype=np.float32),
            "y": np.array([0.1, 0.2], dtype=np.float32),
            "z": np.zeros(2, dtype=np.float32),
            "type": np.array([0, 1], dtype=np.int32),
        })})

        cf = mp_to_cpu(frame)
        assert cf.n_atoms == 2
        np.testing.assert_allclose(cf.get_x(), [1.0, 2.0], rtol=1e-6)
        np.testing.assert_allclose(cf.get_y(), [0.1, 0.2], rtol=1e-6)

    def test_mass_transferred(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": np.zeros(3), "y": np.zeros(3), "z": np.zeros(3),
            "type": np.array([0, 1, 2], dtype=np.int32),
            "mass": np.array([39.948, 1.008, 12.011], dtype=np.float32),
        })})

        cf = mp_to_cpu(frame)
        np.testing.assert_allclose(cf.get_mass(), [39.948, 1.008, 12.011],
                                   rtol=1e-5)

    def test_box_transferred(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": [0.0], "y": [0.0], "z": [0.0],
            "type": np.array([0], dtype=np.int32),
        })})
        frame.box = mp.Box(np.eye(3) * 10.0)

        cf = mp_to_cpu(frame)
        np.testing.assert_allclose(cf.get_box(), np.eye(3) * 10.0, rtol=1e-6)

    def test_velocities_transferred(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": [0.0], "y": [0.0], "z": [0.0],
            "type": np.array([0], dtype=np.int32),
            "vx": np.array([0.5], dtype=np.float32),
            "vy": np.array([-0.5], dtype=np.float32),
            "vz": np.zeros(1, dtype=np.float32),
        })})

        cf = mp_to_cpu(frame)
        np.testing.assert_allclose(cf.get_vx(), [0.5], rtol=1e-6)

    def test_bonds_transferred(self):
        frame = mp.Frame({
            "atoms": mp.Block({
                "x": np.zeros(4), "y": np.zeros(4), "z": np.zeros(4),
                "type": np.zeros(4, dtype=np.int32),
            }),
            "bonds": mp.Block({
                "atomi": np.array([0, 1], dtype=np.int32),
                "atomj": np.array([1, 2], dtype=np.int32),
                "type": np.zeros(2, dtype=np.int32),
            }),
        })
        cf = mp_to_cpu(frame)  # should not raise

    def test_missing_atoms_raises(self):
        frame = mp.Frame()
        with pytest.raises(ValueError, match="atoms"):
            mp_to_cpu(frame)

    def test_missing_fields_raises(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": [0.0], "y": [0.0],
            # missing z and type
        })})
        with pytest.raises(ValueError, match="missing required fields"):
            mp_to_cpu(frame)


class TestCpuToMp:
    def test_writeback_updates_frame(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": np.array([0.0, 1.0], dtype=np.float32),
            "y": np.zeros(2, dtype=np.float32),
            "z": np.zeros(2, dtype=np.float32),
            "type": np.zeros(2, dtype=np.int32),
        })})

        cf = mp_to_cpu(frame)
        # Simulate: modify positions in CpuFrame
        cf.set_positions(
            np.array([10.0, 20.0], dtype=np.float32),
            np.zeros(2, dtype=np.float32),
            np.zeros(2, dtype=np.float32),
        )
        cpu_to_mp(cf, frame)

        np.testing.assert_allclose(frame["atoms"]["x"], [10.0, 20.0],
                                   rtol=1e-6)

    def test_writeback_sets_metadata(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": [0.0], "y": [0.0], "z": [0.0],
            "type": np.array([0], dtype=np.int32),
        })})

        cf = mp_to_cpu(frame)
        cpu_to_mp(cf, frame)

        assert "energy" in frame.metadata
        assert "kinetic" in frame.metadata
