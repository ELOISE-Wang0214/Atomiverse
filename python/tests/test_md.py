"""Tests for the OOP MD driver with mp.Frame."""

import numpy as np
import pytest
import molpy as mp

import atomiverse as av


class TestComponentConstruction:
    def test_velocity_verlet(self):
        assert av.VelocityVerlet(0.001) is not None

    def test_init_velocity(self):
        assert av.InitVelocity(300.0, 42) is not None

    def test_nose_hoover(self):
        assert av.NoseHoover(300.0, 100.0, 0.001) is not None

    def test_lj(self):
        assert av.LJ(0.0104, 3.4, 8.5) is not None

    def test_neighbor_list(self):
        assert av.NeighborList(8.5) is not None

    def test_neighbor_list_with_skin(self):
        assert av.NeighborList(8.5, 0.5) is not None

    def test_thermo(self):
        assert av.Thermo(100) is not None

    def test_harmonic_bond(self):
        assert av.HarmonicBond([100.0], [1.54]) is not None

    def test_harmonic_angle(self):
        assert av.HarmonicAngle([60.0], [1.91]) is not None

    def test_multi_harmonic_dihedral(self):
        assert av.MultiHarmonicDihedral([[1.0, -1.0, 0.0, 0.0, 0.0]]) is not None


class TestMDDriverBuild:
    def test_empty_build(self):
        sim = av.MD()
        sim.build()

    def test_dt_shorthand(self):
        sim = av.MD(dt=0.001)
        sim.build()

    def test_add_components(self):
        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(10.0, 42))
        sim.add(av.NoseHoover(300.0, 100.0, 0.001))
        sim.build()

    def test_on_returns_phase_builder(self):
        sim = av.MD()
        assert sim.on("force") is not None

    def test_phase_builder_chain(self):
        sim = av.MD(dt=0.001)
        sim.on("force").add(av.NeighborList(8.5), order=-10) \
                       .add(av.LJ(0.0104, 3.4, 8.5))
        sim.build()

    def test_invalid_phase_raises(self):
        sim = av.MD()
        with pytest.raises(ValueError, match="Unknown phase"):
            sim.on("nonexistent")

    def test_invalid_driver_component_raises(self):
        sim = av.MD()
        with pytest.raises(TypeError, match="Cannot add"):
            sim.add(av.LJ(0.01, 3.4, 8.5))

    def test_invalid_phase_component_raises(self):
        builder = av.MD().on("force")
        with pytest.raises(TypeError, match="Cannot add"):
            builder.add(av.VelocityVerlet(0.001))


class TestMDRun:
    def test_velocity_init_sets_velocities(self, argon_frame):
        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(10.0, 42))
        sim.run(argon_frame, 1)
        vx = np.asarray(argon_frame["atoms"]["vx"])
        assert not np.allclose(vx, 0.0), "Velocities should be initialized"

    def test_kinetic_energy_positive_after_init(self, argon_frame):
        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(10.0, 42))
        sim.run(argon_frame, 1)
        assert argon_frame.metadata["kinetic"] > 0.0

    def test_positions_change_after_md(self, argon_frame):
        x_before = np.asarray(argon_frame["atoms"]["x"]).copy()
        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(10.0, 42))
        sim.run(argon_frame, 10)
        x_after = np.asarray(argon_frame["atoms"]["x"])
        assert not np.allclose(x_before, x_after), "Positions should move"

    def test_lj_md_energy_finite(self, argon_frame, full_lj_md):
        full_lj_md.run(argon_frame, 5)
        assert np.isfinite(argon_frame.metadata["energy"])
        assert np.isfinite(argon_frame.metadata["kinetic"])

    def test_lj_md_forces_finite(self, argon_frame, full_lj_md):
        full_lj_md.run(argon_frame, 1)
        fx = np.asarray(argon_frame["atoms"]["fx"])
        assert np.all(np.isfinite(fx)), "Forces should be finite"

    def test_auto_build(self, argon_frame):
        """run() should call build() automatically."""
        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(10.0, 42))
        sim.run(argon_frame, 1)

    def test_missing_mass_raises(self):
        frame = mp.Frame({"atoms": mp.Block({
            "x": [0.0], "y": [0.0], "z": [0.0],
            "type": np.array([0], dtype=np.int32),
            # no mass!
        })})
        frame.box = mp.Box(np.eye(3) * 10.0)
        sim = av.MD(dt=0.001)
        with pytest.raises(ValueError, match="mass"):
            sim.run(frame, 1)

    def test_summary_returns_string(self, argon_frame):
        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(10.0, 42))
        sim.enable_timing()
        sim.run(argon_frame, 1)
        assert isinstance(sim.summary(), str)

    def test_run_zero_steps(self, argon_frame):
        sim = av.MD(dt=0.001)
        sim.run(argon_frame, 0)

    def test_repr(self):
        sim = av.MD(dt=0.001)
        assert "not built" in repr(sim)
        sim.build()
        assert "built" in repr(sim)
