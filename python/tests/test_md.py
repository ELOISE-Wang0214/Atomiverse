"""Tests for MD driver and component lifecycle."""

import numpy as np
import pytest

import atomiverse as av


class TestComponentConstruction:
    def test_velocity_verlet(self):
        vv = av.VelocityVerlet(0.001)
        assert vv is not None

    def test_init_velocity(self):
        iv = av.InitVelocity(300.0, 42)
        assert iv is not None

    def test_nose_hoover(self):
        nh = av.NoseHoover(300.0, 100.0, 0.001)
        assert nh is not None

    def test_lj(self):
        lj = av.LJ(0.0104, 3.4, 8.5)
        assert lj is not None

    def test_neighbor_list(self):
        nl = av.NeighborList(8.5)
        assert nl is not None

    def test_neighbor_list_with_skin(self):
        nl = av.NeighborList(8.5, 0.5)
        assert nl is not None

    def test_thermo(self):
        th = av.Thermo(100)
        assert th is not None

    def test_harmonic_bond(self):
        hb = av.HarmonicBond([100.0], [1.54])
        assert hb is not None

    def test_harmonic_angle(self):
        ha = av.HarmonicAngle([60.0], [1.91])  # ~109.5° in radians
        assert ha is not None

    def test_multi_harmonic_dihedral(self):
        mhd = av.MultiHarmonicDihedral([[1.0, -1.0, 0.0, 0.0, 0.0]])
        assert mhd is not None


class TestMDDriverBuild:
    def test_empty_build(self):
        sim = av.MD()
        sim.build()

    def test_add_velocity_verlet(self):
        sim = av.MD()
        sim.add_velocity_verlet(av.VelocityVerlet(0.001))
        sim.build()

    def test_on_phase_returns_builder(self):
        sim = av.MD()
        builder = sim.on(av.Phase.force)
        assert builder is not None

    def test_phase_builder_chain(self):
        sim = av.MD()
        builder = sim.on(av.Phase.force)
        builder.add_neighbor_list(av.NeighborList(8.5), -10).add_lj(
            av.LJ(0.0104, 3.4, 8.5)
        )
        sim.build()


class TestMDRun:
    def test_velocity_init_sets_velocities(self, argon_frame):
        sim = av.MD()
        sim.add_velocity_verlet(av.VelocityVerlet(0.001))
        sim.add_init_velocity(av.InitVelocity(10.0, 42))
        sim.build()
        sim.run(argon_frame, 1)
        # After at least 1 step, velocities should be nonzero
        vx = argon_frame.get_vx()
        assert not np.allclose(vx, 0.0), "Velocities should be initialized"

    def test_kinetic_energy_positive_after_init(self, argon_frame):
        sim = av.MD()
        sim.add_velocity_verlet(av.VelocityVerlet(0.001))
        sim.add_init_velocity(av.InitVelocity(10.0, 42))
        sim.build()
        sim.run(argon_frame, 1)
        assert argon_frame.kinetic > 0.0

    def test_positions_change_after_md(self, argon_frame):
        x_before = argon_frame.get_x().copy()
        sim = av.MD()
        sim.add_velocity_verlet(av.VelocityVerlet(0.001))
        sim.add_init_velocity(av.InitVelocity(10.0, 42))
        sim.build()
        sim.run(argon_frame, 10)
        x_after = argon_frame.get_x()
        assert not np.allclose(x_before, x_after), "Positions should move"

    def test_lj_md_energy_finite(self, argon_frame, full_lj_md):
        full_lj_md.build()
        full_lj_md.run(argon_frame, 5)
        assert np.isfinite(argon_frame.energy)
        assert np.isfinite(argon_frame.kinetic)

    def test_lj_md_forces_finite(self, argon_frame, full_lj_md):
        full_lj_md.build()
        full_lj_md.run(argon_frame, 1)
        fx = argon_frame.get_fx()
        assert np.all(np.isfinite(fx)), "Forces should be finite"

    def test_summary_returns_string(self, argon_frame):
        sim = av.MD()
        sim.add_velocity_verlet(av.VelocityVerlet(0.001))
        sim.add_init_velocity(av.InitVelocity(10.0, 42))
        sim.enable_timing()
        sim.build()
        sim.run(argon_frame, 1)
        summary = sim.summary()
        assert isinstance(summary, str)

    def test_run_zero_steps(self, argon_frame):
        sim = av.MD()
        sim.add_velocity_verlet(av.VelocityVerlet(0.001))
        sim.build()
        sim.run(argon_frame, 0)  # should not raise


class TestPhaseEnum:
    def test_all_phases_accessible(self):
        phases = [
            av.Phase.initial_integrate,
            av.Phase.pre_force,
            av.Phase.force,
            av.Phase.post_force,
            av.Phase.final_integrate,
            av.Phase.end_of_step,
        ]
        assert len(phases) == 6

    def test_phases_are_distinct(self):
        phases = [
            av.Phase.initial_integrate,
            av.Phase.pre_force,
            av.Phase.force,
            av.Phase.post_force,
            av.Phase.final_integrate,
            av.Phase.end_of_step,
        ]
        assert len(set(phases)) == len(phases)
