"""Python integration tests for the MD driver running on a molrs.Frame.

Covers rust-frame-lowering-03 ac-002: a Python MD run starts from a
molrs.Frame end to end — the driver borrows the molrs.Frame, advances it,
and writes positions / velocities / energy back into that same frame.

Also asserts ac-001's structural deliverable from the Python side: no
CpuFrame type is exposed and atomiverse._convert no longer exists.
"""

import importlib

import numpy as np
import pytest
import molrs

import atomiverse as av


# ── Component construction ──────────────────────────────────────────


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


# ── Driver build ────────────────────────────────────────────────────


class TestMDDriverBuild:
    def test_dt_shorthand_builds(self, full_lj_md):
        full_lj_md.build()

    def test_phase_builder_chain(self):
        sim = av.MD(dt=0.001)
        sim.on("force").add(av.NeighborList(8.5), order=-10) \
                       .add(av.LJ(0.0104, 3.4, 8.5))
        sim.build()

    def test_invalid_phase_raises(self):
        sim = av.MD()
        with pytest.raises(ValueError, match="Unknown phase"):
            sim.on("nonexistent")

    def test_repr(self):
        sim = av.MD(dt=0.001)
        assert "not built" in repr(sim)
        sim.on("force").add(av.NeighborList(8.5), order=-10) \
                       .add(av.LJ(0.0104, 3.4, 8.5))
        sim.build()
        assert "built" in repr(sim)


# ── ac-002: MD on a molrs.Frame, end to end ─────────────────────────


class TestMDRunOnMolrsFrame:
    def test_positions_written_back(self, argon_frame, full_lj_md):
        """Positions update in the SAME molrs.Frame after run()."""
        x_before = np.asarray(argon_frame["atoms"]["x"]).copy()
        full_lj_md.run(argon_frame, 10)
        x_after = np.asarray(argon_frame["atoms"]["x"])
        assert not np.allclose(x_before, x_after), \
            "positions should move and be written back into the molrs.Frame"

    def test_velocities_written_back(self, argon_frame, full_lj_md):
        """InitVelocity-seeded velocities land back in the molrs.Frame."""
        full_lj_md.run(argon_frame, 1)
        vx = np.asarray(argon_frame["atoms"]["vx"])
        assert np.any(vx != 0.0), \
            "velocities should be initialized and written back"

    def test_energy_finite(self, argon_frame, full_lj_md):
        """Run yields a finite total energy in the molrs.Frame."""
        full_lj_md.run(argon_frame, 5)
        assert "scalars" in argon_frame.keys()
        energy = float(np.asarray(argon_frame["scalars"]["energy"])[0])
        kinetic = float(np.asarray(argon_frame["scalars"]["kinetic"])[0])
        assert np.isfinite(energy)
        assert np.isfinite(kinetic)
        assert kinetic > 0.0

    def test_forces_finite(self, argon_frame, full_lj_md):
        """Forces written back are finite."""
        full_lj_md.run(argon_frame, 1)
        fx = np.asarray(argon_frame["atoms"]["fx"])
        assert np.all(np.isfinite(fx))

    def test_auto_build(self, argon_frame, full_lj_md):
        """run() compiles the plan automatically."""
        # full_lj_md is not explicitly built; run() must build it.
        full_lj_md.run(argon_frame, 1)

    def test_run_zero_steps(self, argon_frame, full_lj_md):
        full_lj_md.run(argon_frame, 0)

    def test_missing_required_field_raises(self):
        """A molrs.Frame without atom positions fails loud."""
        frame = molrs.Frame()
        sim = av.MD(dt=0.001)
        sim.on("force").add(av.NeighborList(8.5), order=-10) \
                       .add(av.LJ(0.0104, 3.4, 8.5))
        with pytest.raises(RuntimeError, match="atoms"):
            sim.run(frame, 1)

    def test_non_frame_argument_raises(self, full_lj_md):
        """Passing a non-molrs.Frame object fails with a clear message."""
        with pytest.raises(RuntimeError, match="molrs.Frame"):
            full_lj_md.run(object(), 1)


# ── ac-001: no CpuFrame / no _convert on the Python surface ─────────


class TestFrameBindingRemoved:
    def test_no_cpuframe_on_extension(self):
        from atomiverse._binding import _cpp
        assert not hasattr(_cpp, "CpuFrame")

    def test_no_frame_attr_on_package(self):
        # The old `atomiverse.Frame = _cpp.CpuFrame` alias is gone.
        assert not hasattr(av, "Frame")

    def test_convert_module_gone(self):
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module("atomiverse._convert")
