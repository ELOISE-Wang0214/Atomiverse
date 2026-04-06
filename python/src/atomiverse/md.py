"""OOP Molecular Dynamics driver — accepts ``mp.Frame``, converts internally."""

from __future__ import annotations

import molpy as mp

from atomiverse._binding import _cpp
from atomiverse._convert import mp_to_cpu, cpu_to_mp

# ---------------------------------------------------------------------------
# Phase name → C++ enum
# ---------------------------------------------------------------------------
_PHASE_MAP: dict[str, object] = {
    "initial_integrate": _cpp.Phase.initial_integrate,
    "pre_force": _cpp.Phase.pre_force,
    "force": _cpp.Phase.force,
    "post_force": _cpp.Phase.post_force,
    "final_integrate": _cpp.Phase.final_integrate,
    "end_of_step": _cpp.Phase.end_of_step,
}

# ---------------------------------------------------------------------------
# Component dispatch tables
# ---------------------------------------------------------------------------
_DRIVER_DISPATCH: dict[type, str] = {
    _cpp.VelocityVerlet: "add_velocity_verlet",
    _cpp.InitVelocity: "add_init_velocity",
    _cpp.NoseHoover: "add_nose_hoover",
}

_PHASE_DISPATCH: dict[type, str] = {
    _cpp.LJ: "add_lj",
    _cpp.NeighborList: "add_neighbor_list",
    _cpp.HarmonicBond: "add_harmonic_bond",
    _cpp.HarmonicAngle: "add_harmonic_angle",
    _cpp.MultiHarmonicDihedral: "add_multi_harmonic_dihedral",
    _cpp.Thermo: "add_thermo",
}


class PhaseBuilder:
    """Fluent builder returned by :meth:`MD.on`.

    Example::

        sim.on("force").add(av.NeighborList(cutoff=8.5), order=-10) \\
                       .add(av.LJ(epsilon=0.01, sigma=3.4, cutoff=8.5))
    """

    def __init__(self, cpp_builder: object) -> None:
        self._cpp = cpp_builder

    def add(self, component: object, order: int = 0) -> PhaseBuilder:
        """Add *component* at this phase with the given execution *order*."""
        method_name = _PHASE_DISPATCH.get(type(component))
        if method_name is None:
            supported = ", ".join(t.__name__ for t in _PHASE_DISPATCH)
            raise TypeError(
                f"Cannot add {type(component).__name__} to a phase. "
                f"Supported: {supported}"
            )
        getattr(self._cpp, method_name)(component, order)
        return self


class MD:
    """Molecular dynamics driver (CPU backend).

    Parameters
    ----------
    dt : float, optional
        Timestep in fs.  If provided, a ``VelocityVerlet`` integrator
        is added automatically.

    Example::

        import molpy as mp
        import atomiverse as av

        frame = mp.Frame({"atoms": mp.Block({
            "x": x, "y": y, "z": z, "mass": m, "type": t,
        })})
        frame.box = mp.Box(np.eye(3) * 30.0)

        sim = av.MD(dt=0.001)
        sim.add(av.InitVelocity(temperature=10.0, seed=42))
        sim.on("force").add(av.NeighborList(cutoff=8.5), order=-10)
        sim.on("force").add(av.LJ(epsilon=0.01, sigma=3.4, cutoff=8.5))
        sim.run(frame, 1000)
    """

    def __init__(self, dt: float | None = None) -> None:
        self._cpp = _cpp.MD()
        self._built = False
        if dt is not None:
            self.add(_cpp.VelocityVerlet(dt))

    # ------------------------------------------------------------------
    # Component registration
    # ------------------------------------------------------------------

    def add(self, component: object) -> MD:
        """Add a component that has default phase bindings.

        Supported: ``VelocityVerlet``, ``InitVelocity``, ``NoseHoover``.
        For phase-specific components use :meth:`on`.
        """
        method_name = _DRIVER_DISPATCH.get(type(component))
        if method_name is None:
            supported = ", ".join(t.__name__ for t in _DRIVER_DISPATCH)
            raise TypeError(
                f"Cannot add {type(component).__name__} directly. "
                f"Use sim.on(phase).add() instead.  "
                f"Driver-level components: {supported}"
            )
        getattr(self._cpp, method_name)(component)
        self._built = False
        return self

    def on(self, phase: str) -> PhaseBuilder:
        """Return a :class:`PhaseBuilder` for explicit phase binding.

        *phase* is one of: ``"initial_integrate"``, ``"pre_force"``,
        ``"force"``, ``"post_force"``, ``"final_integrate"``,
        ``"end_of_step"``.
        """
        cpp_phase = _PHASE_MAP.get(phase)
        if cpp_phase is None:
            valid = ", ".join(f'"{p}"' for p in _PHASE_MAP)
            raise ValueError(
                f"Unknown phase {phase!r}. Valid phases: {valid}"
            )
        self._built = False
        return PhaseBuilder(self._cpp.on(cpp_phase))

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

    def build(self) -> None:
        """Compile the execution plan (called automatically by ``run``)."""
        self._cpp.build()
        self._built = True

    def run(self, frame: mp.Frame, n_steps: int) -> None:
        """Run *n_steps* of MD on *frame*.

        Converts ``mp.Frame`` to internal representation, runs the
        simulation, then writes positions/velocities/forces/energy back
        into *frame* in-place.

        Raises
        ------
        ValueError
            If *frame* is missing required fields (atoms with x/y/z/type).
            If ``mass`` is not set in the atoms block (required for MD).
        """
        if "atoms" in frame and "mass" not in frame["atoms"]:
            raise ValueError(
                "MD requires 'mass' in atoms block. "
                "Set frame['atoms']['mass'] before calling run()."
            )

        cf = mp_to_cpu(frame)
        if not self._built:
            self.build()
        self._cpp.run(cf, n_steps)
        cpu_to_mp(cf, frame)

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    def enable_timing(self) -> None:
        """Enable per-component timing."""
        self._cpp.enable_timing()

    def summary(self) -> str:
        """Return a timing summary string."""
        return self._cpp.summary()

    def __repr__(self) -> str:
        state = "built" if self._built else "not built"
        return f"MD({state})"
