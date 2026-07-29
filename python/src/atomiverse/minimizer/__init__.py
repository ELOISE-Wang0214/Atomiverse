"""CPU geometry minimization — categorical setters matching C++.

Mirrors ``cpu::minimizer::Minimize`` (headers under ``cpu/minimize/``)::

    opt = Minimize()
    opt.set_optimizer(LBFGS(memory=5, alpha0=0.02, max_ls=8))
    opt.add_criterion(ForceConverged(ftol=0.5, …))  # optional
    opt.set_evaluator(PairLJCut(…))
    opt.set_neighbor(NeighborList(…))
    opt.build()
    result = opt.run(frame, max_iter)
"""

from __future__ import annotations

from atomiverse._binding import _cpp

_PHASE_MAP: dict[str, object] = {
    "search_direction": _cpp.MinimizePhase.search_direction,
    "line_search": _cpp.MinimizePhase.line_search,
    "force": _cpp.MinimizePhase.force,
    "convergence_check": _cpp.MinimizePhase.convergence_check,
    "end_of_step": _cpp.MinimizePhase.end_of_step,
}

_EVALUATOR_DISPATCH: dict[type, str] = {
    _cpp.PairLJCut: "set_evaluator_pair_lj_cut",
}

_PHASE_DISPATCH: dict[type, str] = {
    _cpp.PairLJCut: "add_pair_lj_cut",
    _cpp.NeighborList: "add_neighbor_list",
}


class PhaseBuilder:
    """Fluent builder returned by :meth:`Minimize.on` (escape hatch)."""

    def __init__(self, cpp_builder: object) -> None:
        self._cpp = cpp_builder

    def add(self, component: object, order: int = 0) -> PhaseBuilder:
        method_name = _PHASE_DISPATCH.get(type(component))
        if method_name is None:
            supported = ", ".join(t.__name__ for t in _PHASE_DISPATCH)
            raise TypeError(
                f"Cannot add {type(component).__name__} to a minimizer phase. "
                f"Supported: {supported}"
            )
        getattr(self._cpp, method_name)(component, order)
        return self


class Minimize:
    """Geometry-minimization driver for the CPU backend.

    Parameters
    ----------
    convergence : ForceConverged, optional
        Stop policy. If omitted, C++ ``install_defaults()`` attaches a
        default ``ForceConverged`` at build time.
    device : str, optional
        Execution backend. v1 accepts only ``"cpu"``.
    """

    def __init__(
        self,
        convergence: object | None = None,
        device: str = "cpu",
    ) -> None:
        if device != "cpu":
            raise ValueError(
                f"device={device!r}: GPU minimization is deferred; "
                "only device='cpu' is supported"
            )
        self.device = device
        self._cpp = _cpp.Minimize()
        self._built = False
        if convergence is not None:
            self.add_criterion(convergence)

    def set_optimizer(self, component: object) -> Minimize:
        """Install the optimizer (``Minimize::set_optimizer``)."""
        if type(component) is not _cpp.LBFGS:
            raise TypeError(
                f"set_optimizer expects LBFGS, got {type(component).__name__}"
            )
        self._cpp.set_optimizer(component)
        self._built = False
        return self

    def set_evaluator(self, component: object) -> Minimize:
        """Register an energy/force source (``Driver::set_evaluator``)."""
        method_name = _EVALUATOR_DISPATCH.get(type(component))
        if method_name is None:
            supported = ", ".join(t.__name__ for t in _EVALUATOR_DISPATCH)
            raise TypeError(
                f"set_evaluator does not support {type(component).__name__}. "
                f"Supported: {supported}"
            )
        getattr(self._cpp, method_name)(component)
        self._built = False
        return self

    def set_neighbor(self, component: object) -> Minimize:
        """Install the neighbor list (``Driver::set_neighbor``)."""
        if type(component) is not _cpp.NeighborList:
            raise TypeError(
                f"set_neighbor expects NeighborList, got "
                f"{type(component).__name__}"
            )
        self._cpp.set_neighbor(component)
        self._built = False
        return self

    def add_criterion(self, criterion: object) -> Minimize:
        """Attach a termination criterion (``Driver::add_criterion``)."""
        if type(criterion) is not _cpp.ForceConverged:
            raise TypeError(
                f"add_criterion expects ForceConverged, got "
                f"{type(criterion).__name__}"
            )
        self._cpp.add_criterion(criterion)
        self._built = False
        return self

    def on(self, phase: str) -> PhaseBuilder:
        """Escape hatch for explicit phase binding."""
        cpp_phase = _PHASE_MAP.get(phase)
        if cpp_phase is None:
            valid = ", ".join(f'"{p}"' for p in _PHASE_MAP)
            raise ValueError(f"Unknown minimizer phase {phase!r}. Valid: {valid}")
        self._built = False
        return PhaseBuilder(self._cpp.on(cpp_phase))

    def build(self) -> None:
        """Compile the execution plan (called automatically by ``run``)."""
        self._cpp.build()
        self._built = True

    def run(self, frame: object, max_iter: int) -> object:
        """Optimize ``frame`` for at most ``max_iter`` iterations.

        Returns a :class:`MinimizeResult` with loop-control metadata.
        Geometry / energy / forces are written back into *frame*.
        """
        if not self._built:
            self.build()
        return self._cpp.run(frame, max_iter)

    def __repr__(self) -> str:
        state = "built" if self._built else "not built"
        return f"Minimize({state})"


LBFGS = _cpp.LBFGS
ForceConverged = _cpp.ForceConverged
MinimizeResult = _cpp.MinimizeResult

__all__ = [
    "Minimize",
    "LBFGS",
    "ForceConverged",
    "MinimizeResult",
    "PhaseBuilder",
]
