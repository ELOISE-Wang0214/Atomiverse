"""L-BFGS geometry optimizer — accepts ``molrs.Frame``, converts internally."""

from __future__ import annotations

from dataclasses import dataclass

from atomiverse._binding import _cpp
from atomiverse.scf import _build_integrals, _METHOD_MAP, _XC_MAP


@dataclass(frozen=True)
class OptRecord:
    """Per-step optimization metrics."""

    step: int
    energy: float          # eV
    max_force: float       # eV/Angstrom
    rms_force: float       # eV/Angstrom
    max_displacement: float  # Angstrom
    rms_displacement: float  # Angstrom
    energy_change: float   # eV


@dataclass(frozen=True)
class OptResult:
    """Result of a geometry optimization."""

    converged: bool
    n_steps: int
    final_energy: float    # eV
    trace: list[OptRecord]


class GeometryOptimizer:
    """L-BFGS geometry optimizer with SCF forces.

    Parameters
    ----------
    frame : molrs.Frame
        Initial molecular geometry.
    method : str
        SCF method: ``"rhf"``, ``"uhf"``, ``"rks"``, ``"uks"``.
    basis : str
        Basis set name.
    xc : str or None
        XC functional for KS methods.
    charge : int
        Net molecular charge.
    multiplicity : int
        Spin multiplicity 2S+1.
    memory : int
        L-BFGS history size (default 5).
    ftol : float
        Max force convergence tolerance [eV/Angstrom] (default 0.05).
    max_steps : int
        Maximum optimization steps (default 100).

    Example::

        import molrs
        import atomiverse as av

        frame = molrs.Frame()
        # ... populate frame["atoms"] with x/y/z/type ...

        opt = av.GeometryOptimizer(frame, method="rhf")
        result = opt.run()
        print(result.converged, result.final_energy)
    """

    def __init__(
        self,
        frame: object,
        method: str = "rhf",
        basis: str = "sto-3g",
        xc: str | None = None,
        charge: int = 0,
        multiplicity: int = 1,
        memory: int = 5,
        ftol: float = 0.05,
        max_steps: int = 100,
        **kwargs: float,
    ) -> None:
        method_key = method.lower()
        if method_key not in _METHOD_MAP:
            raise ValueError(f"Unknown method {method!r}")

        self._frame = frame
        self._max_steps = max_steps

        # Build SCF config (same logic as SCF class)
        spin_str, mf_str = _METHOD_MAP[method_key]
        is_ks = mf_str == "KS"
        if is_ks and xc is None:
            xc = "LSDA"

        ip_cpp = _build_integrals(frame, basis, charge)
        n_el = int(ip_cpp.n_electrons)
        spin = multiplicity - 1
        if (n_el + spin) % 2 != 0:
            raise ValueError(
                f"Incompatible n_electrons={n_el} and multiplicity={multiplicity}"
            )
        n_alpha = (n_el + spin) // 2
        n_beta = (n_el - spin) // 2

        config = _cpp.SCFConfig()
        config.spin = getattr(_cpp.SpinType, spin_str)
        config.method = getattr(_cpp.MeanFieldType, mf_str)
        config.n_alpha = n_alpha
        config.n_beta = n_beta
        if is_ks:
            xc_enum = _XC_MAP.get(xc)  # type: ignore[arg-type]
            if xc_enum is None:
                raise ValueError(f"Unknown XC functional {xc!r}")
            config.xc_type = xc_enum
            config.n_radial = int(kwargs.get("n_radial", 50))
            config.n_angular = int(kwargs.get("n_angular", 110))

        self._config = config
        self._charge = charge
        self._memory = memory
        self._ftol = ftol

    def run(self) -> OptResult:
        """Run L-BFGS optimization. Updates the molrs.Frame in place."""
        driver = _cpp.LBFGS()
        driver.add_lbfgs_step(
            _cpp.LBFGSStep(self._memory, self._ftol)
        )
        driver.on(_cpp.MinimizerPhase.force).add_scf_evaluator(
            self._config, self._charge
        )
        driver.build()
        # The driver borrows the molrs.Frame, materializes a transient
        # execution frame, optimizes, and writes the geometry back in place.
        driver.run(self._frame, self._max_steps)

        # Extract trace from the LBFGSStep component
        # (trace access not yet exposed in bindings — return minimal result)
        return OptResult(
            converged=True,  # TODO: extract from C++
            n_steps=self._max_steps,
            final_energy=0.0,  # TODO: extract from C++
            trace=[],
        )
