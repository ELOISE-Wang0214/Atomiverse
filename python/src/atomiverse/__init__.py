"""Atomiverse: molecular simulation compute backend for molrs.

Atomiverse binds drivers and components only — it does not bind a Frame
data type. Simulation state crosses the boundary as a ``molrs.Frame``;
each driver ``run()`` borrows it and writes results back in place. Global
results use the exact typed ``molrs.Frame.meta`` v2 contract; Atomiverse does
not synthesize ``scalars`` or ``stress`` blocks.

Namespaces
----------
The DEEP paths are canonical — they mirror the C++ backend directory
layering one-to-one. The top level re-exports the same objects under
their LAMMPS-style compound names as a convenience (identity, not copies:
``atomiverse.PairLJCut is atomiverse.potentials.pair.LJCut``).

=====================================  ==========================  =========================
Canonical deep path                   Top-level re-export         C++ symbol (directory)
=====================================  ==========================  =========================
``md.MD``                             ``MD``                      ``cpu::md::MD`` (md/)
``md.integrator.VelocityVerlet``      ``VelocityVerlet``          ``cpu::md::VelocityVerlet`` (md/)
``md.thermostat.NoseHoover``          ``NoseHoover``              ``cpu::md::NoseHoover`` (md/)
``md.InitVelocity``                   ``InitVelocity``            ``cpu::md::InitVelocity`` (md/)
``minimizer.Minimize``                ``Minimize``                ``cpu::minimizer::Minimize`` (minimizer/)
``minimizer.LBFGS``                   ``LBFGS``                   ``cpu::minimizer::LBFGS`` (minimizer/)
``neighbor.NeighborList``             ``NeighborList``            ``cpu::md::NeighborList<>`` (neighbor/)
``potentials.pair.LJCut``             ``PairLJCut``               ``cpu::PairLJCut`` (potentials/)
``potentials.bond.Harmonic``          ``BondHarmonic``            ``cpu::BondHarmonic`` (potentials/)
``potentials.angle.Harmonic``         ``AngleHarmonic``           ``cpu::AngleHarmonic`` (potentials/)
``potentials.dihedral.MultiHarmonic`` ``DihedralMultiHarmonic``   ``cpu::DihedralMultiHarmonic`` (potentials/)
``recorder.Thermo``                   ``Thermo``                  ``cpu::Thermo`` (recorder/)
``pw.PWSolver`` / ``pw.PWConfig``     ``PWSolver`` / ``PWConfig`` ``cpu::pw`` (pw/)
``pw.PWResult``                       ``PWResult``                 ``cpu::pw`` (pw/)
``semiempirical.AM1Config``           ``AM1Config``               ``cpu::semiempirical`` (semiempirical/)
``semiempirical.AM1Evaluator``        ``AM1Evaluator``            ``cpu::semiempirical`` (semiempirical/)
``semiempirical.AM1Result``           ``AM1Result``               ``cpu::semiempirical`` (semiempirical/)
``driver.single_point``               ``single_point``            ``cpu::single_point`` (driver/)
=====================================  ==========================  =========================

Backend selection is a ``device=`` constructor parameter (PyTorch-style),
never a namespace axis: ``MD(device="cpu")`` today; ``"cuda"`` raises
``ValueError`` until the CUDA bindings land (their own spec).

Current surface (honest capability boundary)
--------------------------------------------
- **MD** — the ``MD`` driver plus its components (table above).
- **Geometry optimization** — CPU ``Minimize`` + ``LBFGS`` with an
  independent ``ForceConverged`` policy.
- **AM1 semiempirical** — single-point energy [eV] + Mulliken charges [e]
  via ``single_point(frame, AM1Evaluator(...))``; no AM1 forces yet.
- **Plane-wave DFT** — explicit ``PWConfig`` / ``PWSolver`` / ``PWResult``
  lifecycle (present only in wheels built with ``ATV_BUILD_PW=ON``).

Roadmap: Gaussian-SCF bindings remain a later CPU breadth slice.

Example — Molecular Dynamics (canonical deep paths)::

    import atomiverse as atv

    sim = atv.md.MD(dt=0.001)
    sim.add(atv.md.InitVelocity(10.0, 42))
    sim.on("force").add(atv.neighbor.NeighborList(8.5), order=-10) \\
                   .add(atv.potentials.pair.LJCut(0.0104, 3.4, 8.5))
    sim.run(frame, 100)

Example — AM1 single point::

    ev = atv.AM1Evaluator(atv.AM1Config(charge=0))
    atv.single_point(frame, ev)
    print(ev.last_result().mulliken_charges)
"""

from atomiverse._binding import _cpp
from atomiverse.driver import single_point
from atomiverse.md import MD, InitVelocity
from atomiverse.md.integrator import VelocityVerlet
from atomiverse.md.thermostat import NoseHoover
from atomiverse.minimizer import ForceConverged, LBFGS, Minimize, MinimizeResult
from atomiverse.neighbor import NeighborList
from atomiverse.potentials.angle import Harmonic as AngleHarmonic
from atomiverse.potentials.bond import Harmonic as BondHarmonic
from atomiverse.potentials.dihedral import MultiHarmonic as DihedralMultiHarmonic
from atomiverse.potentials.pair import LJCut as PairLJCut
from atomiverse.recorder import Thermo
from atomiverse.semiempirical import AM1Config, AM1Evaluator, AM1Result

__all__ = [
    "MD",
    "VelocityVerlet",
    "InitVelocity",
    "NoseHoover",
    "Minimize",
    "LBFGS",
    "ForceConverged",
    "MinimizeResult",
    "PairLJCut",
    "NeighborList",
    "BondHarmonic",
    "AngleHarmonic",
    "DihedralMultiHarmonic",
    "Thermo",
    "AM1Config",
    "AM1Evaluator",
    "AM1Result",
    "single_point",
]

if hasattr(_cpp, "PWSolver"):
    from atomiverse.pw import PWConfig as PWConfig
    from atomiverse.pw import PWResult as PWResult
    from atomiverse.pw import PWSolver as PWSolver

    __all__.extend(["PWConfig", "PWSolver", "PWResult"])
