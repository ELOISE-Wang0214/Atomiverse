"""Atomiverse: molecular simulation compute backend for molrs.

Atomiverse binds drivers and components only — it does not bind a Frame
data type. Simulation state crosses the boundary as a ``molrs.Frame``;
each driver ``run()`` borrows it and writes results back in place.

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
``neighbor.NeighborList``             ``NeighborList``            ``cpu::md::NeighborList<>`` (neighbor/)
``potentials.pair.LJCut``             ``PairLJCut``               ``cpu::PairLJCut`` (potentials/)
``potentials.bond.Harmonic``          ``BondHarmonic``            ``cpu::BondHarmonic`` (potentials/)
``potentials.angle.Harmonic``         ``AngleHarmonic``           ``cpu::AngleHarmonic`` (potentials/)
``potentials.dihedral.MultiHarmonic`` ``DihedralMultiHarmonic``   ``cpu::DihedralMultiHarmonic`` (potentials/)
``recorder.Thermo``                   ``Thermo``                  ``cpu::Thermo`` (recorder/)
``pw.pw_scf`` / ``pw.PWResult``       ``pw_scf`` / ``PWResult``   ``cpu::pw`` (pw/)
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
- **AM1 semiempirical** — single-point energy [eV] + Mulliken charges [e]
  via ``single_point(frame, AM1Evaluator(...))``; no AM1 forces yet.
- **Plane-wave DFT** — ``pw_scf`` / ``PWResult`` (requires a wheel built
  with ``ATV_BUILD_PW=ON``).

Roadmap: Gaussian-SCF and geometry-optimization bindings are future work
with their own spec — deliberately NOT exposed rather than exposed broken.

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

from atomiverse.driver import single_point
from atomiverse.md import MD, InitVelocity
from atomiverse.md.integrator import VelocityVerlet
from atomiverse.md.thermostat import NoseHoover
from atomiverse.neighbor import NeighborList
from atomiverse.potentials.angle import Harmonic as AngleHarmonic
from atomiverse.potentials.bond import Harmonic as BondHarmonic
from atomiverse.potentials.dihedral import MultiHarmonic as DihedralMultiHarmonic
from atomiverse.potentials.pair import LJCut as PairLJCut
from atomiverse.pw import PWResult, pw_scf
from atomiverse.recorder import Thermo
from atomiverse.semiempirical import AM1Config, AM1Evaluator, AM1Result

__all__ = [
    "MD",
    "VelocityVerlet",
    "InitVelocity",
    "NoseHoover",
    "PairLJCut",
    "NeighborList",
    "BondHarmonic",
    "AngleHarmonic",
    "DihedralMultiHarmonic",
    "Thermo",
    "pw_scf",
    "PWResult",
    "AM1Config",
    "AM1Evaluator",
    "AM1Result",
    "single_point",
]
