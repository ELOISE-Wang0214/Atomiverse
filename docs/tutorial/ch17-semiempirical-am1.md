# Chapter 17: Semiempirical AM1 — Cheap Quantum Forces

Gaussian SCF and GPW give transferable electronic structure at a price.
Classical force fields are fast but parameterized. **Semiempirical** methods
sit in the middle: a minimal valence electronic Hamiltonian with element
parameters fitted once, then reused for organic chemistry at a cost between
classical MD and RHF.

Atomiverse ships **AM1** (Austin Model 1) on the **NDDO** stack. The
architecture is the same as every other method in this book:

```text
Preset (Evaluate / MD / Minimize / Hessian)
  └── AM1Evaluator  :  cpu::Semiempirical   (kind)
         └── AM1Solver : NDDOSolver         (inner SCF + forces)
```

There is **no** `SemiempiricalDriver`. DFTB and xTB, if added later, get their
own Solver + Evaluator + kind — not a seat on `NDDOSolver`.

The example to read alongside this chapter is
[`examples/15_am1.cpp`](https://github.com/MolCrafts/Atomiverse/blob/dev/examples/15_am1.cpp).

## What AM1 Computes

| Quantity | Units | Where |
|----------|-------|--------|
| Total energy \(E_\mathrm{elec}+E_\mathrm{nuc}\) | eV | `frame.meta.get(kEnergy)` |
| Heat of formation \(\Delta H_f\) (MOPAC-style) | kcal/mol | `AM1Result::heat_of_formation_kcal_mol` |
| Mulliken-like charges | e | `kAtom_charge` |
| Nuclear forces | eV/Å | `kAtom_fx/fy/fz` when `with_gradient(true)` |

Heat of formation follows the openmopac definition:

\[
\Delta H_f = (E_\mathrm{tot} - \sum_A e_\mathrm{isol}(A))\cdot 23.061
  + \sum_A E_\mathrm{atom}(A)
\]

with \(e_\mathrm{isol}\) from the AM1 parameter table and \(E_\mathrm{atom}\)
the experimental gas-phase atomic heats at 298 K.

Forces are the **SCF analytic gradient**: at a stationary density \(P\),
displace nuclei at fixed \(P\) and differentiate the energy functional
(Hellmann–Feynman for variational NDDO). They match re-SCF energy finite
differences; regressions live under `regressions/semiempirical/`.

## Single-point energy and forces

```cpp
#include "atomiverse/cpu/driver/evaluate.h"
#include "atomiverse/cpu/semiempirical/am1_evaluator.h"
#include "atomiverse/facade/stage.h"

// Frame: type = atomic number Z, positions [Å], cubic box for MD later.
cpu::Frame frame = /* water or small organics */;

cpu::Evaluate eval;
eval.on(evaluate::Phase::evaluate)
    .add(cpu::semiempirical::AM1Evaluator{}.with_gradient(true));
eval.build();
eval.run(frame);

double E_eV = frame.meta.get(ATV::kEnergy);
// charges in kAtom_charge; forces in kAtom_fx/fy/fz
```

`default_bindings()` is empty on purpose: like SCF and PW, you **must** bind
with `.on(phase).add(...)`. Energy-only single points:

```cpp
eval.on(evaluate::Phase::evaluate)
    .add(cpu::semiempirical::AM1Evaluator{}.with_gradient(false));
```

## Geometry optimization

Same `Minimize` preset as classical or RHF — only the force component changes:

```cpp
cpu::minimizer::Minimize opt;
opt.set_optimizer(cpu::minimizer::LBFGS{5, 0.02, 8});
opt.add_handler(cpu::ForceConverged{0.5, 0.3, 0.5, 0.3, 1e-3});
opt.on(minimizer::Phase::force)
    .add(cpu::semiempirical::AM1Evaluator{}.with_gradient(true));
opt.build();
opt.run(frame, 50);
// opt.converged(), opt.stop_reason()
```

Pin: `regressions/semiempirical/am1_water_minimize`.

## Born–Oppenheimer MD

```cpp
cpu::md::MD sim;
sim.set_integrator(cpu::md::VelocityVerlet{/*dt [fs]=*/0.25});
sim.add(cpu::md::InitVelocity{/*T [K]=*/10.0, /*seed=*/7u});
sim.on(md::Phase::force)
    .add(cpu::semiempirical::AM1Evaluator{}.with_gradient(true));
sim.build();
sim.run(frame, 20);
```

**Molecular unwrap.** AM1 is not a periodic electronic method. MD may wrap
atoms across the box; `AM1Evaluator` re-assembles an isolated fragment
(anchor on atom 0, min-image) before the SCF, then restores wrapped
coordinates — the same idea as Gaussian SCF AIMD. Always set physical masses
on `kAtom_mass` (amu).

Pin: `regressions/semiempirical/am1_water_bomd`.

## Inner solver (no Driver)

For charges or ΔHf without a plan:

```cpp
cpu::semiempirical::AM1Solver solver;
solver.build(frame);
const auto& res = solver.run(frame);
// res.total_energy_ev, res.heat_of_formation_kcal_mol, res.mulliken_charges
solver.compute_forces(frame); // fills res.forces_ev_a [eV/Å]
```

## Limits (honest)

- Closed-shell RHF-shaped AM1 only; s/p valence (d-orbital atoms throw).
- Parameter set is AM1 only (no PM6/PM7 in-tree yet — same `NDDOSolver` later).
- Hessian = finite difference of these forces via `cpu::vib::Hessian`, not
  analytic second derivatives.
- Force evaluation reuses the SCF density (warm start across MD/opt steps).
- DFTB / xTB are **not** this chapter — separate method families.

## Where next

- Classical many-body metals/solids: [Chapter 16](ch16-many-body-potentials.md).
- Swap AM1 for SCF on the same Minimize/MD plan: [Chapter 15](ch15-cross-method.md).
- Roadmap census: [Semiempirical](../roadmap/semiempirical.md).
