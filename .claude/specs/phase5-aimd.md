---
title: Phase 5 — Ab initio molecular dynamics
status: approved
created: 2026-04-14
landed: 2026-05-06
---

# Phase 5 — AIMD

## Summary

AIMD = `cpu::md::MD` driver composed with `cpu::scf::SCFEvaluator` bound
to `md::Phase::force`. The integration loop is unchanged from classical
MD; only the force evaluator is different. There is no `AIMDDriver`
class and there is not going to be one.

This spec was held in draft until the Phase 1–4 prerequisite stack
matured enough that AIMD could be exercised against the same
`SCFEvaluator` that single-point and geometry-optimization workflows
use. With `cpu-driver-taxonomy` (status `done`) consolidating
`SCFEvaluator` as the unified force-providing component for all three
CPU driver categories, AIMD reduces to wiring + validation, not new
infrastructure.

The original gate (Phases 1–4 all `status: done`) was a precaution
against AIMD spawning a parallel god-class driver during early
implementation. That risk is now moot: `SCFEvaluator` already exists
and already composes cleanly with the MD driver. Phase 1's remaining
open items (f-shell support, additional default basis sets) are
orthogonal to AIMD on H₂ / H₂O / STO-3G workloads. Phases 3 (hybrid
DFT) and 4 (vibrational analysis) are independent science capabilities
that AIMD does not depend on. The gate is therefore released and the
spec is promoted to `approved`.

## Binding rules (load-bearing — do not remove)

- **No new driver.** AIMD uses the existing `cpu::md::MD` driver
  composed with `SCFEvaluator` at `md::Phase::force`. This composition
  already works today for classical MD; AIMD differs only in which
  force evaluator is bound.
- **No `AIMDDriver` / `AIMDRunner` / `AIMDWorkflow` class.** If the
  existing MD driver is insufficient, that is a bug in the MD driver —
  fix the driver, don't add a parallel class.
- **Thermostats, barostats, and integrators come from `cpu::md::*`
  components already in the tree.** If a specific ensemble needs a
  new component, it is added as a standalone `ComponentBase`, bound by
  the user.
- Reference comparison for accuracy/stability work: i-PI / CP2K
  reference trajectories or AIMD-derived thermodynamic averages from
  published benchmarks on the same system (documented per-test).
  PySCF does not ship a first-class AIMD driver, so PySCF parity is
  not a meaningful target for AIMD.

## Realized composition

```cpp
// NVE
cpu::scf::SCFEvaluator scf{config};
scf.predictor(/*order=*/3);

cpu::md::MD md;
md.add(cpu::md::VelocityVerlet{0.5});
md.on(md::Phase::force).add(std::move(scf));
md.build();
md.run(frame, n_steps);

// NVT (Nose-Hoover) — add one component
md.add(cpu::md::NoseHoover{T_target, tau, dt});

// NVT (Langevin) — replace the thermostat component
md.add(cpu::md::FixLangevin{T_target, damp, dt, seed});

// NPT — add NoseHoover + Berendsen Barostat
md.add(cpu::md::NoseHoover{T_target, tau_T, dt});
md.add(cpu::md::Barostat{P_target, tau_P, dt});
```

`SCFEvaluator::launch` reassembles the molecule via minimum-image
relative to atom 0 before each integral build, then restores the wrapped
positions before returning. This is the fix for the historical
"wrap × non-minimum-image SCF" hazard: `VelocityVerlet::wrap_orthorhombic`
wraps each atom independently into the primary box, but the SCF integral
builder reads raw `frame.atoms.x/y/z`. Without the unwrap, an atom that
wrapped across the boundary would appear ~one box length away from its
partners and the SCF would collapse onto the unphysical RHF ionic limit.
The unwrap assumes a single isolated fragment per SCF call, matching
current CPU SCF scope; multi-fragment SCF (Ewald, dimer / multimer
clusters in one cell) will need a generalised anchoring scheme.

## Acceptance — landed

| Criterion | Where |
|---|---|
| 1-step single-point and 5-step BOMD on H₂ for RHF / UHF / KS-LSDA / KS-PBE | `tests/cpu/scf/test_scf_evaluator.cpp` (`SCFEvaluatorTest`, `BOMDTest`) |
| Energy / momentum / bond-oscillation conservation under NVE | same file (`BOMDTest::RHF{EnergyConservation,BondOscillation,MomentumConservation}`) |
| ASPC predictor warm-start on H₂ NVE | same file (`BOMDTest::RHFWithASPC`) |
| NVT (Nose-Hoover) composition with SCFEvaluator | same file (`BOMDTest::RHFNoseHooverNVT`) |
| NVT (Langevin) composition — multiple Force-capable components | same file (`BOMDTest::RHFLangevinNVT`) |
| NPT (Berendsen Barostat) composition | same file (`BOMDTest::RHFBerendsenNPT`) |
| 200-step NVE energy drift bound on H₂ | `tests/cpu/scf/test_aimd_conservation.cpp::AIMDConservation::RHFNVEEnergyDrift200Steps` |
| ASPC state-leakage check (force from propagation = force from fresh single-point at same geometry, ≤ 1e-5 eV/Å) | same file (`AIMDConservation::ASPCMatchesFreshSinglePointMidTrajectory`) |
| Chained-builder rvalue safety (`SCFEvaluator{cfg}.predictor(3)` directly into `add()`) | `tests/cpu/scf/test_scf_evaluator.cpp::BOMDTest::ChainedBuilder*` |
| Wrap × non-min-image SCF recovery (atoms on box boundary, MD wrap → SCF unwraps via min-image) | `tests/cpu/scf/test_scf_evaluator.cpp::BOMDTest::BoundaryAtomsRecoverAfterWrap` |
| User-facing tutorial example | `examples/tutorial/06_aimd.cpp` |
| Per-MD-step wall-time measurement (NoPredict / ASPC3 / SinglePoint baseline) | `benchmarks/bench_aimd_step.cpp` |
| Tutorial documentation | `docs/tutorial/06-aimd.md` |

## Acceptance — deferred to a follow-up benchmark spec

These targets need offline-only runtime budgets and external reference
data, so they do not belong in CI:

- 10 ps liquid water NVE total-energy drift ≤ 1e-4 eV/atom/ps at dt = 0.5 fs.
- Pair-correlation function $g(r)$ for liquid water matching CPMD / CP2K
  references within typical reproducibility (~0.05 in peak height).
- Per-step AIMD wall time within 1.2× of a single fresh gradient on a
  production-scale system (water box, peptide). The
  `bench_aimd_step` target prints this ratio for $\mathrm{H}_2$, where
  the predictor's bookkeeping overhead dominates and the ratio is
  worse-than-1; this is a known small-system trait, not a regression
  in the predictor.

These should be captured in a separate `aimd-water-benchmark.md` (or
similar) when someone has the cycles to run them on representative
hardware. Cross-link this spec from there.

## Out of scope for this spec

- Path-integral MD (PIMD).
- Metadynamics / collective-variable dynamics.
- Constrained dynamics (SHAKE / RATTLE-style algebraic constraints).
- Nonadiabatic (mixed quantum-classical) dynamics.
- Open-shell AIMD via UHF — spin-contamination drift is an active
  research problem and does not belong in a first-pass AIMD spec.
- xyz + velocity restart from checkpoint. Not strictly AIMD-specific;
  belongs in a recorder/restart spec.

## Risks (live)

- **Drift from inconsistent gradients / energies.** Any inconsistency
  between SinglePoint energy and the gradient at the same geometry,
  even at the 1e-8 Ha level, will manifest as a slow non-conservation
  drift in NVE. The state-leakage check pins this in CI.
- **ASPC breakdown at large dt.** Predictor order vs step size is a
  tunable trade-off. Default is order 3 with corrector tol 1e-4; if
  the corrector exceeds its budget the evaluator falls back to a full
  SCF for that step (see `aspc_predictor.h`). Document tuning in the
  tutorial and in the small benchmark output.
- **PBC-wrap × non-minimum-image SCF.** Resolved by the
  in-`SCFEvaluator::launch` unwrap. Pinned by the
  `BoundaryAtomsRecoverAfterWrap` regression test. Future
  multi-fragment SCF (Ewald or multiple non-bonded molecules sharing
  one cell) will need a generalised anchoring scheme; the current
  anchor-on-atom-0 unwrap is single-fragment only.
