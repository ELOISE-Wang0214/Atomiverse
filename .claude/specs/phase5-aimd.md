---
title: Phase 5 — Ab initio molecular dynamics (placeholder)
status: draft
created: 2026-04-14
---

# Phase 5 — AIMD

## Summary

**Gated placeholder.** AIMD is out-of-scope for active work until Phases 1–4 are `status: done`. This document fixes the rules of engagement so that when the gate opens, the implementation does not regress into god-classes or bypass the Driver/Component model.

## Binding rules (apply here)

- **No new driver.** AIMD uses the existing `cpu::md::MD` driver composed with `SCFEvaluator` at `md::Phase::force`. This composition already works today for classical MD; AIMD differs only in which force evaluator is bound.
- **No `AIMDDriver` / `AIMDRunner` / `AIMDWorkflow` class.** If the existing MD driver is insufficient, that is a bug in the MD driver — fix the driver, don't add a parallel class.
- Thermostats, barostats, and integrators come from `cpu::md::*` components already in the tree. If a specific ensemble needs a new component, it is added as a standalone `ComponentBase`, bound by the user.
- PySCF does not ship a first-class AIMD driver. Ground truth becomes **i-PI / CP2K reference trajectories** or AIMD-derived thermodynamic averages from published benchmarks on the same system (documented per-test).

## Prerequisites (hard gate — all must be `done`)

1. **Phase 1** — benchmark suite, recorders, default level of theory (AIMD will reuse `TrajectoryWriter` unchanged).
2. **Phase 2** — `SinglePoint` driver (AIMD uses `SCFEvaluator` → `SinglePoint` internally per step).
3. **Phase 3** — at least PBE/def2-SVP (AIMD with hybrid DFT is optional; most use cases are pure functionals due to cost).
4. **Phase 4** — vibrational analysis available for validating integration stability (ZPE check, energy drift vs characteristic vibrational period).

## Expected composition (post-gate, for illustration only)

```cpp
// NVE — existing components, no new glue
cpu::md::MD md;
md.add(cpu::md::InitVelocity{/*T=*/300.0});
md.add(cpu::md::VelocityVerlet{/*dt=*/0.5 /*fs*/});
md.on(md::Phase::force).add(
    cpu::scf::SCFEvaluator{cpu::scf::default_config(integrals)}
        .predictor(/*order=*/3));   // ASPC warm-start existing
md.on(md::Phase::end_of_step)
   .add(cpu::recorder::TrajectoryWriter{"water_aimd.xyz", 10})
   .add(cpu::Thermo{100});
md.build();
md.run(frame, 1000);

// NVT — add NoseHoover component, already exists
md.on(md::Phase::post_force).add(cpu::md::NoseHoover{/*T=*/300.0, /*tau=*/100.0});
```

## Planned acceptance (deferred; finalized when phase opens)

- **Stability:** 10 ps water NVE simulation total-energy drift ≤ 1e-4 eV/atom/ps at dt=0.5 fs.
- **Accuracy:** pair correlation function g(r) for liquid water matches published CPMD/CP2K references within typical reproducibility (~0.05 in peak height).
- **Performance:** per-step AIMD wall time matches a non-SCF-dominated gradient cost (≤ 1.2× equivalent single gradient) — validates ASPC predictor efficacy.
- **Ground-truth check:** single-step energy + forces during AIMD propagation match a stand-alone `SinglePoint` + gradient on the same frame to ≤ 1e-8 Ha / 1e-7 Ha/Bohr (catches state leakage across steps).

## Expected scope (deferred; refine at gate)

- NVE, NVT (Nose-Hoover), NPT (Parrinello-Rahman via existing `Barostat`).
- Restart from an xyz + velocity checkpoint (new recorder component, if needed).
- ASPC predictor tuning guide for AIMD workloads.
- Not in scope for the first AIMD spec: path-integral MD, metadynamics, constrained dynamics, nonadiabatic coupling.

## Risks to plan for at gate time

- **Drift from inconsistent gradients / energies:** any inconsistency between `SinglePoint` energy and gradient (even at 1e-8 Ha level) is a showstopper. Phase 2 acceptance must pin this.
- **ASPC breakdown at large dt:** predictor order vs step size trade-off. Documented tuning in the AIMD spec.
- **Open-shell AIMD:** spin-contamination drift in UHF-AIMD is an active research problem. Out of scope for first implementation.

## Exit gate

Spec gets promoted from `draft` to `approved` only when Phases 1–4 are `done`. At that point, this file is rewritten in full with concrete acceptance numbers.
