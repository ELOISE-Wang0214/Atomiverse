---
slug: esmm-contract-04-forces-requested
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: driver / pipeline-contract
queue: esmm-campaign
deferred: false
depends_on: [esmm-contract-02]
unblocks: [rhf-in-mc, qmmm-in-mc, esmm-qmmm-03-electrostatic-rewrite]
grilled: true
---

# Spec: `esmm-contract-04-forces-requested` — energy-only phases, plumbed and typed

**Slug:** `esmm-contract-04-forces-requested`
**Domain:** CPU driver pipeline (force-request policy + MC warm-start guard)
**Queue:** esmm-campaign (Wave 1 — composition contract)
**Status:** approved. Grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

`LaunchContext::forces_requested` was designed for MC energy-only trials (its
own comment, `component.h:41-43`) and **never plumbed**: it defaults `true`
and no driver ever sets it false (`runtime.cpp:136, 210, 248`). Meanwhile
GPW/PW/AM1/RM1 evaluators and `HarmonicBias` already guard force assembly on
it. This spec compiles the policy into the Plan so MC's `sampling::Phase::energy`
runs every source energy-only — a QM engine in MC skips its entire gradient
for free.

Second deliverable: the **ASPC-under-MC guard**. Transaction rollback restores
geometry but cannot notify an evaluator, so a rejected move's density stays in
the ASPC history. Plain last-density reuse is correctness-neutral (SCF
converges to the geometry's fixed point to `tol` regardless of guess); the
**predictor** is the hazard — its corrector runs at loosened `corrector_tol_`
on an extrapolation of a time-continuous history that accept/reject jumps
violate, making energies history-dependent and breaking detailed balance.
Minimal correct thing: a typed throw.

Ships regression **R3** — the first QM-in-MC case in the tree.

## Domain basis

Policy carrier: the **Plan slot**, stamped once at build — no per-step
branching, no new Context concept for components (they keep reading
`ctx.forces_requested`). Nesting: an outer energy-only veto must survive into
a nested driver's own slots, so the AND rides a Context channel
(`ExecServices`).

Rejected alternatives: widening `DriverNode::run_nested_erased` with a
parameter (public-interface churn; a nested driver still needs the Context to
propagate into its own slots); an `OutcomePolicy → component` rollback
notification hook (cross-cutting machinery with one consumer, and the
predictor's loosened corrector would still bias energies even with
notification).

## Design

- `core/plan.h` `StageSlot`: `bool forces_requested = true;`
- `cpu/driver/driver.h`: new `protected virtual bool
  phase_wants_forces(PhaseEnum) const { return true; }`; at the end of
  `cpu::Driver::build()` (after `Base::build()`, `:214`):
  `for (auto& s : this->plan_.slots) s.forces_requested =
  phase_wants_forces(static_cast<PhaseEnum>(s.phase));`
- `cpu/component/context.h` `ExecServices`: `bool forces_requested = true;` —
  the nesting channel.
- `runtime.cpp` `step` / `step_plan` / `run_slot_`: construct each
  `LaunchContext` with `ctx.forces_requested = slot.forces_requested &&
  context_.exec.forces_requested`; around a nested-driver dispatch,
  save `context_.exec.forces_requested`, AND in the slot flag, restore after —
  `step_plan` reads the *nested* plan's slots, so the AND is what propagates
  the outer veto through arbitrary nesting depth.
- `mc_base.h`: `bool phase_wants_forces(sampling::Phase p) const override
  { return p != sampling::Phase::energy; }`
- `mean_field_method.cpp` `run_nested_erased`: `opts.forces = compute_force &&
  ctx.exec.forces_requested;` (drop the `(void)ctx`), plus the guard:
  `use_predictor_ && ctx.mc.proposal != nullptr` → `std::runtime_error`
  ("ASPC predictor is MD-only — rejected MC moves poison the density history
  and break detailed balance; disable predictor() for sampling").

Bootstrap interaction: none — `step_once(kBootstrapStep = −1)` goes through
`Runtime::step`, the slot flag applies identically (bootstrap is energy-only
too; the outcome policy seeds from `kEnergy`). Frequency-gated ops already
skip negative steps. GCMC insert/delete resizes `nlocal`; ZeroForces reads
`owned_count(ctx)` per launch — adapts.

## Files to create or modify

- `include/atomiverse/core/plan.h` — StageSlot flag
- `include/atomiverse/cpu/driver/driver.h` — virtual + build stamp
- `include/atomiverse/cpu/component/context.h` — ExecServices flag
- `src/cpu/driver/runtime.cpp` — the three construction sites + nested
  save/AND/restore
- `include/atomiverse/cpu/montecarlo/mc_base.h` — override
- `src/cpu/scf/mean_field_method.cpp` — `run_nested_erased` consult + ASPC
  guard (**gated file** — see Notes)
- `tests/cpu/driver/` — U4 slot-flag observation test
- `tests/cpu/scf/` — U6 ASPC-under-MC throw test
- NEW `regressions/mc/mc_bias_es/{reg_mc_bias_es.cpp,REFERENCE.md}` (R3)
- `regressions/CMakeLists.txt` — register R3 (**gated file**)

## Tasks

- [ ] RED (U4): stub observer bound at `sampling::Phase::energy` records
  `ctx.forces_requested == false`; the same stub at `md::Phase::force`
  records `true` (fails today: both true)
- [ ] RED (U6): `MeanFieldMethod` with `predictor()` invoked via
  `run_nested_erased` under a Context whose `mc.proposal` is non-null throws
  the typed error (fails today: silently proceeds)
- [ ] Implement plan.h / driver.h / context.h / runtime.cpp plumbing
- [ ] McBase override; MFM consult + guard
- [ ] GREEN both units + full ctest (MD/Minimize/Evaluate slots must all still
  report true — assert in U4)
- [ ] R3: `regressions/mc/mc_bias_es` — MC Metropolis at T with a nested RHF
  H₂ energy + `HarmonicBias` umbrella on the bond CV; sample the bond-length
  distribution; pin ⟨(r−r₀)²⟩ against the analytic stiff-spring limit
  k_BT/k_eff (bias-dominated regime; k_eff = k_bias + k_RHF with the RHF
  curvature at r₀ pinned offline) to the few-percent tolerance recorded in
  REFERENCE.md; ≤ 60 s
- [ ] REFERENCE.md: derive the stiff-spring variance analytically; pin the
  RHF/STO-3G curvature at r₀; state Metropolis + equipartition citations
  (Frenkel & Smit)
- [ ] Register R3 in regressions CMake (gated)

## Testing strategy

**Unit:** U4 (slot flag per phase, nested AND propagation via a stub nested
driver), U6 (typed guard). RED first, stubs only.

**Regression (R3):** exercises in one case: node-aware MC validation
(esmm-contract-02), energy-only trials (no SCF gradient — this is what keeps
it ≤ 60 s), per-sweep zeroing, and the ASPC guard staying off-path
(no predictor). Pin is analytic, not a snapshot.

## Out of scope

- Evaluator `=`→`+=` (esmm-contract-05) — R3's nested RHF is a sole source
  after ZeroEnergy, correct under either contract
- A rollback-notification hook for components (rejected; recorded here)
- Incremental/local MC energies

## Notes for implementers

- HARD GATE: before touching `src/cpu/scf/mean_field_method.cpp` or
  `regressions/CMakeLists.txt`, verify `git status` shows them clean. As of
  2026-07-26 the operator's libcint-removal/C++20 work (~760 files) is
  uncommitted and includes both. If dirty: STOP the spec and report — do not
  rebase over uncommitted operator work.
- Classical potentials ignore `forces_requested` (harmless — MC discards
  forces and per-sweep zeroing lands in esmm-contract-01). Do not "fix" them
  here.
- Keep the guard message actionable (name `predictor()` and the sampling
  context).
