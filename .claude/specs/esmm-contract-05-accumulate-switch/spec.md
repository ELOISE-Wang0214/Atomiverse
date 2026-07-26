---
slug: esmm-contract-05-accumulate-switch
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: driver / pipeline-contract
queue: esmm-campaign
deferred: false
depends_on: [esmm-contract-01]
unblocks: [es-mm-cobinding, dual-es-cobinding, esmm-qmmm-03-electrostatic-rewrite, esmm-cpmd-01-electronic-home]
grilled: true
---

# Spec: `esmm-contract-05-accumulate-switch` — one accumulation contract for every source

**Slug:** `esmm-contract-05-accumulate-switch`
**Domain:** CPU force pipeline (the `=`→`+=` unification + `kEnergyElectronic`)
**Queue:** esmm-campaign (Wave 1 — composition contract)
**Status:** approved. Grilled via plan-mode design review + operator Q&A, 2026-07-26.

## Summary

The engine's composition contract is split: classical
`Potential`/`Bias`/DFTB2 accumulate into `kEnergy`/`kAtom_f*`/`kStress`, while
**five** electronic-structure evaluators overwrite — GPW (7 branches), MFM,
PW, AM1, RM1. Two ES sources clobber each other; ES+MM co-binding silently
drops whichever ran first. This is the root blocker of the operator's
"electronic structure × molecular mechanics 自由组合" directive. This spec
flips every in-plan source to `+=`, adds the `kEnergyElectronic` per-term key
so the energy decomposition invariant heals, and ships regression **R1** (the
first ES+MM co-bound case).

## Domain basis — the contract (binding, campaign-wide)

Two channel classes:

- **Accumulator channels** — `kEnergy`, every `kEnergy*` per-term key,
  `kAtom_fx/fy/fz`, `kStress` (virial part). *An in-plan `launch()` /
  `run_nested_erased()` source only ever `+=`; the hosting driver owns
  zeroing* (landed by `esmm-contract-01`).
- **State channels** — sole-producer, last-writer-wins `=`: `kAtom_charge`
  (AM1/RM1 Mulliken), `kStressKinetic` (Virial), `kDipole`, wavefunction
  members.

**Standalone entries keep overwrite**: `MeanFieldMethod::run(Frame&)`,
`forces()`, `mp2_forces()`, `ccsd_forces()` — "the mutating entry always
writes" is their documented Frame-out contract, and overwrite *is* the correct
self-zeroing for a sole author. Rejected alternative: standalone zero-then-`+=`
(changes public API behavior for no gain — a standalone call cannot have a
co-author).

Equivalence argument: after a driver's Zero*, `+=` ≡ `=` for every sole-source
configuration. That equivalence holds in MD/Evaluate/Minimize/Phonon/MC today
and in Hessian/PathDriver/stress **only after** `esmm-contract-01` — hence the
dependency. All existing tests/regressions must pass unchanged.

## Design — exact edit sites

| File | Change |
| ---- | ------ |
| `src/cpu/scf/mean_field_method.cpp` (+ private `DriveOpts` in the header, `:730` block) | `DriveOpts::accumulate = false` flag; `drive_scf_core_` epilogue (`:958-985`) branches `=` vs `+=` on it; `run_nested_erased` (`:283-294`) sets `accumulate = true`. `mp2_forces`/`ccsd_forces` untouched. |
| `src/cpu/gpw/gpw_evaluator.cpp` | Minimal-diff `=`→`+=`: `kEnergy` at `:399, 526, 572, 727, 874, 1010, 1145`; force-span blocks after each (`:435-437, 537-539, 577-579, 755-757, 934-936, 1025-1027, 1171-1173`); `kStress` `:1317`. `+= kEnergyElectronic` beside each `kEnergy`. **No branch restructuring** (smallest textual footprint against the CPMD sibling specs). |
| `src/cpu/planewave/pw_evaluator.cpp` | `:36, 44` energy, `:57-64` forces, `:37-39, 71-73` stress → `+=`; add key. |
| `src/cpu/semiempirical/am1_evaluator.cpp` | `:87` energy, `:103-105` forces → `+=` (+ key). **Charges stay `=`** (state channel — document at the site). |
| `src/cpu/semiempirical/rm1_evaluator.cpp` | Same as AM1 (`:78` energy + forces below). |
| `src/cpu/semiempirical/dftb2.cpp` | Already `+=`; add `+= kEnergyElectronic` at `:302`. |
| `include/atomiverse/core/frame_schema.h` | `inline constexpr MetaKey<double> kEnergyElectronic{"energy_electronic", true};` (after `kEnergyChargingWork`). |
| `src/cpu/frame/frame.cpp` | Pre-insert in `Frame(int)` (`:~484`). |
| `src/cpu/component/zero.cpp` | `ZeroEnergy` clears the new key (`:48-53` block); rewrite the invariant comment to the honest form: `kEnergy == Σ per-term + unclassified (Bias)` — Bias staying term-less is pre-existing, **named** debt. |

Single key rationale: co-bound ES sources legitimately merge into one
"electronic share of kEnergy" (Δ-learning, QM/MM's QM part). Rejected:
per-kind keys (`kEnergyGPW`, …) — key explosion with no consumer. Meta is
open-schema through molrs (`src/cpu/io/molrs_meta.h:55-58` iterates
`serializable_entries()`): the new serializable key round-trips with **zero
molrs coordination**.

## Files to create or modify

The table above, plus:

- NEW `tests/cpu/evaluate/test_source_cobinding.cpp` (U7)
- `tests/cpu/component/` ZeroEnergy test extension (U8)
- NEW `regressions/architecture/lj_rhf_coupled/{reg_lj_rhf_coupled.cpp,REFERENCE.md}` (R1)
- `regressions/CMakeLists.txt` — register R1 (**gated**)
- Audit pass over `tests/cpu/{gpw,planewave,semiempirical,scf}`: any evaluator
  unit test that launches twice on one frame without an intervening Zero gets
  an explicit zero or a fresh frame.

## Tasks

- [ ] RED (U7): Evaluate with a `+=` stub potential AND a nested RHF stub-able
  source — after run, `kEnergy == E_stub + E_RHF` and
  `kEnergyElectronic == E_RHF` (fails today: RHF overwrote the stub)
- [ ] RED (U8): `ZeroEnergy` clears `kEnergyElectronic`
- [ ] Register the meta key + pre-insert + ZeroEnergy clear + honest comment
- [ ] Flip MFM (accumulate flag), GPW (7+7+1 sites), PW, AM1, RM1; add the
  per-term writes; DFTB2 key only
- [ ] Grep sweep review (the anti-miss checklist): `grep -n "meta.get(kEnergy)
  *=" src/cpu/{gpw,planewave,semiempirical,scf}` and force-span `=` writes —
  every remaining `=` must be a standalone entry or a state channel, each
  justified in the PR text
- [ ] Test audit: double-launch evaluator unit tests fixed
- [ ] GREEN + full ctest + the AIMD/BOMD regression subset (sole-source
  equivalence proof)
- [ ] R1: `regressions/architecture/lj_rhf_coupled` — one `Evaluate`, two Ar
  atoms under `PairLJCut` (literature ε/σ) + nested RHF H₂ placed beyond the
  Ar cutoff; assert `kEnergy == 4ε[(σ/r)¹²−(σ/r)⁶] + E_RHF` to 1e-6 with
  `E_RHF(H₂/STO-3G) = −1.117506 Ha` (Szabo & Ostlund, Table 3.11) and
  `kEnergyElectronic == E_RHF` alone
- [ ] REFERENCE.md: cite Szabo & Ostlund and the LJ parameter source; state
  that the SUM is the pin of the `+=` contract
- [ ] Register R1 (gated)

## Testing strategy

**Unit:** U7 (co-binding sums; also route one co-bound config through a
GPW-shaped evaluator if buildable cheaply, else the grep checklist covers
GPW), U8 (zeroing of the key). RED first.

**Regression (R1):** both terms individually literature-pinned; the sum pins
the contract. ≤ 60 s, public API only.

**Risk gate:** a missed GPW site silently reverts to overwrite *and passes
every sole-source test* — the grep sweep task is mandatory, not advisory.

## Out of scope

- Bias per-term key (named debt, recorded in zero.cpp comment)
- Test-local overwrite stubs (`MullerBrown`, `ConstHydrostaticStress`) — valid
  sole sources, out of the engine contract
- QM/MM coupling internals (Wave 2) and CPMD branches' restructuring (Wave 3)

## Notes for implementers

- HARD GATE: before touching `src/cpu/scf/mean_field_method.cpp`,
  `src/cpu/gpw/gpw_evaluator.cpp`, or any CMakeLists, verify `git status`
  shows them clean. As of 2026-07-26 the operator's libcint-removal/C++20
  work (~760 files) is uncommitted and includes them. If dirty: STOP the spec
  and report — do not rebase over uncommitted operator work.
- Sole-source equivalence is bit-exact only because Zero* precedes at order
  −1000/−999/−998; a user binding a source at order < −999 breaks it — note
  in the contract docs, no code guard.
- Coordinate with Wave 3: the CPMD specs rebase over the GPW diff — keep it
  minimal-diff (`=` token → `+=` token), resist any tidying.
