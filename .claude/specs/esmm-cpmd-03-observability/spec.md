---
slug: esmm-cpmd-03-observability
status: code-complete
created: 2026-07-26
revised: 2026-07-26
domain: electronic-propagation / cpmd
queue: esmm-campaign
deferred: false
depends_on: [esmm-cpmd-02-timestep-protocol]
unblocks: [esmm-cpmd-04, esmm-cpmd-06]
grilled: true
---

# Spec: `esmm-cpmd-03-observability` — conserved quantity + adiabaticity on Frame

**Slug:** `esmm-cpmd-03-observability`
**Domain:** CPMD observability (CPU)
**Status:** approved (grilled via plan-mode design review + operator Q&A, 2026-07-26)

## Summary

CPMD's health signals — fictitious electronic kinetic energy T_e (the
adiabaticity monitor) and the Car–Parrinello conserved energy
E_cons = E_KS + T_nuclear + T_e + E_thermostat — exist today only as
frame-invisible evaluator readbacks (`last_electronic_ke_ev()`,
`last_cp_energy_ev()`). No Observer/Recorder can watch CP drift; a Thermo
line cannot show whether `cpmd_mu`/dt are safe for the system at hand. Put
them on Frame meta and assemble E_cons into the existing `kConserved` column
Thermo already prints.

## Design

### Meta keys — one-liners in `core/frame_schema.h`

```cpp
/// Fictitious CP electronic kinetic energy T_e [eV] (CPMD evaluators; adiabaticity monitor).
inline constexpr MetaKey<double> kEnergyElectronKinetic{"energy_electron_kinetic", true};
/// Nose–Hoover electronic thermostat energy [eV] (0 for NVE electrons).
inline constexpr MetaKey<double> kEnergyElectronThermostat{"energy_electron_thermostat", true};
```

Serializable (meta is open-schema through molrs — zero coordination). **Not**
pre-inserted by `Frame(int)`: producers `meta.insert()` with a `contains()`
guard (the `kDipole` pattern, `observe/multipole.cpp:97-98`). Every GPW CPMD
branch (Γ RKS/UKS, multi-k) writes both keys on every launch (0 on the seed
step). Readbacks stay.

### New component `cpu::md::CpmdConserved` (`cpu/dynamics/cpmd_conserved.h`)

```cpp
/// @brief CpmdConserved — end-of-step observer assembling the Car–Parrinello
///        conserved energy into frame.meta[kConserved] [eV]:
///        E_cons = kEnergy + kKinetic + kEnergyElectronKinetic
///               + kEnergyElectronThermostat.
///        Default binding md::Phase::end_of_step @ -90 (after KineticEnergy
///        @ -100, before Reduce @ -50 / Publish / Thermo) so the Thermo
///        "Conserved" column shows CP drift and adiabaticity transfer.
///        Throws at launch when the electronic keys are absent (no
///        CPMD-enabled evaluator ran — fail fast, no fallback).
///        Do not combine with integrators that own kConserved (NPT / CSVR).
class CpmdConserved : public cpu::Component { /* name(), default_bindings(), launch() */ };
```

Rationale (decided, do not relitigate): the evaluator cannot write E_cons
itself — at force phase `kKinetic` is half-step stale (KineticEnergy runs at
end_of_step @ -100). `kConserved` already has extended-system owners
(`nph_integrator.cpp:295`, `nose_hoover_npt.cpp:886`,
`parrinello_rahman_npt.cpp:622`, `fix_csvr.cpp:161`) — therefore
**explicit wiring, not auto-injection**, with a documented do-not-combine
(CPMD+NPT is not a product path). Drift termination needs no library
criterion: `Criterion::verdict(LoopState)` sees the frame, so a user-space
`StopWhenDrift` reading `kConserved` is ~15 lines — the tutorial page gets
the snippet.

### Regression rewire

`regressions/dynamics/cpmd-water` gains `CpmdConserved` and pins the relative
`kConserved` drift envelope directly off the meta key (replacing the
readback-arithmetic assembly it does today). Existing tolerance unchanged.

## Files to create or modify

- `include/atomiverse/core/frame_schema.h` (two meta keys)
- NEW `include/atomiverse/cpu/dynamics/cpmd_conserved.h`,
  `src/cpu/dynamics/cpmd_conserved.cpp`
- `src/cpu/gpw/gpw_evaluator.cpp` (write both keys on every CPMD launch —
  all Γ RKS/UKS + multi-k branches)
- NEW `tests/cpu/dynamics/test_cpmd_conserved.cpp`
- `regressions/dynamics/cpmd-water/` (+ REFERENCE.md note)
- `docs/tutorial/methods/cpmd-water.md` (drift-criterion user-space snippet)
- `CMakeLists.txt` (one new source)

## Tasks

- [x] RED: `test_cpmd_conserved.cpp` — launch with pre-seeded meta values
      asserts the sum lands on kConserved; launch with the electronic keys
      absent throws; default binding is end_of_step @ -90
      (target test_cpu_md_cpmd_conserved; 5 tests incl. thermostat-key-alone
      hole; name pinned "cpmd_conserved"; throw messages must name key+cpmd)
- [x] Register the two meta keys; GPW CPMD branches insert-guarded write both
      each launch (0 on seed step)
      (5 publication sites via shared helper; bit-for-bit with readbacks by
      construction; also fixed stale last_etherm_ev_ on re-build seed — named)
- [x] Implement `CpmdConserved` (component + default_bindings + launch)
- [x] Rewire `cpmd-water` regression onto `CpmdConserved` + meta-key drift pin
      (pin block byte-identical; binary output identical pre/post; producer
      test H2 fixture ULP-verified; REFERENCE μ=400→2000 doc rot fixed)
- [x] Tutorial snippet: user-space `StopWhenDrift` criterion
      (compile-verified -fsyntax-only; stale hand-assembly route replaced;
      sibling regressions' readback-assembly NAMED as follow-up)
- [x] Scoped build + ctest + cpmd regressions green
      (full ctest 377/377 — incl. the previously-transient uhf_hessian —
      + 5/5 cpmd regressions PASS)

## Testing strategy

**Unit:** `test_cpmd_conserved.cpp` (sum, missing-key throw, binding).
A `test_gpw_cpmd.cpp` case asserts both meta keys exist and are finite after
a CPMD launch, and that T_e matches `last_electronic_ke_ev()` bit-for-bit.
**Regression:** `cpmd-water` (existing case, upgraded assertion source) — no
new regression; the feature's witness is the existing drift pin now flowing
through Frame meta.

## Out of scope

- SCF-side key writes (esmm-cpmd-04 reuses the same keys)
- A provides/consumes capability bit for kConserved multi-writer detection
  (documented-not-gated; heavier design, operator may revisit)
- Thermo column additions (Conserved column already prints)

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/gpw/gpw_evaluator.{h,cpp}`,
  `src/cpu/scf/mean_field_method.{h,cpp}`, `scf/internal/gradient.*`, or any
  CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and
  includes them. If dirty: STOP and report.
- Meta write pattern: `if (!f.meta.contains(kEnergyElectronKinetic)) f.meta.insert(...)` then `set` — never assume pre-insertion.
- eV conversions single-sourced through `units::` (T_e is Hartree-native in
  the propagator).

> Hygiene (/mol:simplify 2026-07-27): dead public accessor last_cp_energy_ev
> deleted (0 call sites; taught the retired hand-assembly); stale CP comment
> updated to the meta channel; 2 format drifts fixed; header wiring example
> corrected to add() (on().add defaults to order 0, after Reduce). NAMED:
> GPWEvaluator::launch 1027-line sprawl (pre-existing, route /mol:refactor);
> 2 rule captures for /mol:note (toolchain gap 4th re-arm; open-schema meta
> producer pattern now in 6 files).
