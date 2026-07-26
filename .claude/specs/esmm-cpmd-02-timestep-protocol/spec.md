---
slug: esmm-cpmd-02-timestep-protocol
status: approved
created: 2026-07-26
revised: 2026-07-26
domain: electronic-propagation / cpmd
queue: esmm-campaign
deferred: false
depends_on: [esmm-cpmd-01-electronic-home]
unblocks: [esmm-cpmd-03, esmm-cpmd-04]
grilled: true
---

# Spec: `esmm-cpmd-02-timestep-protocol` — dt single-sourcing + cpmd_dt deletion

**Slug:** `esmm-cpmd-02-timestep-protocol`
**Domain:** MD timestep protocol (CPU)
**Status:** approved (grilled via plan-mode design review + operator Q&A, 2026-07-26)

## Summary

CPMD requires the electronic step to run on the SAME clock as the nuclear
integrator, but today `GPWEvaluator::cpmd_dt` duplicates `VelocityVerlet`'s
private `dt_` on faith ("must match the nuclear dt" is a docs-only contract —
no runtime check is possible). Kill the duplicated truth with an opt-in
protocol pair pushed by `cpu::md::MD` at build time, and **delete `cpmd_dt`**
(house rule: no back-compat). Also land the missing
`cpmd && compute_stress → throw` gate: today `compute_stress=true` is
**silently ignored** on every CPMD branch (they return before the stress
block at `gpw_evaluator.cpp:~1300`, and no build gate exists).

## Design

### New `include/atomiverse/cpu/component/timestep_protocol.h`

CUDA-free header beside `neighbor_protocol.h` / `long_range_protocol.h`:

```cpp
/// @brief Supply side: a component that owns the MD integration step.
class TimestepProvider {
  public:
    virtual ~TimestepProvider() = default;
    /// @brief Nuclear MD timestep [fs]; must be finite and > 0.
    virtual double timestep() const = 0;
};

/// @brief Demand side: a component/method that must know the nuclear dt.
class TimestepConsumer {
  public:
    virtual ~TimestepConsumer() = default;
    /// @brief True when this consumer cannot run without a pushed dt
    ///        (e.g. `cpmd == true`). Gates the missing-provider throw.
    virtual bool timestep_required() const = 0;
    /// @brief Receive the integrator timestep [fs] (before component build).
    /// @throws std::invalid_argument when dt_fs is not finite and > 0.
    virtual void set_timestep(double dt_fs) = 0;
};
```

`Component` gains the probe accessors beside `as_neighbor_consumer()` etc.:
`virtual const TimestepProvider* as_timestep_provider() const { return nullptr; }`
and `virtual TimestepConsumer* as_timestep_consumer() { return nullptr; }`.
Nested method drivers (`MeanFieldMethod`, a `DriverNode`) are reached by
`dynamic_cast<TimestepConsumer*>` in the push walk (the `HF` cross-cast
precedent, `md.h:384-387`).

### Providers (one-line accessor override each)

`VelocityVerlet`, `NoseHooverNPT`, `NphIntegrator`, `ParrinelloRahmanNPT`,
`NVELimit`, `Brownian` — each implements `double timestep() const` returning
its ctor dt [fs].

### Push site — end of `MD::install_defaults()` (`md.h:395`)

Runs before any component `build()`, so pushed dt is visible to build-time
validation. Walk components via `for_each_pending_mut`; walk nested driver
nodes via a new mutable kind-walk `for_each_pending_kind_mut` (beside
`driver.h:317`). Rules:

- ≥2 providers whose dt differ by >1e-12 **relative** → throw naming both.
- ≥1 consumer with `timestep_required()` and 0 providers → throw.
- Otherwise push dt to every consumer.

Evaluate/Minimize need **zero changes**: no push happens there, and the
consumer's own build throws — `GPWEvaluator::build`:
`cpmd && !(cpmd_dt_fs_ > 0)` → "cpmd requires an MD integrator publishing the
timestep (TimestepProvider); CPMD outside MD is meaningless". (The same gate
lands on `MeanFieldMethod::build_on` in esmm-cpmd-04.)

### `cpmd_dt` deletion — verified 16-file touch list

`gpw_evaluator.h` (member + docs; `GPWEvaluator` gains `TimestepConsumer`
with `timestep_required() ⇔ cpmd`, private `double cpmd_dt_fs_ = 0`),
`gpw_evaluator.cpp` (gate at `:194-197` deleted; three `advance(...,
cpmd_dt, ...)` call sites `:453,:553,:557,:590` → `cpmd_dt_fs_`),
`tests/cpu/gpw/test_gpw_cpmd.cpp` (negative-dt guard test replaced by
set_timestep-validation + missing-provider tests), 5 regression sources +
REFERENCE.md (`cpmd-water`, `cpmd-electronic-thermostat`, `cpmd-kpoints`,
`cpmd-uks`, `cpmd-rattle-water` — dt now lives only on `VelocityVerlet`),
docs `tutorial/methods/cpmd-{water,rattle-water,kpoints,electronic-thermostat}.md`,
`docs/roadmap/molecular-dynamics.md`.

### Missing stress gate (same file, same change)

`GPWEvaluator::build`: `cpmd && compute_stress` → typed throw ("CPMD stress
path is not productized; variable-cell CP is deferred") + guard test.
Today the user's stress request is silently dropped — the exact failure mode
the gate lattice exists to prevent.

## Files to create or modify

- NEW `include/atomiverse/cpu/component/timestep_protocol.h`
- NEW `tests/cpu/dynamics/test_timestep_protocol.cpp`
- `include/atomiverse/cpu/component/component.h` (two probe accessors)
- `include/atomiverse/cpu/driver/driver.h` (`for_each_pending_kind_mut`)
- `include/atomiverse/cpu/dynamics/{velocity_verlet.h,nose_hoover_npt.h,nph_integrator.h,parrinello_rahman_npt.h,nve_limit.h,brownian.h}`
- `include/atomiverse/cpu/dynamics/md.h` (push block)
- `include/atomiverse/cpu/gpw/gpw_evaluator.h`, `src/cpu/gpw/gpw_evaluator.cpp`
- `tests/cpu/gpw/test_gpw_cpmd.cpp`
- `regressions/dynamics/{cpmd-water,cpmd-electronic-thermostat,cpmd-kpoints,cpmd-uks,cpmd-rattle-water}/` (+ REFERENCE.md)
- `docs/tutorial/methods/cpmd-*.md` ×4, `docs/roadmap/molecular-dynamics.md`

## Tasks

- [ ] RED: `test_timestep_protocol.cpp` — conflicting providers throw naming
      both; required consumer + no provider throws; push received by a stub
      consumer; non-required consumer with no provider is fine
- [ ] Land the protocol header + Component accessors + `for_each_pending_kind_mut`
- [ ] Implement the 6 integrator providers + the MD push block
- [ ] GPWEvaluator: TimestepConsumer; delete `cpmd_dt`; migrate 3 advance sites;
      build gate for missing dt; build gate for `cpmd && compute_stress`
- [ ] Update `test_gpw_cpmd.cpp` guards (delete stale, add new — same change)
- [ ] Update the 5 cpmd regressions + REFERENCE.md + 4 tutorial pages + roadmap
- [ ] Scoped build + ctest + the 5 cpmd regressions green

## Testing strategy

**Unit:** `test_timestep_protocol.cpp` (above); `test_gpw_cpmd.cpp` gains
set_timestep validation (non-finite/≤0 throws), Evaluate+cpmd-evaluator
build-throw, stress-gate throw. **Regression:** the 5 existing cpmd cases
re-pass with dt sourced only from `VelocityVerlet` — their physics pins are
the no-change gate. No new regression (protocol is infrastructure; its
feature witness is the cpmd cases running on a single dt truth).

## Out of scope

- SCF-side CPMD members (esmm-cpmd-04 lands `MeanFieldMethod`'s consumer)
- Making thermostats (`NoseHoover`, `FixCSVR`, `FixLangevin`) providers —
  would cross-validate thermostat dt but changes throw behavior of shipping
  setups; wants a separate operator decision
- Variable-cell CP physics (only the throw gate lands here)

## Notes for implementers

- **HARD GATE:** before touching `src/cpu/gpw/gpw_evaluator.{h,cpp}`,
  `src/cpu/scf/mean_field_method.{h,cpp}`, `scf/internal/gradient.*`, or any
  CMakeLists, verify `git status` shows them clean. As of 2026-07-26 the
  operator's libcint-removal/C++20 work (~760 files) is uncommitted and
  includes them. If dirty: STOP and report.
- Grep-verify deletion: `grep -rn "cpmd_dt" include/ src/ tests/ regressions/ docs/`
  must return zero hits after this spec.
- Guard-test hygiene (house rule): every gate added/removed updates its
  EXPECT_THROW twin in the same change.
