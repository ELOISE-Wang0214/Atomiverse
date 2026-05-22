---
title: CUDA driver taxonomy — 02 MD reparent and component tagging
status: approved
created: 2026-05-22
follows: cuda-driver-taxonomy-01-foundations
---

# CUDA driver taxonomy — 02 MD reparent and component tagging

## Summary

Reparents the CUDA molecular-dynamics driver `cuda::md::MD` from the raw
`Driver<md::Phase, cuda::device_tag>` template onto the category base
`cuda::DynamicsDriver` delivered by link 01-foundations, so the CUDA MD
driver inherits exactly the same outer-loop semantics, `outer_status`
handling, and build-time component validation that `cpu::md::MD` already
has. As part of the reparent it fixes a latent design bug: the
`ZeroForces` / `ZeroEnergy` / `ZeroKinetic` infrastructure components
currently hardcode `md::Phase` bindings, which `cuda::minimizer::CG`/`SD`
work around by reaching into `md::ZeroForces` and re-binding it. After
this spec the Zero* components are phase-agnostic and each concrete
driver binds them to its own phase enum via explicit `on(phase).add()`,
unblocking 03-minimizer. Every CUDA MD-domain component gains explicit
`role()` / `capabilities()` overrides matching the CPU MD components, so
`DynamicsDriver::validate_components()` accepts a well-formed MD setup
(a `VelocityVerlet` integrator plus a `Force`-capable evaluator) and
throws on a malformed one. No physics changes — the existing CUDA MD
regression suite is the behaviour-preservation gate.

## Domain basis

Pure structural refactor — no equations introduced or modified.
Physical correctness is preserved by construction: the integrator,
thermostat and force-kernel code paths are untouched; only the driver
inheritance, component phase binding, and `role()`/`capabilities()`
tags change. Correctness is verified by the existing CUDA MD regression
tests — NVE energy conservation (`test_cuda_md`) and Nosé-Hoover
target-temperature convergence (`test_cuda_md_thermostat`) — which must
pass unchanged.

## Design

**Reparent.** At HEAD-of-01 `cuda::md::MD` reads
`class MD : public Driver<Phase>` (the interim 2-arg form). This spec
changes it to `class MD : public cuda::DynamicsDriver`. After the
reparent, `MD` inherits the category base's continuable `run()` (full
step budget, honours `outer_status.stop_requested`, ignores
`outer_status.converged`) and `validate_components()`. `MD` itself no
longer names `Driver<...>` directly anywhere — neither in its base-class
clause nor in member signatures.

**Phase-agnostic Zero* components.** Today `ZeroForces` / `ZeroEnergy` /
`ZeroKinetic` (`cuda/md/init.cuh`) are bound by `MD`'s constructor via
`on(Phase::pre_force).add(ZeroForces{}, -1000)` etc. — but the
`cuda::minimizer::CG`/`SD` drivers also depend on `md::ZeroForces`, which
forces a `md`-namespace dependency from the minimizer domain and an
explicit re-bind at `minimizer::Phase::force`. Mirroring the CPU
resolution (`cpu/md/init.h`), the Zero* classes stay phase-agnostic:
they expose no `default_bindings()` pointing at a fixed `md::Phase`, and
each concrete driver adds them with explicit `on(phase).add()` for ITS
own phase enum. `cuda::md::MD` keeps binding them at `md::Phase` in its
constructor; 03-minimizer will bind the same classes at
`minimizer::Phase`. No `default_bindings()` override means
`add()`-without-`on()` is rejected for these components — that is
intentional, exactly as on CPU.

**Component tagging.** The CUDA force-evaluator bases (`cuda::Pair`,
`cuda::Bond`, `cuda::Angle`, `cuda::Dihedral`) already declare
`role()==Evaluator` and `capabilities()==Energy|Force|Stress`, so
`cuda::PairLJCut` and the bonded evaluators are already correctly
tagged through their bases — no change needed there. The work in this
spec is the MD-domain components that currently inherit the
`ComponentBase` defaults (`role()==None`, `capabilities()==None`):

| Component | `role()` | `capabilities()` |
|---|---|---|
| `cuda::md::VelocityVerlet` | `Integrator` | `None` |
| `cuda::md::NoseHoover` | `Thermostat` | `None` |
| `cuda::md::Barostat` | `Barostat` | `None` |
| `cuda::md::Deform` | `Integrator` | `None` |
| `cuda::md::InitVelocity` | `None` (structural) | `None` |
| `cuda::md::KineticEnergy` | `None` (structural) | `None` |
| `cuda::md::ZeroForces`/`ZeroEnergy`/`ZeroKinetic` | `None` (structural) | `None` |

These match the CPU MD components verbatim (`cpu/md/velocity_verlet.h`,
`nose_hoover.h`, `barostat.h`, `deform.h`, `init.h`). The CPU `Barostat`
also overrides `reads()==Stress`; CUDA `ComponentBase` does not yet
expose a `reads()` virtual, so the CUDA `Barostat` declares only
`role()==Barostat` here — the `reads()` mirror is out of scope and left
to 05-scrub-and-docs or a later CUDA auto-injection spec. `KineticEnergy`
writes `frame.kinetic`, which is not a `Capability` bit, so it stays
`capabilities()==None` (matching CPU). The Zero* components keep
`capabilities()==None`: zeroing an accumulator is not "producing" the
quantity, and `DynamicsDriver` validation must not mistake a `ZeroForces`
for the `Force`-capable evaluator the integrator needs.

**Validation wiring.** `cuda::DynamicsDriver::validate_components()`
(from 01) mirrors CPU exactly: it requires ≥1 component with
`role()==Integrator` and ≥1 component with `Capability::Force`. After
tagging, a `VelocityVerlet` + a `PairLJCut` satisfy both; an MD setup
with no integrator throws at `build()`.

## Files to create or modify

- `include/atomiverse/cuda/md/md.cuh` — reparent base class to
  `cuda::DynamicsDriver`; update the constructor's explicit
  `on(phase).add()` calls for the phase-agnostic Zero* components;
  refresh the class docstring.
- `include/atomiverse/cuda/md/velocity_verlet.cuh` — add
  `role() == ComponentRole::Integrator` override.
- `include/atomiverse/cuda/md/nose_hoover.cuh` — add
  `role() == ComponentRole::Thermostat` override.
- `include/atomiverse/cuda/md/barostat.cuh` — add
  `role() == ComponentRole::Barostat` override.
- `include/atomiverse/cuda/md/deform.cuh` — add
  `role() == ComponentRole::Integrator` override.
- `include/atomiverse/cuda/md/init.cuh` — drop any hardcoded
  `default_bindings()` from `ZeroForces`/`ZeroEnergy`/`ZeroKinetic`,
  make them phase-agnostic; confirm `InitVelocity` role stays `None`;
  refresh docstrings.
- `include/atomiverse/cuda/md/kinetic_energy.cuh` — confirm `role()`
  stays `None`; docstring note that it writes a non-`Capability` field.
- `src/cuda/md/velocity_verlet.cu` — adjust if the `role()` override is
  defined out-of-line (header-inline preferred, matching CPU).
- `src/cuda/md/nose_hoover.cu` — same as above for `NoseHoover`.
- `src/cuda/md/barostat.cu` — same as above for `Barostat`.
- `src/cuda/md/deform.cu` — same as above for `Deform`.
- `src/cuda/md/init.cu` — remove `default_bindings()` definitions for
  the Zero* components if they live in the `.cu`.
- `src/cuda/md/kinetic_energy.cu` — no functional change expected;
  touch only if a `role()` definition lives here.
- `tests/cuda/md/test_md.cu` — add the RED structural + validation
  tests (is-a `DynamicsDriver`, throw-on-missing-integrator,
  Zero* phase-agnostic).
- `tests/cuda/md/test_thermostat.cu` — confirm it still compiles and
  passes against the reparented `MD`.

## Tasks

- [ ] Write failing tests for the MD reparent and tagging in `tests/cuda/md/test_md.cu` — assert `std::is_base_of_v<cuda::DynamicsDriver, cuda::md::MD>`, that building an `MD` with no `Integrator` component throws, and that `ZeroForces`/`ZeroEnergy`/`ZeroKinetic` return an empty `default_bindings()`
- [ ] Make `ZeroForces`/`ZeroEnergy`/`ZeroKinetic` phase-agnostic in `include/atomiverse/cuda/md/init.cuh` (and `src/cuda/md/init.cu`) by removing hardcoded `md::Phase` `default_bindings()`
- [ ] Reparent `cuda::md::MD` to `cuda::DynamicsDriver` in `include/atomiverse/cuda/md/md.cuh` and bind the phase-agnostic Zero* components via explicit `on(md::Phase::...).add()` in the constructor
- [ ] Add `role()` overrides to `cuda::md::VelocityVerlet`, `NoseHoover`, `Barostat`, `Deform` (`velocity_verlet.cuh`, `nose_hoover.cuh`, `barostat.cuh`, `deform.cuh`) matching the CPU MD tags
- [ ] Confirm `InitVelocity` (`init.cuh`) and `KineticEnergy` (`kinetic_energy.cuh`) keep `role()==None`/`capabilities()==None` and document why
- [ ] Verify `cuda::md::MD` builds through `cuda::DynamicsDriver::validate_components()` with a `VelocityVerlet` + `PairLJCut` and that `test_cuda_md_thermostat` still compiles against the reparented driver
- [ ] Run full check + test suite (`ctest --test-dir build -R cuda_md --output-on-failure`)

## Testing strategy

- **Happy path** — `cuda::md::MD` with `VelocityVerlet` + `PairLJCut` +
  `NeighborList` builds clean through `DynamicsDriver::validate_components()`
  and runs the requested step budget.
- **Edge case — missing integrator** — building an `MD` containing a
  `PairLJCut` (Force-capable) but no `Integrator`-role component throws
  at `build()` from `validate_components()`.
- **Edge case — missing force writer** — building an `MD` with a
  `VelocityVerlet` but no `Force`-capable evaluator throws at `build()`.
- **Structural** — compile-time `std::is_base_of_v<cuda::DynamicsDriver,
  cuda::md::MD>` holds; `cuda::md::MD` does not name `Driver<...>`
  directly; Zero* components return an empty `default_bindings()`.
- **Behaviour preservation (domain validation)** — the existing CUDA MD
  regression suite passes unchanged: `test_cuda_md` (NVE energy
  conservation / drift bound) and `test_cuda_md_thermostat`
  (Nosé-Hoover target-temperature convergence). Run via
  `ctest --test-dir build -R cuda_md --output-on-failure`. Requires a
  CUDA toolkit and a CUDA-capable device.

## Out of scope

- CPU-style categorical-setter API on `cuda::DynamicsDriver` /
  `cuda::md::MD` (`set_pair` / `set_integrator` / `set_thermostat` / …) —
  chain decision: the CUDA category bases provide only the
  loop / validation / `outer_status` core.
- Mirroring CPU `Barostat::reads()==Stress` and the corresponding
  `VirialStress` auto-injection — CUDA `ComponentBase` has no `reads()`
  virtual; deferred.
- `cuda::minimizer` (CG/SD) re-binding the now phase-agnostic Zero*
  components at `minimizer::Phase` — that is link 03-minimizer.
- `cuda::scf::SCF` reparent — link 04-scf.
- Removal of the legacy `ATV::md::*` namespace aliases / final
  docs scrub — link 05-scrub-and-docs.
- Any change to integrator / thermostat / force-kernel numerics.
