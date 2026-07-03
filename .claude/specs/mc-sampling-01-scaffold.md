---
title: MC sampling scaffold — sampling::Phase, Mover role, Transaction, SamplingDriver
status: in-progress
created: 2026-07-03
---

> **Impl note (2026-07-03):** all scaffold code is written and compiles cleanly
> — the test TU `tests/cpu/test_mc_scaffold.cpp` front-end-compiles to an object
> with the project's real flags (`-DATV_BUILD_TESTS`, xtensor/gtest includes),
> exercising every new symbol incl. the `CPUDriver<sampling::Phase>` template
> instantiation. The full build + `ctest` run is **blocked by a pre-existing
> `dev`↔molrs bridge gap unrelated to this change**: `src/cpu/recorder/xyz.cpp`
> and `src/cpu/observer/trajectory_metrics.cpp` call molrs symbols
> (`read_frame_zarr_first`, `Impl::obj`) that exist in no molrs ref, so
> `libatomiverse` will not link. Code-type criteria are verified by
> inspection + isolated compile; runtime criteria stay pending until the molrs
> bridge is available.

# MC sampling scaffold — sampling::Phase, Mover role, Transaction, SamplingDriver

## Summary

This is sub-spec **01 of 3** in the "v1 canonical NVT Metropolis Monte Carlo sampler" chain; **02** adds the concrete `cpu::mc::MC` driver plus `TranslateMove`, and **03** adds physics validation (acceptance-ratio and Boltzmann-distribution checks). Sub-spec 01 lands only the **dormant machinery**: the `sampling::Phase` execution skeleton, a `ComponentRole::Mover` structural tag, a pure-data `cpu::mc::Transaction` changeset, two new non-owning `cpu::Context` handles, and the `cpu::SamplingDriver` category base with an inline-Metropolis `run()` shape and a `validate_components()` gate. Nothing is user-runnable yet — there is no concrete driver and no move — but every symbol compiles, is Doxygen-documented, and is covered by fast reference-free unit tests, so the fifth outer driver category (Sampling, per CLAUDE.md's "Adding a fifth user-facing category") is scaffolded green ahead of the move implementation in 02.

## Domain basis

The sampler targets the canonical (NVT) ensemble by the **Metropolis** criterion (Metropolis et al. 1953, DOI 10.1063/1.1699114; Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., ch. 3). For a proposed trial state a move is accepted with probability

```
acc = min(1, exp(-ΔE / (k_B T))),   ΔE = E_new - E_old
```

Any move with `ΔE ≤ 0` is accepted unconditionally (the `exp` argument is ≥ 0 so `min(1, ·) = 1`). Because v1 uses only **symmetric** proposals (equal forward/reverse proposal probability), the proposal terms cancel in the acceptance ratio and detailed balance reduces to this **energy-only** form; asymmetric-proposal generalizations are out of scope for the whole chain. Thermal energy is `k_B·T [eV] = kB_ev_per_k · T[K]` with `kB_ev_per_k = 8.617333262e-5 eV/K` from `core/units.h`.

**Bootstrap invariant (correctness):** the driver's running current energy `e_cur_` MUST be seeded from **one full-Frame energy evaluation of the initial configuration before the first move**. If `e_cur_` is left at zero (or any stale value), the first `ΔE` is wrong and detailed balance is violated for the opening move; every subsequent `e_cur_` inherits the error. This spec fixes the seam and the invariant; the numerical proof lives in sub-spec 03.

## Design

Five coordinated changes, all dormant until 02 supplies a `Mover` and a concrete driver:

- **`sampling::Phase` (in `facade/stage.h`).** New `namespace ATV::sampling` holding `enum class Phase { propose, energy, end_of_step, COUNT }` plus a `phase_names(Phase)` overload returning the three names in enum order — matching the exact shape of the sibling `md::Phase` / `vib::Phase` enums (every phase always walked, `COUNT` sentinel, per-value doc comments). There is **deliberately no `accept` phase**: accept/reject is not a per-phase component slot but the driver's outer-loop control inside `run()` (documented on `SamplingDriver`), exactly as MD's step-budget termination is not a phase.

- **`ComponentRole::Mover` (in `facade/component_role.h`).** New structural-axis role inserted after `Barostat` (keeping the Integrator/Optimizer/Thermostat/Barostat cluster contiguous), doc-commented "Proposes a trial state change for a sampling driver (TranslateMove, ...)." Orthogonal to `capabilities()` as with every other role.

- **`cpu::mc::Transaction` (new `cpu/mc/transaction.h`, `namespace ATV::cpu::mc`).** A standalone public **pure-data** struct: the per-move changeset — moved-atom indices (`LocalInt`), and per-moved-atom old and new `(x,y,z)` positions (`StorageReal`, the Frame coordinate scalar type). It carries **no** methods that reach into neighbor lists or the Frame — consumers derive their own views (CLAUDE.md: "hot-path data structures carry only fundamental physical quantities"). It exists to **apply** a trial move and to **roll it back** on rejection. A validating constructor (checking the index/position array lengths agree) is provided; there is **no** factory function. It is **never** a member of `cpu::Frame` (Frame-purity) — it lives on the driver and is published through `Context`.

- **`cpu::Context` wiring (in `cpu/component/context.h`).** Two new non-owning handles: `mc::Transaction* proposal = nullptr;` and `std::mt19937_64* rng = nullptr;`, with `mc::Transaction` forward-declared exactly as `md::PairList` is at lines 16-17. Documented like `pair_list`: published by the sampling driver during setup, ignored by non-sampling drivers.

- **`cpu::SamplingDriver` (new `cpu/driver/sampling_driver.h`).** Category base `: public CPUDriver<sampling::Phase>`, owning `mc::Transaction proposal_`, `std::mt19937_64 rng_`, and `AccReal e_cur_`. Its `run(cpu::Frame&, int n_sweeps)` clones `DynamicsDriver`'s fixed-budget shape (setup once; per-sweep `mutable_outer_status().reset()` → `step_once(idx)` → break only on `outer_status().stop_requested`; `converged` ignored) with the Sampling additions: (a) **before** the loop, seed `e_cur_` from one initial full-Frame energy evaluation (bootstrap invariant); (b) each sweep, after `step_once`, read `E_trial = frame.scalars.energy`, apply inline Metropolis `min(1, exp(-ΔE/(k_B T)))` against `e_cur_` using the owned `rng_`, then either **commit** (`e_cur_ = E_trial`) or **roll back** the moved-atom positions from `proposal_.old_*`. A `// ponytail:` seam comment marks where a swappable acceptance policy (REMD/GCMC) would later plug in. `validate_components()` clones `DynamicsDriver`'s (dynamics_driver.h:115): require ≥ 1 component with `role()==ComponentRole::Mover` **and** ≥ 1 with `Capability::Energy`, else `throw_missing(...)`. The continuable-run bookkeeping and the `molrs::FrameRef` overload are **stubbed/omitted** in 01 (no concrete driver exists to exercise them) — kept minimal but compiling. Members `proposal_`/`rng_`/`e_cur_` live on the base; publishing `&proposal_`/`&rng_` onto `Context` happens in the setup override.

  *Alternative considered — rejected for v1:* fold accept/reject inside the `Mover` component so the driver is a dumb loop and nothing lands on `Context` (the librarian-preferred shape). Rejected here in favor of the user's explicit inline-in-driver directive; revisit when a swappable acceptance policy (REMD/GCMC) makes the driver seam worthwhile.

## Files to create or modify

- `include/atomiverse/facade/stage.h` — add `namespace ATV::sampling` with `Phase` + `phase_names`.
- `include/atomiverse/facade/component_role.h` — add `ComponentRole::Mover`.
- `include/atomiverse/cpu/mc/transaction.h` (new) — `cpu::mc::Transaction` pure-data struct.
- `include/atomiverse/cpu/component/context.h` — forward-declare `mc::Transaction`; add `proposal` + `rng` fields.
- `include/atomiverse/cpu/driver/sampling_driver.h` (new) — `cpu::SamplingDriver` base.
- `tests/cpu/test_mc_scaffold.cpp` (new) — reference-free unit tests + test-only stub `Mover` / stub Energy evaluator.
- `tests/CMakeLists.txt` — register `test_mc_scaffold` via `atv_add_cpp_test`.

## Tasks

- [x] Write failing tests for the MC scaffold in `tests/cpu/test_mc_scaffold.cpp` (enum shape, Mover role, Transaction round-trip, Context pointers, `validate_components` throws twice) and register it in `tests/CMakeLists.txt`
- [x] Add `ATV::sampling::Phase` + `phase_names(Phase)` in `include/atomiverse/facade/stage.h` with per-value Doxygen and the no-`accept`-phase note
- [x] Add `ComponentRole::Mover` after `Barostat` in `include/atomiverse/facade/component_role.h` with its doc comment
- [x] Implement the pure-data `cpu::mc::Transaction` (validating constructor, no factory) in `include/atomiverse/cpu/mc/transaction.h` with Doxygen `@brief`/`@param` and units
- [x] Wire `mc::Transaction* proposal` and `std::mt19937_64* rng` (forward-declared `mc::Transaction`) into `cpu::Context` in `include/atomiverse/cpu/component/context.h`
- [x] Implement `cpu::SamplingDriver` in `include/atomiverse/cpu/driver/sampling_driver.h` — bootstrap-seeded fixed-budget `run()`, inline Metropolis commit/rollback, `// ponytail:` seam, and `validate_components()` requiring Mover + Capability::Energy
- [x] Add Doxygen `@brief`/`@param`/`@return` (with units) to every new public symbol and confirm the no-factory / Frame-purity constraints hold
- [ ] Run full check + test suite — **BLOCKED** by the pre-existing `dev`↔molrs bridge gap (see Impl note); isolated TU compile is green

## Testing strategy

All tests are fast, reference-free CPU unit tests (`atv_add_cpp_test`), no physics:

- **Happy path — enum shape:** `sampling::Phase` has values `propose, energy, end_of_step, COUNT` in that order; `static_cast<int>(Phase::COUNT) == 3`; `phase_names(Phase{}).size() == 3` and equals `{"propose","energy","end_of_step"}`.
- **Happy path — role tag:** `ComponentRole::Mover` exists and is distinct from `Barostat`/`Integrator` (compile + value check).
- **Happy path — Transaction round-trip:** construct a `Transaction` with N moved indices and matching old/new position arrays; read back indices, old `(x,y,z)`, new `(x,y,z)` unchanged.
- **Edge — Transaction validation:** the validating constructor rejects mismatched index/position array lengths (throws).
- **Happy path — Context defaults:** a default-constructed `cpu::Context` has `proposal == nullptr` and `rng == nullptr`.
- **Edge — validate_components gate:** `SamplingDriver::validate_components()` throws when no `Mover` is registered, and throws when a `Mover` is present but no component advertises `Capability::Energy` — exercised with a test-only stub `Mover` (role `Mover`, capabilities `None`) and a stub Energy evaluator.
- **Build gate:** the dormant machinery compiles green under `-DATV_BUILD_TESTS=ON` and the full `ctest` gate passes (implicit in the `runtime: full check + test suite`).

Explicitly **not** tested in 01: accept/reject dynamics, rollback correctness under a real move, bootstrap-energy numerics — no move exists yet; those land in 02/03.

## Out of scope

- The concrete `cpu::mc::MC` driver and `TranslateMove` component (→ 02).
- Any accept/reject dynamics, rollback-under-real-move, or bootstrap-energy test (→ 02/03).
- All physics validation: acceptance-ratio targets, Boltzmann-distribution / detailed-balance checks (→ 03).
- Local (per-move) ΔE evaluation — v1 reads the full-Frame `scalars.energy`.
- Other move types (volume, rotation, insertion/deletion), asymmetric proposals.
- CUDA backend (`cuda::` sampling).
- A swappable acceptance policy object (REMD/GCMC) — only the `// ponytail:` seam is placed now.
