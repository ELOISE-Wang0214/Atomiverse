---
title: MC Sampling 02 — Concrete MC driver + TranslateMove Mover
status: approved
created: 2026-07-03
---

# MC Sampling 02 — Concrete MC driver + TranslateMove Mover

## Summary
Deliver the first usable Metropolis Monte Carlo sampler on the CPU backend: a concrete `cpu::mc::MC` driver plus a single `TranslateMove` Mover component. After this sub-spec a user can assemble a working NVT MC run — `cpu::mc::MC mc{T, seed}; mc.set_pair(PairLJCut{...}); mc.set_neighbor(...); mc.set_mover(TranslateMove{delta}); mc.set_thermo(100); mc.build(); mc.run(frame, n_sweeps)` — reusing the existing `PairLJCut` energy evaluator with no new evaluator and no factory functions. The driver owns the RNG stream, the proposal transaction, and the cached accepted energy `e_cur_`, and drives the propose→energy→Metropolis→commit/rollback loop provided by the `SamplingDriver` base from sub-spec 01. This sub-spec proves core correctness (bootstrap, exact rollback, energy caching, Metropolis exactness, determinism) with fast, reference-free unit tests; the physics distribution gates (NIST LJ) land in sub-spec 03.

## Domain basis
Single-particle Metropolis Monte Carlo in the canonical (NVT) ensemble. Each trial proposes a symmetric displacement of one atom and accepts it with probability

  acc(o→n) = min(1, exp(−ΔE / k_B T)),  ΔE = E_trial − E_cur

Because the translation proposal is symmetric (forward and reverse proposal densities are equal for a uniform box `[−δ, δ]³`), the acceptance reduces to the energy-only Metropolis criterion and satisfies detailed balance with the Boltzmann distribution ∝ exp(−E/k_B T) (Metropolis et al. 1953, DOI 10.1063/1.1699114; Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed., ch. 3). `k_B` is taken from `core/units.h` (eV/K); temperature `T` in K, energies in eV — consistent with `frame.scalars.energy`.

Full-Frame ΔE decision (design rationale, not an optimization): `E_trial` is a full recompute over the mutated frame via the registered `Capability::Energy` evaluator, and `ΔE = E_trial − e_cur_` where `e_cur_` is cached from the last accepted state. Computing ΔE from a full recompute — rather than a local sub-Frame delta — is correct for *all* potentials including long-range / reciprocal-space (`kspace`) terms, where a local single-atom energy is undefined. That generality is the reason v1 recomputes the whole frame. On rejection the move is undone in place from `Transaction.old_pos`; the whole Frame is never copied.

## Design
- `cpu::mc::MC` (namespace `ATV::cpu::mc`) subclasses `SamplingDriver` (added in 01) → `CPUDriver<sampling::Phase>`. It mirrors the structure of `cpu::md::MD`:
  - Constructor binds always-on infra: `cpu::ZeroEnergy{}` at the **head** of `sampling::Phase::energy` at order **−999**, so the energy evaluator accumulates onto a freshly zeroed `frame.scalars.energy` each trial (exactly as MD binds `ZeroEnergy` at `pre_force −999`).
  - Local categorical setter `set_mover(...)` registers the Mover into `sampling::Phase::propose` via the inherited `set_singular_default_` helper (singular slot, idempotent on same-name re-call, throws on a different Mover). `set_pair` / `set_neighbor` / `set_observer` / `set_record` / `set_thermo` are inherited unchanged from `CPUDriver`.
  - `MCConfig` (or the driver's own validating constructor `MC{T, seed}`) carries temperature `T` [K] and RNG seed; both validated in the constructor (`T > 0`), **no factory function**.
  - Owns `std::mt19937_64 rng_` (seeded from `seed`), `mc::Transaction proposal_`, `AccReal e_cur_`, and publishes `&proposal_` and `&rng_` onto `cpu::Context` at setup (the `Context::proposal` / `Context::rng` handles added by 01).
  - `run(frame, n_sweeps)` uses the `SamplingDriver` base loop: bootstrap `e_cur_` from an initial full energy evaluation, then per sweep propose→energy→inline Metropolis (using `rng_`, `T`, `ΔE`)→commit (keep, update `e_cur_`) or rollback (restore from `proposal_.old_pos`).
  - `validate_components()` is provided by `SamplingDriver` (01): requires exactly one `role()==Mover` and ≥1 `capabilities()` containing `Capability::Energy`.
- `cpu::mc::TranslateMove` Mover component: `role()==ComponentRole::Mover`, `capabilities()==Capability::None` (it mutates positions; it is not an Energy/Force writer). Config: displacement half-width `δ` [Å] via a validating constructor (`δ > 0`), **no factory**. Bound to `sampling::Phase::propose`. On `launch(ctx)`: draw a uniform atom index `i ∈ [0, N)` from `ctx.runtime.rng`, draw a per-axis uniform displacement in `[−δ, δ]` for x, y, z, mutate `ctx.frame` in place, and record the changeset (moved index; old x,y,z; new x,y,z) into `*ctx.runtime.proposal`, clearing/reusing the `Transaction` each call. Single-atom move for v1. All randomness flows through the shared `ctx.runtime.rng` stream for reproducibility.
- Lifecycle: RNG/proposal/`e_cur_` are driver members published to `Context` at `setup`; no per-launch allocation. The Transaction is reused (cleared then filled) each proposal.

## Files to create or modify
- `include/atomiverse/cpu/mc/mc.h` (new) — `cpu::mc::MC` driver (header-only inline, mirroring `md.h`).
- `include/atomiverse/cpu/mc/translate_move.h` (new) — `cpu::mc::TranslateMove` Mover declaration.
- `src/cpu/mc/translate_move.cpp` (new) — `TranslateMove::launch` implementation.
- `CMakeLists.txt` (modify) — add `src/cpu/mc/translate_move.cpp` to the `atomiverse` source list (next to the CPU MD components block near line 167).
- `tests/cpu/mc/test_mc_move.cpp` (new) — core-correctness unit tests.
- `tests/CMakeLists.txt` (modify) — register `atv_add_cpp_test(test_cpu_mc_move cpu/mc/test_mc_move.cpp)`.

## Tasks

- [ ] Write failing tests for MC core correctness — bootstrap, exact rollback, e_cur consistency, Metropolis exactness, determinism, acceptance-ratio trend, and build-validation throws (tests/cpu/mc/test_mc_move.cpp)
- [ ] Implement TranslateMove Mover (role==Mover, caps==None, δ-validating ctor) in include/atomiverse/cpu/mc/translate_move.h + src/cpu/mc/translate_move.cpp
- [ ] Implement cpu::mc::MC driver in include/atomiverse/cpu/mc/mc.h — ZeroEnergy head-bind at sampling::Phase::energy −999, set_mover via set_singular_default_, T/seed-validating ctor, owns rng_/proposal_/e_cur_ published to Context
- [ ] Add Doxygen (@brief/@param with units/@return) to every public symbol in mc.h and translate_move.h, plus the assembled-usage example comment in mc.h
- [ ] Register src/cpu/mc/translate_move.cpp in CMakeLists.txt and tests/cpu/mc/test_mc_move.cpp in tests/CMakeLists.txt
- [ ] Run full check + test suite

## Testing strategy
Happy path
- Assembled `cpu::mc::MC` with `TranslateMove` + `PairLJCut` + `NeighborList` builds and runs `run(frame, n_sweeps)` without error.

Core correctness (reference-free, fast)
- Bootstrap: after `build()` / first evaluation, `e_cur_` equals a fresh full-Frame energy of the initial configuration to 1e-12.
- e_cur consistency (accept): after an accepted move, `e_cur_` equals a fresh full recompute of the new state to 1e-12.
- Exact rollback (reject): after a rejected move, `e_cur_` is unchanged AND every frame position is bitwise/1e-15 identical to the pre-move state (restored from `Transaction.old_pos`).
- Metropolis exactness: a constructed downhill/degenerate move (ΔE ≤ 0) is accepted unconditionally (probability 1); at `T → 0` only non-uphill moves are accepted.
- Determinism: identical seed ⇒ identical accept/reject sequence and identical final trajectory; a different seed produces a different sequence.
- Acceptance-ratio sanity: the accepted-fraction decreases monotonically as `δ` grows, and a mid-range `δ` lands in a broad band (~20–60%); assert the monotonic trend and the loose band only — never a point value (tuning heuristic, not a theorem).

Edge cases
- `set_mover` re-called with a different Mover throws; `MC{T, seed}` with `T ≤ 0` throws; `TranslateMove{δ}` with `δ ≤ 0` throws.
- Build validation: omitting the Mover, or omitting any `Capability::Energy` evaluator, throws at `build()` (exercises `SamplingDriver::validate_components` through the concrete driver).

## Out of scope
- Physics distribution gates and the NIST LJ end-to-end validation → sub-spec 03.
- Local / incremental ΔE (v1 always full-Frame recompute) and adaptive-δ tuning.
- Additional move types (rotate / swap / insert / delete), GCMC, REMD.
- CUDA backend, swappable acceptance policy, and external-potential paths.
