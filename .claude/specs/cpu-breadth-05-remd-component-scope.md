---
title: CPU breadth 05 — REMD and explicit per-replica component scope
slug: cpu-breadth-05-remd-component-scope
status: approved
execution: queued
last_progress_update: 2026-07-16
science_required: true
created: 2026-07-15
acceptance: cpu-breadth-05-remd-component-scope.acceptance.md
depends_on: cpu-breadth-04-molpy-forcefield-pipeline
---

# CPU breadth 05 — REMD and explicit per-replica component scope

## Summary

Implement the first real budget-terminating replica method: CPU replica
exchange molecular dynamics (REMD). It must use `IdentityReplication`, one
Frame input/output, an outer replica loop, nested MD propagation, and adjacent
temperature exchanges.

The feature exists to force an architecture answer that NEB intentionally
left open: stateful components cannot be shared across replicas, while
`ComponentBase` has no general clone operation. This link chooses explicit
per-replica program construction and instance scope. It does not add a blanket
`clone()` virtual or require every component to carry a map keyed by replica.

## Progress (2026-07-16)

**Execution: queued after link 04 — 0/9 tasks and 0/10 active acceptance
criteria verified. GPU remains deferred.**

No REMD implementation is claimed. Existing ReplicaSet, IdentityReplication,
SubDriver, and NEB code are evidence inputs only; ReplicaExecutor,
per-replica program construction, independent stateful component trees,
exchange transactions, deterministic continuation, and acceptance tests remain
to be implemented.

## Domain basis

For configurations with potential energies `U_i`, `U_j` at inverse
temperatures `beta_i`, `beta_j`, the adjacent replica-exchange probability is

`p = min(1, exp[(beta_i-beta_j)(U_i-U_j)])`.

Only potential energy enters this exponent. Temperature slots remain fixed.
On an accepted configuration swap, coordinates/configuration state exchange;
velocities are rescaled by the square root of the temperature ratio so the
receiving slot has the correct kinetic distribution. Thermostat state remains
with its temperature slot. References: Hukushima & Nemoto, JPSJ 65, 1604
(1996), and Sugita & Okamoto, CPL 314, 141 (1999).

The public regression uses a one-particle analytic harmonic/double-well model
with four overlapping temperatures. It checks temperature marginals and a
configuration round trip with a fixed seed. Units follow the CPU MD contract:
energy [eV], temperature [K], time [fs].

## Recon / current evidence

- `ReplicaSet` and `IdentityReplication` already model one Frame expanded to N
  parameterized replicas; `IdentityReplication` has no production caller.
- `PathDriver` contains correct per-replica/band-wide scheduling for NEB but
  hard-codes a convergence loop and one shared component Plan.
- `SubDriver` correctly isolates nested loop status on a borrowed Runtime, but
  one SubDriver owns only one inner component set.
- NeighborList, thermostat, RNG, and optimizer components hold mutable state.
  Reusing one instance sequentially across replica Frames is correct only for
  stateless components and can invalidate caches or correlate random streams.
- `ComponentBase` has no clone contract, and many components are intentionally
  non-copyable. Adding a dead clone method to every component is not justified.
- Path setup currently learns replica count after ordinary component build,
  which has already caused first-launch allocation for band workspaces.

## Architecture questions

1. What is the explicit instance scope of a component: shared, per loop, or
   per replica?
2. Can N owned component/driver trees execute sequentially on one CPU Runtime
   without constructing N Backends or leaking child status into the parent?
3. Can NEB and REMD share the replica scheduler while retaining different
   termination and coupling policies?
4. Can replica count/capacity be known before build so every per-replica
   workspace is allocated outside launch?

## Design

### Per-replica program builder, not ComponentBase clone

Introduce an OOP builder/program protocol in `cpu/replica/` whose instance owns
the construction strategy for one nested executable node. At outer setup,
after `ReplicaExpansion` fixes the replica descriptors, it constructs N owned
MD driver/component trees exactly once. Each concrete tree has independent
integrator, thermostat, RNG, NeighborList, evaluator scratch, handlers, and
Plan. Stateless immutable parameter tables may be shared explicitly.

The protocol is a class with constructor-established configuration, not a
`make_remd()` free function and not a `std::function` hidden in Context. It
returns executable nodes at the node boundary; each child preserves typed MD
phases internally. A user provides one program configuration, never N Frames or
N manually assembled drivers.

All children borrow one outer CPU Runtime through nested execution. They do not
create N domain decompositions, communicators, timing registries, or Backends.
The first implementation is single-rank sequential CPU; replica-parallel MPI
groups are a later CPU extension.

### Shared replica executor

Extract from `PathDriver` the generic scheduling primitive that:

- owns/borrows a `ReplicaSet` generated from one Frame;
- makes replica descriptors and count available before child build;
- executes one per-replica node/stage and one band-wide `couple` stage;
- restores active Frame and loop services on every normal/exception exit;
- collapses through the expansion strategy.

PathDriver continues to attach `AllReplicasConverged`; REMD attaches a budget
and exchange diagnostics. REMD is a preset/method in `cpu/path/`, not a new
Driver category. Existing NEB source-level public usage and regression remain
unchanged.

### Temperature slots and exchange

`IdentityReplication(temperatures)` defines finite, positive, strictly
increasing temperature slots. Each exchange interval runs the nested MD child
for every replica, then an exchange component/policy bound at the generic
`path::Phase::couple` attempts adjacent pairs. Schedules alternate even and odd
edges, so one replica exchanges at most once per round.

RNG streams are deterministically derived from a master seed and domain label:
one independent stream per propagator plus one exchange stream. `random_device`
and a process-global generator are forbidden. Diagnostics expose attempted and
accepted counts per edge, exchange history, and temperature/configuration
round trips.

On rejection both Frames and slot state remain bitwise unchanged. On
acceptance configuration state swaps, velocities rescale exactly, thermostat
state stays with the fixed temperature slot, and every configuration-dependent
cache is explicitly invalidated/rebuilt at its declared lifecycle boundary.

### Continuation and collapse

The REMD object is continuable: its first run expands and builds the replica
programs; subsequent advance calls on the same object retain the band, RNG, and
statistics. The caller receives the collapsed lowest-temperature slot in the
same Frame after each call. Replacing the input system requires an explicit
reset; a changed system fingerprint cannot silently reuse an old band.

`run(frame, budget)` remains the simple public entry. Internally, chunked
advance and one-shot advance with the same seed/budget must produce identical
state. Destroy releases every child exactly once.

## Files to create or modify

- `include/atomiverse/cpu/replica/replica_executor.h` and implementation
- `include/atomiverse/cpu/replica/replica_program.h` and implementation
- `include/atomiverse/cpu/path/remd.h` and `src/cpu/path/remd.cpp`
- `cpu::PathDriver` refactor onto the shared executor without public API change
- exchange policy/component, diagnostics, deterministic RNG derivation
- unit tests under `tests/cpu/path/` or `tests/cpu/replica/`
- `regressions/path/remd_harmonic/` (new, with `REFERENCE.md`)
- lifecycle/performance benchmark and CMake registration
- `.claude/notes/notes.md`

## Tasks

- [ ] Write failing instance-scope tests that expose shared thermostat/RNG/
      NeighborList state and assert distinct child/component addresses.
- [ ] Write exchange formula, rejection transaction, accepted swap/rescale,
      even/odd schedule, RNG derivation, and continuation tests.
- [ ] Write the public harmonic/double-well REMD regression before extracting
      or generalizing replica scheduling.
- [ ] Extract the minimal shared ReplicaExecutor and keep the existing NEB
      public regression source unchanged.
- [ ] Implement the per-replica program builder and N owned nested MD trees on
      one borrowed CPU Runtime; do not add `ComponentBase::clone()`.
- [ ] Implement fixed temperature slots, exchange diagnostics, cache
      invalidation, lowest-slot collapse, reset, and chunked continuation.
- [ ] Move all N-dependent child/workspace allocation to setup/build and add
      memory/rebuild counters.
- [ ] Run unit, REMD regression, unchanged NEB regression, determinism,
      allocation, and performance gates.
- [ ] Record the local instance-scope mechanism and evidence in
      `.claude/notes/notes.md`; link 07 decides whether it becomes a wider protocol.

## Testing strategy

- **Unit**: builder/child ownership, nested status isolation, fixed-energy
  exchange probability, exact reject/accept state changes, schedule, RNG,
  continuation/reset, cache invalidation, and lifecycle counters.
- **Regression**: four-temperature one-particle analytic model, public API,
  marginal energy and at least one deterministic round trip, ≤60 seconds.
- **Non-regression**: existing CI-NEB regression and SubDriver isolation tests
  pass without source changes.
- **Performance**: setup and O(N) memory reported separately; after build,
  propagation/exchange allocates zero unless an explicit cache rebuild occurs.

## Deferred — GPU

- CUDA child programs, device-resident exchanges, RNG, and cache invalidation.
- single-/multi-GPU REMD and CPU/GPU statistical/performance parity.
- multi-GPU execution through the shared CUDA Driver/Plan.

GPU work is a future campaign. Replica-parallel CPU MPI communicator groups
are also outside this first link, but remain CPU backlog rather than GPU work.

## Out of scope

- PIMD, solute tempering, Hamiltonian exchange, multidimensional ladders.
- User-supplied `vector<Frame>` or one Backend per replica.
- A mandatory clone interface or replica-indexed mutable state in every
  component.
- Replica-parallel MPI, checkpoint file format, or distributed exchange.
