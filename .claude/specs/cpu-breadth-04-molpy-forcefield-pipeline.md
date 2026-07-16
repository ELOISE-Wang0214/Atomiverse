---
title: CPU breadth 04 — molpy/molrs ForceField to native Atomiverse Plan
slug: cpu-breadth-04-molpy-forcefield-pipeline
status: approved
execution: queued
last_progress_update: 2026-07-16
science_required: true
created: 2026-07-15
acceptance: cpu-breadth-04-molpy-forcefield-pipeline.acceptance.md
depends_on: cpu-breadth-03-charged-npt-kspace
---

# CPU breadth 04 — molpy/molrs ForceField to native Atomiverse Plan

## Summary

Deliver the first real ecosystem consumer pipeline:

`molpy build/typify -> molrs Frame + ForceField -> Atomiverse native CPU
components -> run -> same Frame -> molrs analysis`.

The initial fixture is two distorted, flexible TIP3P water molecules. It is
large enough that bond, angle, intermolecular LJ, and intermolecular Coulomb
terms are all non-zero, while deliberately avoiding torsion, improper, RATTLE,
and KSpace in this link. The user does not insert integer types or register
force components by hand.

This probe tests ownership, type resolution, units, mixing, special-bonds,
fail-closed style support, and build/launch separation. It must compile to
ordinary Atomiverse components; a hidden monolithic "ForceField component"
would evade the architecture question and is forbidden.

## Progress (2026-07-16)

**Execution: queued as the next implementation link — 0/10 tasks and 0/11
active acceptance criteria verified. GPU remains deferred.**

The prerequisite Frame/Element hard cut is present. Work specific to this link
has not started: the molrs ForceField snapshot/CXX bridge, CPU PairCoulCut,
deterministic label-to-ID compiler, owned `CompiledForceField`, unit conversion
table, public two-water regression, and release chain all remain open.

## Domain basis

molrs is the single owner of ForceField style/type semantics, combination
rules, `special_bonds`, and canonical parameter units. For the accepted
snapshot those units are explicitly carried, not inferred from a Python
constructor string. Atomiverse converts once at compile/build time to its
engineering units.

The convention mismatch is load-bearing:

- molrs harmonic terms use `1/2 k (x-x0)^2`;
- Atomiverse BondHarmonic/AngleHarmonic use `K (x-x0)^2`;
- therefore `K_ATV = 1/2 k_molrs * kcal_mol_to_eV` for the accepted molrs
  kcal/mol parameterization;
- LJ epsilon converts kcal/mol -> eV; sigma/r0 remain Å; angle equilibrium and
  derivatives use radians; charges remain elementary-charge units.

The first supported style matrix is exactly:

- atom style sufficient for typed charged atoms (`atom/full` semantics);
- `bond:harmonic`;
- `angle:harmonic`;
- `pair:lj/cut`;
- `pair:coul/cut`.

CPU `PairCoulCut` is delivered by this link. Substituting `coul/long` or PME for
a ForceField that requested `coul/cut` is a forbidden fallback. Scientific
acceptance uses molrs Potentials as an independent Rust-kernel oracle, centered
force finite differences, and a pinned LAMMPS validation case.

## Recon / current evidence

- Before frame-block 04, molpy re-exported `molrs.Frame` as `molpy.Frame`.
  That compatibility surface is deleted by the prerequisite; this link uses
  the direct `molrs.Frame` returned by current molpy builders. Atomiverse has a
  smoke test that manually inserts an integer `type`, registers
  LJ/NeighborList, and runs ten MD steps.
- molrs owns the ForceField hierarchy and exports a Python
  `_ffi_forcefield_capsule()` carrying a `molrs.ForceFieldRef`. molrs-cxxapi
  does not yet enable the `molrs-ffi/ff` feature or bridge that handle.
- The capsule documentation says shared/no-deep-copy, but its current Python
  implementation clones the ForceField before wrapping it. This spec chooses
  and documents immutable snapshot semantics rather than shared mutation.
- Canonical molrs topology/type columns can be strings; Atomiverse
  materialization currently expects dense integer `type` values.
- Atomiverse does not consume the generic molrs ForceField and has no plain CPU
  `coul/cut` pair style.
- Existing special-bonds handling is insufficient evidence for arbitrary
  non-zero 1-4 scaling, and long-range reciprocal exclusion correction is not
  implemented. The first fixture exercises 1-2/1-3 exclusion and explicitly
  fails on unsupported 1-4/kspace styles.

## Architecture questions

1. Can molrs remain the sole schema/ForceField owner while Atomiverse receives
   a stable, read-only, backend-neutral compile snapshot?
2. Can string labels become deterministic engine-local dense IDs without
   mutating or replacing canonical Frame columns?
3. Can one compiled program install ordinary native components into Evaluate,
   MD, and Minimize without preset-specific translation branches?
4. Can all FFI, string lookup, unit conversion, and combination work stop at
   compile/build so the hot loop is native and allocation-free?

## Design

### Three-repository ownership

- **molpy** builds/typifies through public user APIs and owns only Python
  convenience. It does not define a parallel ForceField or simulation kernel.
- **molrs** owns Frame/ForceField, style parsing, type labels, mixing rules,
  special-bonds, units, and the immutable resolved snapshot.
- **Atomiverse** consumes the snapshot once, resolves dense CPU types, converts
  units, constructs native NeighborList/Pair/Bond/Angle components, and owns
  their runtime state.

Atomiverse never depends on molpy. The optional cross-repo Python test imports
molpy as a producer, then calls Atomiverse's public adapter.

### Versioned immutable ForceField snapshot

molrs-cxxapi adds a CXX wrapper around `molrs_ffi::ForceFieldRef` and enables
the `ff` feature. Under the Python GIL it resolves an immutable snapshot
containing:

- API/schema version and ForceField fingerprint;
- category/style and style parameters;
- ordered type labels and fully resolved type/pair parameter tables, including
  the selected mixing rule;
- canonical unit identifiers;
- `special_bonds.lj[3]` and `.coul[3]`;
- deterministic mapping information needed to resolve Frame topology/type
  labels, with every missing label reported by block, row, label, and style.

The `ForceFieldRef` is `!Send`; neither it nor Python objects enter Driver state
or worker threads. Snapshot creation may copy compact parameter tables. Later
ForceField mutation does not affect a compiled program; users explicitly
recompile. Documentation and tests are corrected to state this snapshot
semantics.

### Native compiled program, not a composite evaluator

An Atomiverse `CompiledForceField`/equivalently named non-Driver,
non-Component object owns:

- a Frame-ingress dense type resolver/fingerprint;
- unit-normalized parameter tables;
- ordinary native components: exactly one shared NeighborList, PairLJCut,
  PairCoulCut, BondHarmonic, and AngleHarmonic for the first fixture.

It installs owned components through the general Driver/phase-traits
registration path established in link 01. It has no `if (MD)` / `if
(Minimize)` materialization branches; Evaluate, MD, and Minimize select their
typed evaluator phase at the call boundary. Production hot paths contain no
molrs Potentials call, CXX call, string lookup, JSON/dict conversion, or
HashMap resolution.

`cpu::materialize(source, resolved_types)` consumes the resolver once while
leaving the source Frame's string columns unchanged. On reuse, a changed
ForceField or type/topology fingerprint fails with an explicit "recompile"
error rather than running stale parameters.

### Fail-closed supported set

Every ForceField style is inspected. Unsupported styles are collected into one
error listing `category:name` and the supported matrix. Nothing is silently
dropped. In this link, `dihedral:*`, `improper:*`, `pair:coul/long`, `kspace:*`,
non-trivial 1-4 scaling, and constraints are unsupported and must fail before
the Driver builds.

The second planned fixture is OPLS butane, but it is not smuggled into this
slice. Link 07 decides the next CPU style expansion after the minimal pipeline
has evidence.

### Release order

1. molrs releases the ForceField snapshot plus a new CXX API version.
2. Atomiverse pins that released tag/SHA and exact Python package version,
   implements the consumer, and passes the cross-repo contract.
3. molpy pins the same molrs release and lands the public integration example.

Current dirty feature-branch SHAs are not release pins.

## Files to create or modify

- molrs `molrs-ffi`/molrs-cxxapi ForceFieldRef bridge and snapshot schema/tests
- molrs-python capsule semantics docs/tests and exact version metadata
- Atomiverse `cpu/io/` bridge boundary and `cpu/potentials/` compiled program
- CPU `PairCoulCut` component and tests
- Frame materialization type-resolver ingress without source-column mutation
- Python `atomiverse.forcefield` adapter/bindings/docs
- molpy public integration test/example/docs (no Atomiverse hard dependency)
- cross-repo acceptance target and release pins
- one TIP3P pipeline regression/validation fixture and `REFERENCE.md`
- `.claude/notes/notes.md`

## Tasks

- [ ] Write molrs snapshot/capsule contract tests, including version mismatch,
      snapshot-after-mutation semantics, units, mixing, special-bonds, and
      aggregated unsupported styles.
- [ ] Write the two-water public Python acceptance and molrs/LAMMPS energy-force
      oracles before implementing the adapter.
- [ ] Add the ForceFieldRef CXX bridge at its Rust source and increment the CXX
      API version; never hand-edit generated bridge output.
- [ ] Implement CPU PairCoulCut with energy, force, virial, cutoff/dielectric
      parameters, NeighborConsumer protocol, and no launch allocation.
- [ ] Implement deterministic string-label/dense-ID resolution without
      modifying the molrs Frame's canonical columns.
- [ ] Implement the native compiled program and generic owned installation into
      Evaluate, MD, and Minimize.
- [ ] Add all unit conversions to the shared units source; grep production
      adapter code for magic kcal/eV constants.
- [ ] Enforce fail-closed style/fingerprint behavior and prove the hot loop has
      zero FFI/string lookup/allocation after build.
- [ ] Run molrs Rust tests, Atomiverse CPU unit/regression/validation/Python
      tests, and the molpy public integration case in one exact-version
      environment.
- [ ] Release/pin in molrs -> Atomiverse -> molpy order and record the boundary
      decision in `.claude/notes/notes.md`.

## Testing strategy

- **molrs unit**: ForceField snapshot schema, capsule lifetime/snapshot
  semantics, units, resolved tables, and source gates.
- **Atomiverse unit**: PairCoulCut, conversion conventions, dense resolution,
  Plan contents, fail-closed unsupported/fingerprint paths, and no FFI in
  launch.
- **Cross-repo regression**: two distorted flexible waters constructed through
  public molpy APIs; same Frame identity, no manual integer type/component
  setup, four non-zero energy terms, ≤60 seconds.
- **Validation**: per-term energy and all forces versus molrs Potentials and a
  pinned LAMMPS input; centered finite difference independently checks forces.
- **Performance**: compile/build allocations are reported; steady native launch
  performs zero FFI, string/hash lookup, and allocation.

## Deferred — GPU

- Device-resident ForceField snapshot/materialization and CUDA components.
- CUDA PairCoulCut and CPU/CUDA Plan/scientific/performance parity.
- multi-GPU type resolution, special-bonds communication, and migration.

No CUDA/Torch source is touched. CPU `coul/long + PME`, OPLS/improper styles,
and non-trivial 1-4 scaling are also outside this minimal link, but are CPU
backlog rather than GPU-deferred evidence.

## Out of scope

- A molpy subprocess-style Engine subclass or temporary-file input generation.
- JSON/dict/on-disk ForceField IR.
- Calling molrs Potentials from Atomiverse production kernels.
- Rigid water/RATTLE, dihedrals, impropers, KSpace, polarizable styles, or ML
  potentials.
- Mutably sharing a Python ForceField with a running Driver.
