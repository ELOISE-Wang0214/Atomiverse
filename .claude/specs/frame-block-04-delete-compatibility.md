---
title: Frame block 04 — atomic physical cutover and compatibility deletion
slug: frame-block-04-delete-compatibility
status: approved
execution: in_progress
last_progress_update: 2026-07-16
science_required: false
created: 2026-07-15
acceptance: frame-block-04-delete-compatibility.acceptance.md
depends_on: frame-block-03-backend-audit (closed; history in git log)
---

# Frame block 04 — atomic physical cutover and compatibility deletion

## Summary

Complete the signed frame-block verdict in one breaking CPU cutover. After this
link, the physical CPU Frame is exactly `{blocks: string -> Block, meta: string
-> MetaValue, simbox}`. The current `cpu::Atom`, `SystemScalars`, and `Box`
members are not retained behind that API, and no alias, proxy view, deprecated
overload, dual write, old wire decoder, or cross-repo translation keeps the old
representation alive.

This is one atomic link, not another staged migration. Every in-tree CPU
consumer, test, regression, example, Python binding, launcher path, and the
required molrs/molpy boundary changes in the same implementation. Old source,
ABI, Python imports, and serialized metadata are expected to break loudly.

CUDA and Torch physical cutover are explicitly deferred. Their current source
is not edited, their compatibility code is not counted as CPU success, and no
CPU compatibility layer may be added to shelter them.

## Progress / current evidence (2026-07-16)

**Execution: in progress — 9/15 tasks implemented; 4/17 active acceptance
criteria verified. GPU/Torch remains deferred.**

- `cpu::Frame` now physically owns exactly `BlockMap blocks`, `MetaMap meta`,
  and `Box simbox`. The old `cpu::Atom` storage files, named `atoms/scalars/box`
  members, proxy meta views, and Box compatibility operators are deleted.
- CPU production, app, Python, test, regression, example, benchmark, tool, and
  documentation callers were migrated to typed spans/meta and `simbox` in the
  same cut. The checked compatibility inventory, direct-access scrub, and
  registry field-list scrub pass on the current tree.
- Schema-driven native-dtype pack plans now carry generation and canonical
  fingerprints. Runtime columns participate in grow/copy/exchange/broadcast,
  typed runtime metadata preserves exact dtypes, and bigbig `TagInt(2^53+1)`
  is exercised without a `double` funnel.
- molrs now owns typed `MetaValue` metadata across Rust, serialization, C,
  Python, and CXX; Frame schema is exactly v2 and CXX API is exactly v1 with
  frame-block, root `Element`, and AM1-BCC capabilities. Non-current serialized
  Frame versions fail rather than translate.
- molpy imports `Frame`, `Block`, and `Element` directly from molrs where used.
  Its former `molpy.core.frame` module and facade re-exports are absent and
  import-failure tests pin the break.
- TaskSpec is strict v3 with explicit `ComponentId`; the v2 decoder and
  variant-index identity are deleted. Audited solver wrappers, flat numeric
  headers/aliases, constants, and compatibility overloads are also deleted.
- Still blocking closure: explicit stale-handle/ASan mutation fixtures, the
  complete precision/ID/MPI/OpenMP matrix, exact released cross-repository
  pins and full suites, and the resolution/allocation/memory/performance and
  release-documentation gates.

## Binding architecture decision

There is one current API and one current representation. A change replaces the
old form; it never wraps, aliases, decodes, mirrors, or deprecates it.

The following are forbidden in this link and in all later specs:

- `[[deprecated]]` bridges, forwarding aliases, compatibility modules, and old
  namespace re-exports;
- old/new dual storage, dual read, dual write, lazy conversion, and fallback
  key lookup;
- decoding an older TaskSpec, Frame, metadata, checkpoint, or Python payload;
- keeping an old free-function/overload solely to construct the new object;
- silently accepting an old key, dtype, import path, constructor spelling, or
  behavior and translating it to the new one.

An exact current version marker is allowed only to reject every non-current
version. Format-specific canonicalization at a real external I/O boundary
(for example, reading a LAMMPS `q` field into the canonical `charge` column) is
not backward compatibility and remains local to that format reader.

## Design

### Physical CPU Frame is the mapping

Replace facade-over-members storage with three physical regions:

- `blocks`: owns backend-local CPU Blocks and Columns. The canonical `atoms`
  block contains x/v/f, mass, charge, type, tag, image, per-atom topology, and
  runtime-added columns. A Column owns its xtensor storage through the dtype
  variant; there is no separate `ExtraColumn` class and no named xtensor is
  mirrored elsewhere.
- `meta`: owns all non-column bookkeeping through stable typed values/slots and
  mutable typed references, including population counts, topology strides, all
  system scalars and both stress arrays. The engine-local tag-to-row index is a
  typed entry in this same region, tagged never-serialize. There is no second
  bookkeeping aggregate or proxy map.
- `simbox`: owns the canonical cell. The Frame has no `box` alias and callers
  mutate the cell only through the current Box API, which refreshes its cache.

Delete `Frame::atoms`, `Frame::scalars`, `Frame::box`, the proxy-view macros,
and `cpu::Atom` as a Frame storage container. Delete Box's legacy element and
assignment operators, and delete the old `Frame::insert -> Atom::add_column`
forwarder; runtime columns enter through the canonical Block insertion API. A
named accessor such as `cell()` may remain only if it is the canonical current
Box API with a real caller; it must not preserve the old `Mat3` calling
convention. Const and mutable typed resolution must be implemented directly,
with no `const_cast` back through the old storage.

`core/frame_schema.h` remains the only field/schema enumeration, but its rows
no longer carry a C++ member name used to bind `cpu::Atom`. It drives Block
construction, grow/shrink, row copy, migration, halo/exchange payloads,
materialize/write-back, memory accounting, and runtime-added column behavior.
Per-atom CPU topology remains a Block-column layout with runtime strides;
molrs flat topology is converted exactly once at the materialization boundary.

Schema-driven does not mean interpreted inside a row loop. Grow/copy/comm paths
compile narrow column/pack plans by resolving keys and variants once at
build/rebuild; the per-row loop receives typed pointers, strides, and semantic
flags only. A capacity or topology re-stride increments a storage generation
and rebuilds every affected plan before use, so no cached span survives
reallocation. OpenMP workers share immutable plans and write only their
declared disjoint/reduced regions.

Wire plans preserve each native dtype; no generic `double get/set` funnel is
allowed. `StorageReal`, `AccReal`, int8/int32, and 32/64-bit `TagInt` pack and
round-trip bit-exactly, including tags above `2^53`. Before any MPI payload,
ranks collectively compare a canonical schema fingerprint over version,
sorted key, dtype, stride, and semantic flags. Different insertion order for
the same runtime columns yields the same fingerprint/wire order; a true
mismatch fails collectively before communication. Any column-count limit is
validated at registration, never truncated during pack.

### All CPU consumers cut over together

Every real `cpu::Frame` access under `include/atomiverse/cpu`, `src/cpu`,
`tests/cpu`, `regressions`, `examples`, `app`, and `python` moves directly to
typed `span(ColKey)` / typed meta handles / `simbox`. Generic boundary code may
resolve dynamic keys. `build()` may cache logical keys/slot IDs and immutable
plan structure, but never a raw Span/MetaRef across a possible grow, re-stride,
migration, or runtime-column insertion. Raw handles resolve at the safe launch
boundary after structural mutation (including meta insertion), or
generation-check and refresh before use; inner loops receive validated spans
or raw pointers only.

The migration includes communication and nested execution, not just force
kernels: materialize/write-back, broadcast/gather, atom migration, halo
exchange, neighbor construction, MC transactions, replica expansion, path,
minimizers, SCF/PW/GPW evaluators, observers, recorders, TaskSpec, and Python
bindings all use the same physical mapping. Tests and examples are rewritten
to teach only the new API; no test fixture reconstructs the old struct.

### Typed meta cuts over across repositories

molrs changes `Frame.meta` to the one current typed-meta semantic schema in its
Rust core, views, serialization, C/Python/CXX APIs, and tests. Exchangeable
entry kinds and their wire IDs are exact across repositories; backend-local
representations remain independent, and Atomiverse's tag index is explicitly
never serialized. The old `HashMap<String,String>` representation and
stringification fallback are deleted. Serialization writes one new exact
format version and rejects every other version; it has no historical decoder.

Atomiverse materialize/write-back and CXX FFI consume that exact released API.
There is no typed-to-string downgrade. molpy updates its in-tree consumers to
import `molrs.Frame`/`molrs.Block` directly and deletes the compatibility-only
`molpy.core.frame` module and its top-level Frame/Block re-exports. No import
hook or `__getattr__` restores those names. The three repositories release in
dependency order and pin the exact same current contract.

### Delete the remaining CPU compatibility surface

The implementation starts from a checked inventory of production declarations
whose documented purpose is old-call-site, ABI, namespace, behavior, or wire
compatibility. Every such declaration in the CPU scope is removed in this same
link and every in-tree caller moves to its sole replacement. The initial known
set includes:

- `cpu::md::ZeroForces/ZeroEnergy/ZeroStress/ZeroKinetic` namespace aliases;
  `ZeroMomentum` becomes the sole concrete current class and the historical
  `ZeroLinearMomentum` spelling is deleted rather than retained behind an alias;
- all three transient-solver wrappers `scf_solve()`, `pw_scf()`, and
  `gpw_scf()`; in-tree low-level callers construct `SCFSolver`, `PWSolver`, or
  `GPWSolver`, while workflow callers use the matching Evaluator;
- the `rhf_solve`/`uhf_solve`/`uks_solve` adapter layer and its duplicate result
  types; method-specific config builders remain current typed configuration;
- the no-`GridLevel` molecular-grid overload (callers select `GridLevel::SG0`
  explicitly if that current scientific mode is intended), the unused UPF
  `KBProjector::d` diagonal compatibility cache, and the dead `Mat33` alias;
- TaskSpec's `variant.index()` wire identity, the “AM1 must be last” ordering
  constraint, and its compatibility test. The current wire version is bumped,
  uses explicit current `ComponentId` values, and rejects the previous version;
- shared flat compatibility headers/types/constants/log overloads used by CPU
  (`core/constants.h`, `core/real.h`, `core/int.h`, `Real`, `Int`, flat `k*`,
  and old `log_init` overloads). Numerically distinct current constants such as
  QE-rounded conversions move to an explicit owning namespace before the old
  spelling is deleted;
- Python re-export modules, flat aliases, permissive version warnings, and old
  payload/key fallbacks in the Atomiverse/molecular-data boundary.

Scientific methods that are independently named current choices are not
compatibility code. SG0, `LineSearch::Fixed`, explicit linear mixing,
canonical/Löwdin orthogonalization, and method-specific SCF defaults may stay
only as explicitly selected and tested current algorithms—not as implicit
fallbacks chosen to preserve old behavior. `AM1Config`/`AM1Result` remain the
intentional AM1 scientific facade rather than an NDDO migration alias.
Ambiguous cases are resolved by deleting the old spelling and making the
desired current method explicit.

### Performance and failure behavior

- Every accepted hot component declares its exact resolution bound. Debug
  tests run at N=1, 64, and at least 100,000 and prove Block/meta resolution is
  O(1) in atoms and pairs, including nested replay/sub-drivers. The diagnostic
  counter is thread-safe under a deliberately bad OpenMP component, reports
  only after leaving the parallel region, and has no TSAN race. All resolution
  instrumentation compiles out of Release.
- Zero steady allocation is a hard gate for the affected classical CPU
  Frame/comm/neighbor/pair/integrator paths, not a false claim about every
  electronic-structure solver. The interceptor is thread-safe and covers
  `new/new[]`, malloc/calloc/realloc, aligned_alloc, and posix_memalign, reporting
  count and bytes. It covers one/many OpenMP threads, reused/forced neighbor
  rebuild, and pre-reserved two-rank halo/migration. First build or explicit
  capacity growth is measured separately. SCF/PW/GPW may retain existing solver
  allocations but the Frame cutover cannot increase their baseline count/bytes.
- The LJ gate uses two fixed single-rank sentinels: steady neighbor reuse and
  forced reneighbor. The pre-cutover baseline records git SHA, host CPU, ISA,
  compiler/flags, xsimd version/batch width, N, thread/BLAS settings, warmup,
  cadence, and at least ten repeated samples; measured noise must be below 1%.
  Each sentinel independently fails above 2%—no geometric-mean or 8% fallback.
- materialize and write-back are timed separately for atomic and bonded Frames
  (ns/atom, GB/s, allocation/peak bytes). Pack/unpack, broadcast, steady halo,
  halo rebuild, forced migration, and gather are separate 1/2-rank cases. Each
  has a same-runner pre-cutover baseline and a 5% ceiling; absence of a baseline
  is failure, not skip. Plan build/rebuild time is not averaged into steady work.
- Performance runs pin `OMP_PROC_BIND=close`, `OMP_PLACES=cores`, thread count,
  and BLAS threads=1 and record all toolchain/runtime identities. The memory
  ledger separates blocks/meta/simbox, capacity/used bytes, comm scratch, and
  materialize source/destination; peak RSS alone is not sufficient.
- CPU compile-command acceptance preserves language-guarded `-march=native`
  and proves `XTENSOR_USE_XSIMD` is not defined globally. Existing xsimd hot
  translation units consume resolved contiguous spans; the migration does not
  replace them with xtensor expression or type-erased inner loops. CI pins and
  reports the actual xtensor/xsimd pair used by CMake and runs an XC-grid
  `reshape_view` smoke test; it does not hard-code an outdated skill version.
- Because this link removes SCF/PW/GPW functional wrappers and migrates their
  Frame-facing evaluators, it reruns `BM_GBuild_Block_Aspirin` (target no slower
  than 950 ms and <=5% versus the pinned pre-cutover runner) plus
  `BM_GBuild_Block_Aspirin_ccPVDZ` (<=5% versus its pinned baseline), and keeps
  the Rys d-shell numerical cross-check at 1e-7.
- Unknown blocks/columns/meta keys, wrong dtypes, missing required columns, old
  versions, and old Python/C++ spellings fail immediately. They never trigger a
  conversion or alternate path.
- No raw-memory broadcast/serialization of C++ Frame/Block/Box/variant objects
  is permitted; wire payloads are schema-driven current-version data.

### Required source scrub

The link adds a checked scrub target, scoped to real CPU Frame consumers so
unrelated scientific structures named `atoms` do not create false positives.
At minimum it rejects `cpu::Atom`, `cpu/frame/atom.h`,
`ATV_FRAME_META_(INT|REAL)_VIEWS`, `ATV_FIELD_BIND`, `ExtraColumn`, `.extras`,
`const_cast<Frame*>`, `const_cast<Atom*>`, `.scalars.`, and direct
`cpu::Frame` access to named atom columns/counts/topology/capacity/tag-map
methods. Compile-negative public consumers independently prove that
`cpu::Frame` has no `atoms`, `scalars`, or `box` member and that atom-pack has
no `Atom&` overload. The compatibility inventory is a checked data file, not a
one-time reviewer checklist.

## Files to create or modify

- `include/atomiverse/core/frame_schema.h`
- `include/atomiverse/cpu/frame/` and `src/cpu/frame/`
- all real CPU Frame consumers under `include/atomiverse/cpu` and `src/cpu`
- CPU communication, TaskSpec/launcher, Python binding, tests, regressions,
  examples, docs, and source scrub gates
- `AGENTS.md` and architecture/API documentation that currently promise old
  functional shims or describe named Frame storage
- explicitly compatibility-only CPU public headers/implementations discovered
  by the audit
- `/Users/roykid/work/molcrafts/molrs/molrs/`, `molrs-ffi/`, `molrs-capi/`,
  `molrs-cxxapi/`, and `molrs-python/` typed-meta surfaces
- `/Users/roykid/work/molcrafts/molpy/src/molpy/` imports/exports that currently
  route through Frame/Block compatibility aliases
- exact cross-repo version metadata and release notes in all three repositories
- `.claude/notes/notes.md` after acceptance

## Tasks

- [x] Freeze a machine-readable list of old CPU Frame members, compatibility
      declarations, aliases, overloads, old versions, Python exports, and
      cross-repo fallbacks; make the scrub gate fail before implementation.
- [ ] Write compile-fail/import-fail/old-version tests proving the old C++,
      Python, metadata, and wire surfaces are absent—not merely discouraged.
- [x] Implement owning CPU Block/Column and typed-meta storage; delete
      `cpu::Atom` storage, named Frame members, proxy views, and Box shims.
- [x] Move all CPU production, app, Python, test, regression, example, and docs
      call sites atomically to typed handles and `simbox`.
- [x] Make schema traversal the only grow/copy/pack/materialize path and retain
      the explicit flat-topology/per-atom-topology boundary conversion.
- [x] Add native-dtype pack plans, collective canonical schema fingerprinting,
      deterministic runtime-column order, fail-loud registration limits, and
      `TagInt(2^53+1)` wire tests.
- [ ] Add storage-generation invalidation and ASan tests that grow capacity,
      re-stride topology, and insert a runtime column/meta entry between
      launches; stale raw handles must refresh or fail before dereference.
- [ ] Cut molrs to typed meta in core + every binding/serializer and delete the
      string-meta codec/fallback; release one exact new API version.
- [x] Delete molpy Frame/Block compatibility exports and migrate all in-tree
      consumers directly to the exact molrs API.
- [x] Delete the audited CPU compatibility-only aliases, solver wrappers,
      overloads, fallbacks, and behavior-preserving defaults; update callers to
      the one current abstraction.
- [x] Replace TaskSpec variant-order wire IDs with explicit IDs, bump the exact
      current version, and delete all old-version/order constraints and tests.
- [x] Remove CPU use and declarations of shared flat compatibility types,
      headers, constants, and log overloads; move any genuinely distinct
      current numerical constant under its owning method namespace.
- [ ] Run source scrub and the explicit CUDA-off CPU matrix: double-storage
      full domains with MPI off/on; frame/comm tests across supported
      StorageReal × WorkReal and smallbig/bigbig IDs; OpenMP required rather
      than silently optional; plus the supported mixed-storage core profile.
- [ ] Run public regressions/examples and hybrid `mpirun -n 2` tests with
      `OMP_NUM_THREADS=1,4`, molrs Rust/CXX/Python tests, and molpy tests against
      the exact three-repository revision set.
- [ ] Run separate resolution, allocation, memory, LJ reuse/reneighbor,
      materialize/write-back, pack/broadcast/halo/migration/gather, and SCF
      performance gates; record baselines and the breaking surface in
      `.claude/notes/notes.md` and changelogs.

## Testing strategy

- **Source/compile gates**: forbidden old declarations and call syntax are
  absent; small negative consumers fail to compile/import/decode.
- **Unit**: typed Block/meta dtype and ownership, grow/copy/pack, topology
  strides, generation invalidation under ASan, Box cache,
  wrong-key/type/version failures, and zero old aliases.
- **MPI**: schema-driven broadcast, partition, halo, migration, gather, and
  write-back on one and two ranks; native dtype fidelity, collective schema
  fingerprinting, hybrid OpenMP determinism, and no raw aggregate byte copy.
- **Regression**: the existing public CPU suite runs through the new physical
  storage without compatibility fixtures. Scientific reference values do not
  change merely to make the migration pass.
- **Cross-repo**: molrs core/CXX/Python and molpy full suites use the exact
  typed-meta release; old string meta and old imports are negative tests.
- **Performance**: independent resolution, allocation, LJ, lifecycle/comm,
  memory, and SCF block-build gates on a pinned controlled runner.

## Deferred — GPU

- Physical mapping ownership in `cuda::Frame` and `torch::Frame`.
- CUDA/Torch consumer migration, typed meta handles, and removal of their own
  compatibility declarations.
- CPU/CUDA/Torch scientific, memory, transfer, and performance parity.

These are explicit `status: deferred` acceptance work. This CPU link edits no
production file under `include/atomiverse/cuda`, `src/cuda`,
`include/atomiverse/torch`, or `src/torch`. Deferred backends do not justify a
CPU alias, mirror, fallback, conditional old path, or dual representation.
If a deferred backend still includes a deleted shared compatibility header, its
configuration remains explicitly deferred/non-accepted until that backend is
migrated; the CPU symbol is not restored to keep it compiling.

## Out of scope

- New simulation algorithms, force-field styles, or scientific models.
- Automatic checkpoint/schema migration tools or legacy-file converters.
- A second transition spec, deprecation release, compatibility branch, or
  runtime feature flag selecting the old representation.
