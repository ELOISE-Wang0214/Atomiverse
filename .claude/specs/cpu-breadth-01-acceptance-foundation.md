---
title: CPU breadth 01 — truthful acceptance, explicit component lifetime, and exact cross-repo contracts
slug: cpu-breadth-01-acceptance-foundation
status: approved
execution: in_progress
last_progress_update: 2026-07-16
science_required: true
created: 2026-07-15
acceptance: cpu-breadth-01-acceptance-foundation.acceptance.md
depends_on: frame-block-04-delete-compatibility
---

# CPU breadth 01 — truthful acceptance, explicit component lifetime, and exact cross-repo contracts

## Summary

This link makes later CPU probes capable of proving anything. It delivers a
non-vacuous CPU acceptance profile, pins and checks the Atomiverse/molrs CXX
contract, turns performance promises into versioned gates, makes component
borrowing explicit, and adds one out-of-library harmonic-oscillator consumer
that reuses the same user-defined `Energy|Force` evaluator across `Evaluate`,
`MD`, `Minimize`, and `Hessian`.

The probe is deliberately small. It does not add a new simulation method. It
asks whether the public architecture, build matrix, lifecycle, and tests can
observe failures before the remaining links depend on them.

## Domain basis

The conformance case is the isotropic harmonic oscillator

`U(r) = 1/2 k |r-r0|^2`, `F(r) = -k(r-r0)`.

It has exact energy, force, Hessian (`k I`), equilibrium (`r0`), and period
`2*pi*sqrt(m/k)`. Units are Atomiverse engineering units: `k` [eV/Å²], position
[Å], force [eV/Å], mass [amu], and time [fs]. The public regression cites the
classical harmonic-oscillator equations and records why each numerical
tolerance follows from the chosen timestep and finite-difference displacement.

Infrastructure parity is not a scientific oracle. The molrs contract is
accepted by version/ABI identity plus a real Frame round-trip and AM1-BCC
consumer call, never by freezing Atomiverse's own output as a golden value.

## Progress / current evidence (2026-07-16)

**Execution: in progress — 5/11 tasks; 3/11 active acceptance criteria
verified. GPU acceptance remains deferred.**

- **Public harmonic regression shipped** (`reg_harmonic_driver_matrix`):
  one user `Energy|Force` type + `LoopHandler` across Evaluate/MD/Minimize/
  Hessian; ac-005/ac-006 verified.

- Component lifetime is now explicit end to end: `add(rvalue/unique_ptr)` owns,
  `add_ref(lvalue)` borrows, and plain lvalue `add` is deleted on both Driver
  and PhaseBuilder. CPU categorical setters and in-tree callers use the same
  contract; `test_cpu_component_ownership` pins it.
- molrs exports an exact CXX API version/capability mask. Atomiverse checks the
  contract at configure and runtime and mismatch diagnostics name expected and
  found values. The sibling molrs tree is still dirty and unreleased, so no
  release tag/SHA or CI pin is claimed.
- The authoritative lifetime decision and typed Minimize/KSpace follow-on
  evidence are recorded in `.claude/notes/notes.md`.
- Still blocking closure: the public harmonic four-workflow regression,
  cross-repository Frame/Element/AM1-BCC consumer job, CPU core/full manifests,
  `cpu_acceptance` workflow, released molrs pin, measured lifecycle/allocation
  ledgers, versioned performance baselines, and the complete acceptance run.

## Architecture questions

1. Can one external evaluator type serve four CPU workflow shapes without a
   workflow-specific base class or private header?
2. Can component lifetime be read from the registration call rather than
   inferred from C++ value category?
3. Can the CPU build matrix fail when an enabled domain silently registers no
   tests, regressions, or sources?
4. Can a local molrs override remain convenient while a non-current CXX ABI
   fails loudly and CI consumes one exact release revision?
5. Can build, first-launch/rebuild, steady launch, and transfer costs be
   measured separately and compared with a versioned baseline?

## Design

### One authoritative CPU acceptance profile

Add an acceptance workflow/target separate from default `ctest`:

- `cpu-core`: double precision, MPI enabled, unit tests and public regressions;
- `cpu-full`: `ATV_BUILD_RHF=ON`, `ATV_BUILD_PW=ON`,
  `ATV_BUILD_GPW=ON`, `ATV_BUILD_KSPACE=ON`, Python/launcher, tests,
  regressions, and selected validation binaries;
- the existing mixed-storage job remains a compile/unit profile, not a claim
  that unsupported double-only domains were tested;
- both CPU profiles run for `dev` pushes and pull requests targeting `dev`;
- configuration emits and verifies a target/test manifest so a requested
  domain cannot produce a vacuous green build.

`tests/` remains seconds-fast and unit-only. The aggregate `cpu_acceptance`
target invokes public regressions and selected validations outside default
`ctest`.

The delivered clean-room entry is one copy-pasteable login-shell command:

```bash
bash -lc '
  module load buildenv-gcc/2026.03-mpich
  export OPENBLAS_ROOT="$OPENBLASDIR"
  export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu
  export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"
  cmake -S . -B build/cpu-accept -DCMAKE_BUILD_TYPE=Release \
    -DATV_BUILD_CUDA=OFF -DATV_BUILD_TORCH=OFF \
    -DATV_USE_MPI=ON -DATV_BUILD_TESTS=ON \
    -DATV_BUILD_EXAMPLES=ON -DATV_BUILD_REGRESSIONS=ON \
    -DATV_BUILD_RHF=ON -DATV_BUILD_PW=ON -DATV_BUILD_GPW=ON \
    -DATV_BUILD_KSPACE=ON -DATV_BUILD_PYTHON=ON \
    -DATV_BUILD_LAUNCHER=ON -DBLA_VENDOR=OpenBLAS \
    -DATV_MOLRS_PATH="$ATV_MOLRS_PATH"
  cmake --build build/cpu-accept -j$(nproc)
  tools/verify_cpu_test_manifest.py build/cpu-accept
  ctest --test-dir build/cpu-accept --output-on-failure
  cmake --build build/cpu-accept --target cpu_acceptance
'
```

The implementation may wrap this sequence in a checked script, but the script
must print and execute these stages rather than hide a smaller configuration.

### Cross-repo use is an exact-version contract

molrs-cxxapi exports a monotonically increasing `MOLRS_CXX_API_VERSION` and a
small capability set. Atomiverse records the sole expected API version and its
released molrs tag/SHA in one CMake-readable source of truth.
CI checks out that exact revision. `ATV_MOLRS_PATH` remains a developer
override, but the built consumer contract must reject a mismatched API before
simulation work begins.

The contract job compiles and runs three paths: Frame materialize/write-back,
element lookup, and AM1-BCC assignment. It must catch generated-bridge drift;
editing generated CXX output rather than its Rust bridge source is forbidden.
No current dirty sibling checkout is declared the contract baseline.

### Explicit component lifetime

Registration spelling carries ownership:

- `add(std::move(component))` / the owned-pointer overload transfers ownership;
- `add_ref(component)` is the only non-owning registration spelling;
- the same split exists on `PhaseBuilder`;
- implicit lvalue borrowing is removed or made a compile-time error, and all
  in-tree callers migrate intentionally;
- post-run introspection uses an explicit reference registration or a stable
  owned handle, never a pointer recovered from a temporary.

This link does not introduce `ComponentBase::clone()`. Replica instance scope
is deliberately left for link 05, after lifetime itself is unambiguous.

### Public conformance consumer

`regressions/architecture/harmonic_driver_matrix/` defines the evaluator in
the regression source using public headers only and links only `atomiverse`.
The same concrete evaluator instance/type is registered explicitly in:

- `cpu::Evaluate`: exact energy and force at a displaced point;
- `cpu::md::MD`: short NVE trajectory and the analytic period/energy bound;
- `cpu::minimizer::Minimize`: convergence to `r0` under an external handler;
- `cpu::vib::Hessian`: the three Cartesian Hessian eigenvalues equal `k`.

A user-defined `LoopHandler` stops one run by a metric/step rule. No library
source knows the handler or evaluator type. The case is one flat public
consumer and runs in less than 60 seconds.

### Lifecycle and performance ledger

Every later CPU scientific spec inherits six acceptance columns: scientific
oracle, public API regression, lifecycle (`build`/`launch`/`rebuild`), steady
allocation, transfer/synchronization budget, and versioned performance.

The baseline records setup/build, first launch, steady step, explicit rebuild,
and boundary materialization separately. Microbench regressions use a 5%
threshold and end-to-end cases a 10% threshold on repeated medians on the same
runner; a threshold crossing fails the acceptance job and publishes the raw
samples. Allocation-count tests run on the controlled acceptance runner rather
than being silently excluded everywhere.

## Files to create or modify

- `.github/workflows/ci.yml`
- `.github/workflows/acceptance-cpu.yml` (new)
- `cmake/MolrsContract.cmake` (new; exact version/revision source)
- `CMakeLists.txt`, `tests/CMakeLists.txt`, `regressions/CMakeLists.txt`
- `include/atomiverse/facade/driver.h`
- CPU driver registration call sites migrated to explicit lifetime
- `tests/cpu/driver/test_component_ownership.cpp` (new)
- configuration/profile closure tests under `tests/architecture/`
- `regressions/architecture/harmonic_driver_matrix/` (new, with `REFERENCE.md`)
- `benchmarks/REGRESSION.md` and versioned CPU baseline data
- molrs CXX bridge source + consumer-contract test needed for API versioning
- `.claude/notes/notes.md` after acceptance

## Tasks

- [x] Write failing component-lifetime tests for owned rvalues/pointers,
      explicit `add_ref`, and rejection of implicit lvalue borrowing on Driver
      and PhaseBuilder.
- [ ] Write the public harmonic driver-matrix regression before changing the
      registration API; prove that it needs no private/internal header.
- [x] Add molrs CXX API version/capability export and the exact Atomiverse
      build/runtime contract declaration with mismatch gates.
- [ ] Pin CI to the sole released molrs revision after these dirty cross-repo
      changes are released; never fabricate a tag/SHA from a worktree.
- [ ] Add Frame, element, and AM1-BCC cross-repo consumer-contract tests.
- [ ] Add CPU core/full configuration manifests and fail-loud dependency
      closure (`GPW => RHF && PW`, `KSPACE => FFTW`, Python/launcher => MPI and
      interpreter requirements as applicable).
- [ ] Add the `cpu_acceptance` aggregate target and dev/PR workflow; keep
      regressions and validations out of default ctest.
- [x] Implement explicit owned/reference registration and migrate all CPU
      callers without drive-by CUDA/Torch edits.
- [ ] Establish versioned build/steady/rebuild/materialization baselines and
      run allocation gates on the controlled CPU runner.
- [ ] Run CPU core/full builds, unit tests, regressions, selected validation,
      Python tests, launcher parity, and the performance comparison.
- [ ] Record measured counts, baselines, and the lifetime decision in
      `.claude/notes/notes.md`.

## Testing strategy

- **Unit (`tests/`)**: registration ownership, ABI mismatch failure, profile
  dependency closure, and target/test-count assertions. No E2E is added here.
- **Regression (`regressions/`)**: one public harmonic driver-matrix consumer,
  literature equations, ≤60 seconds.
- **Validation (`benchmarks/.../validation`)**: existing external-engine
  comparisons selected by the CPU full profile.
- **Performance (`benchmarks/.../performance`)**: repeated-median baseline,
  allocation ledger, and separated lifecycle phases; never default ctest.
- **Cross-repo**: exact released molrs revision, not the developer's dirty
  sibling tree; local override must pass the same ABI and runtime contract.

## Deferred — GPU

- CUDA/Torch conformance versions of the four-workflow harmonic case.
- Direct CUDA runtime synchronization/D2H/allocation interception.
- CUDA precision/configuration profile closure and performance baselines.

These appear as a deferred acceptance criterion and do not block this link.
No production file under `include/atomiverse/cuda`, `src/cuda`,
`include/atomiverse/torch`, or `src/torch` is modified.

## Out of scope

- New force-field styles, solvers, ensembles, or sampling algorithms.
- A universal serialized Plan/Component IR.
- `ComponentBase::clone()` or replica state policy.
- Fixing every existing scientific benchmark; this link makes their execution
  and result status visible.
- Pinning any current uncommitted molrs/molpy worktree state as a release.
