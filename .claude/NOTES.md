> **Moved to `.claude/notes/notes.md`.** This file kept for backward reference.
> New content is written to `.claude/notes/notes.md` via `/mol:note` or `/av-note`.

# Architecture Notes (legacy copy)

Accumulated architectural decisions and corrections from development iterations.
Agents read this alongside CLAUDE.md. When a note stabilizes, promote it into CLAUDE.md.

---

## [2026-06-06] Async analysis: freeze-then-dispatch invariant + device-ring (async-analysis-snapshot-ring)

**Invariant (freeze-then-dispatch).** Any snapshot delivered to an off-hot-path
consumer (recorder/observer) must be frozen into an immutable, step-stamped copy
**on the stream that mutates the frame**, before the next step can overwrite it;
consumers only ever read the frozen copy, never the live frame. Accuracy
(per-frame consistency + correct step stamp) is unconditional for every
component; completeness (may-drop-frames) is the orthogonal axis, governed by a
single `drop_oldest` policy. The MD hot loop never blocks.

**CPU (link 01).** `cpu::Publish` freezes the frame into a K-deep
`cpu::AnalysisExecutor` host ring with a synchronous OMP memcpy on the hot-path
thread (naturally ordered before the next step). Single uniform path:
`drop_oldest` only (block/grow + the per-component `backpressure()` override
removed), build-time static `mask_union_all`, global ring depth K
(`Config.ring_size`), shared worker pool.

**CUDA (link 02).** The fix lives in `cuda::AsyncBridge`, NOT a new
`DevicePublish` and NOT `cuda::Context.executor`: Thermo/RDF/Mulliken/Density
capture scalars / pair-distances / density-matrix, none of which fit the CPU
`FrameSnapshot` ring, so a shared-host-ring unification was the wrong shape.
Each `AsyncBridge` subclass owns a K-deep device ring + K pinned host slots
(K = compile-time `AsyncBridge::kRingDepth` = 3). `launch()` order:
`freeze_slot` (D2D live→device_ring[slot] on the **compute** stream) →
`cudaEventRecord(pack_events_[slot], main)` → `cudaStreamWaitEvent(defer)` →
`d2h_slot` (device_ring[slot]→pinned, defer stream) →
`cudaLaunchHostFunc(on_data_ready(slot))`. The compute-stream freeze is what
closes the prior cross-stream read-after-write tearing (the old path D2H'd from
live device pointers on the defer stream while the next step raced ahead).
`AsyncRecordBase` gains an internal `device_ring_[K]`+`host_ring_[K]` with
**slot-busy drop_oldest** (its HostWorker queue is unbounded, so it must drop a
capture whose slot is still in flight rather than overwrite it; `AsyncBridge`
needs no explicit drop because `cudaLaunchHostFunc` back-pressures the host
enqueue thread, keeping compute within K steps of defer).

**CUDA execution-model classification — add the device-ring-freeze row.**
Alongside *compute* (GPU-pure launch), *CPU bridge* (D2H + host callback), and
*host-dependent compute*: an `AsyncBridge`/`AsyncRecordBase` recorder/observer
is a **deferred device-ring-freeze CPU bridge** — `launch()` enqueues a
compute-stream D2D freeze + defer-stream D2H + host callback, zero
`cuda::stream_synchronize` on the step loop (guarded by
`tests/cuda/async/test_hot_path_no_sync.cu` and `test_async_bridge_tearing.cu`).

**Memory (link 03).** Ring footprint ≈ K·N·B (device ring) + equal pinned host
ring ⇒ ≈ 2·K·N·B; verified by `benchmarks/bench_async_snapshot_mem`.

---

## [2026-05-22] CUDA driver taxonomy migration — CUDA matches the three-category shape

Landed `.claude/specs/cuda-driver-taxonomy-*.md` (the 5-link
`cuda-driver-taxonomy` chain). Brings the CUDA backend onto the same
category-driven taxonomy CPU adopted in the 2026-05-06 redesign,
replacing the per-algorithm CUDA driver shape.

**The 5-link CUDA migration:**

  - **01** — namespace-normalized the CUDA backend (`ATV::cuda::*`) and
    added the `cuda::CUDADriver<PhaseEnum>` layer + the three category
    bases (`cuda::EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver`),
    parallel to `cpu::CPUDriver<PhaseEnum>`.
  - **02** — `cuda::md::MD` migrated to `cuda::DynamicsDriver`.
  - **03** — `cuda::minimizer::Minimize` (`OptimizeDriver`) with
    `cuda::minimizer::SteepestDescent` as a hot-swappable Optimizer
    component, replacing the `cuda::minimizer::SD`/`CG`
    algorithm-as-driver classes.
  - **04** — CUDA SCF demoted from a `Driver` to the
    `cuda::scf::SCFSolver` inner solver behind the `cuda::scf::SCFEvaluator`
    Evaluator component; the SCF step helpers relocated to
    `cuda::scf::internal::` (`include/atomiverse/cuda/scf/internal/`),
    `cuda::scf::scf.cuh` deleted. `cuda::Evaluate` is the one-pass driver;
    single-point CUDA energy is `cuda::Evaluate + cuda::scf::SCFEvaluator`
    (or the `cuda::single_point` helper) — no `cuda::scf::SinglePoint`
    class.
  - **05** (this link) — scrubbed the now-dangling CUDA references to
    `ATV::scf::Phase`: deleted the dead
    `Driver<scf::Phase, cuda::device_tag>` explicit instantiation in
    `src/cuda/driver/driver.cu` and re-pointed the surviving `@see` /
    usage-snippet doc-comments to the `cuda::scf::SCFSolver` API. A
    `tests/cuda/scf_phase_scrub.sh` grep gate (CTest
    `test_cuda_scf_phase_scrub`) now asserts no `scf::Phase` token
    survives under `include/atomiverse/cuda/` or `src/cuda/`. CLAUDE.md
    and this log updated to the shipped shape.

**Two latent bugs surfaced and fixed in link 01.** The namespace
normalization in 01 replaced the explicit `<float>`/`<double>` CUDA
kernel-launcher instantiations with `<work_t>`-only forced
instantiations. Under the default build config
(`ATV_CUDA_REAL=float`, `ATV_CUDA_WORK=double`) this exposed two
pre-existing precision-tier defects:

  1. **`real_t`/`work_t`/`prec_t` signature drift** (build regression,
     ~25 errors across 6 files). Neighbor-list/SCF kernels were already
     correct, but their launch wrappers and the shared
     `nlist_impl.cuh` declarations had drifted: per-atom SoA/AABB/
     displacement buffers are `real_t` (not `work_t`), cell/inverse-cell
     matrices are `real_t` (consumed by `triclinic.cuh`), and
     `frame.box` is `prec_t`. Launchers and header decls were
     corrected to flow each precision tier uniformly.
  2. **D2H copies sized by the wrong element type.** The RDF observer
     and the XYZ/density recorders copied per-atom / per-pair `real_t`
     device buffers (`frame.pairs.dist`, `frame.atoms.x/y/z`) sized by
     `sizeof(double)`, reading 2x past the end and leaving a sticky
     CUDA "invalid argument" error. Fixed by staging through
     `real_t` pinned buffers copied as `sizeof(real_t)`, then widening
     to `double` for the CPU consumers.

**`ATV::scf::Phase` kept — owned by the CPU plane-wave SCF.** The chain
plan assumed `scf::Phase` (`facade/stage.h`) would be unreferenced once
CUDA SCF became an inner solver and could be deleted. A tree-wide grep
during link 05 drafting showed that is false: the CPU **plane-wave**
SCF (`cpu::pw`) runs as a `Driver<scf::Phase, cpu::device_tag>`
(`src/cpu/driver/driver.cpp`) and seven `cpu::pw::PW_*` components return
`DefaultBinding` values keyed on `scf::Phase::*`. The enum and the CPU
instantiation are kept; only the enum's doc-comment was corrected from
the now-false "Legacy CUDA-only … will retire" to name `cpu::pw` as
its owner.

**Follow-up — migrate `cpu::pw` off `scf::Phase`.** The plane-wave SCF
is the last `Driver<scf::Phase>` outer loop. Moving it onto its own
phase enum / category driver is a genuine architectural refactor with
behavior risk, unrelated to the CUDA migration, and was left out of
this chain. Recorded here as the open follow-up.

## [2026-05-06] CPU driver taxonomy redesign — three categories + SCFSolver inner solver

Landed `.claude/specs/cpu-driver-taxonomy.md` (phases 1-6). Replaces the
single `Driver<PhaseEnum, DeviceTag>` template-with-algorithm-subclasses
scheme with a category-driven layer.

**Three CPU outer driver categories** (decided by termination semantic,
not by algorithm):

  - `cpu::EvaluateDriver` — one pass; concrete: `cpu::Evaluate`
  - `cpu::DynamicsDriver` — full step budget; concrete: `cpu::md::MD`
  - `cpu::OptimizeDriver` — early stop on convergence; concrete:
    `cpu::minimizer::Minimize` (renamed from `cpu::minimizer::LBFGS`)

The underlying `Driver<PhaseEnum, DeviceTag>` template is now a private
plan/dispatch primitive — user code subclasses one of the three
category bases via `cpu::CPUDriver<PhaseEnum>` (the CPU backend
layer, in `include/atomiverse/cpu/driver/cpu_driver.h`). Subclassing
`Driver<>` directly is forbidden in CPU code. CUDA migration will
add a parallel `cuda::CUDADriver<PhaseEnum>` sibling.

**SCF is not an outer category.** `cpu::scf::SCFSolver` is a public
low-level inner-solver class (NOT a Driver subclass), used internally
by `cpu::scf::SCFEvaluator`. Its 5 step helpers
(`HamiltonianBuilder`, `EigenSolver`, `DensityBuilder`, `DIISMixer`,
`DualConvergenceCheck`) live in `cpu::scf::internal::` and never
inherit `cpu::ComponentBase`. Tests may include the internal headers
directly; production facade headers must not. There is no
`SinglePoint` class — single-point energy is `cpu::Evaluate +
cpu::scf::SCFEvaluator(.with_gradient(false))`, plus the optional
`cpu::single_point(frame, evaluator)` helper.

**Algorithm-as-component.** LBFGS / CG / FIRE are Optimizer
components plugged into `Minimize`. NoseHoover / Langevin /
Berendsen are Thermostat / Barostat components plugged into `MD`.
The new `cpu::minimizer::CG` (Polak-Ribière+) is hot-swappable with
`cpu::minimizer::LBFGS` inside the same `Minimize` driver.

**Component tagging — two orthogonal axes**:
  - `role()` — structural slot (`Integrator`, `Optimizer`,
    `Evaluator`, `Thermostat`, `Barostat`, `ConvergenceCheck`,
    `Recorder`, `Observer`, `None`)
  - `capabilities()` — bitset of frame fields written
    (`Energy | Force | Stress | Charge | Dipole | DensityMatrix |
    Orbitals`)

`DynamicsDriver` validates "≥ 1 component with `Capability::Force`"
regardless of role — `FixLangevin` (Thermostat role, Force-capable)
satisfies a pure-Langevin run.

**Outer iteration status channel**:
`cpu::Context::outer_status` is a `{converged, stop_requested,
reason}` struct, owned by the category driver, reset before each
iteration. ConvergenceCheck-tagged components write it; SCF
internals never touch it.

**Refactor was full no-backwards-compat**: `cpu::scf::DIIS` →
`cpu::scf::internal::DIISMixer` (renamed in place; one test caller
migrated). `cpu::minimizer::LBFGSStep` → `cpu::minimizer::LBFGS`.
`cpu::minimizer::LBFGS` driver → `cpu::minimizer::Minimize`. All 6
in-tree call sites updated in one commit per phase. No
`[[deprecated]]` shims.

**CUDA migration deferred**: `cuda::md::MD`, `cuda::minimizer::CG`/
`SD`, `cuda::scf::SCF` keep the per-algorithm shape until
`.claude/specs/cuda-driver-taxonomy.md` lands.

**`scf_solve()`** survives as a ≤30-line shim that constructs a
transient `SCFSolver` and runs it. Decision (vs deletion): retained
because `rhf_solve` / `uhf_solve` / `uks_solve` and the gradient
finite-difference path still call it; revisit during CUDA SCF
migration.

**Phase 5b deferred**: IntegralPack/ShellList split (drop `shells` +
`shell_data` from `IntegralPack` into a separate `ShellList` view)
and the 5 large file splits (`native_integrals.cpp` 1710 lines etc.)
remain TODO. Mechanical work; better as dedicated commits.

**Rule promoted to CLAUDE.md**: the three-category taxonomy, the
`role()`/`capabilities()` axes, and the "SCF is an inner solver, not
an outer category" decision.

## [2026-03-25] DeviceSystem removed — metadata lives in cuda::Frame

DeviceSystem class was deleted. Host-cached metadata (n_atoms, n_bonds, n_angles, n_dihedrals) now lives directly in `cuda::Frame` as plain Int fields with accessor methods. The device-side atom count pointer is `d_n_atoms` (not `n_atoms`). (Note: `DeviceFrame` was later renamed to `cuda::Frame`.)

**Rule**: Contexts (BuildContext, LaunchContext, DestroyContext) hold `cuda::Frame& frame` directly. Components access `ctx.frame.n_atoms()` for host-cached count and `ctx.frame.fx` etc. for device pointers. No wrapper class around `cuda::Frame`.

## [2026-04-06] CPU SCF force evaluators → Component refactoring

RHFForce, UHFForce, KSForce bypass the Driver/Component model — each owns a velocity Verlet loop (170 lines duplicated 3x), ASPC predictor logic (140 lines duplicated), and gradient skeleton (90 lines duplicated). Total ~590 lines of unnecessary duplication.

**Completed refactoring:**
1. ✓ Deduplicated gradient.cpp: extracted `accumulate_W`, `accumulate_common_gradient`, `make_result` helpers; unified FD gradient into single callback-based function
2. ✓ Unified `scf_gradient()` replacing separate rhf_gradient/uhf_gradient/ks_gradient — dispatches HF vs KS based on SCFConfig. Key subtlety: restricted HF stores P_alpha = P_total (factor 2 included), unrestricted stores per-spin densities
3. ✓ Created `SCFEvaluator` (cpu::ComponentBase at md::Phase::force) — single component for all methods, with optional ASPC predictor
4. ✓ Deleted legacy `RHFForce`, `UHFForce`, `KSForce` standalone classes
5. ✓ ASPC predictor lives inside SCFEvaluator (no separate wrapper needed)

**Rule promoted to CLAUDE.md:** "Driver-owns-the-loop" and "No duplication across methods" rules added to Component model section.

**Systematic bugs found during L-BFGS implementation:**
1. CUDA Driver classes use `Driver<Phase>` (1 arg) but template requires 2 `<PhaseEnum, DeviceTag>`. Latent bug — compiles only because CUDA sources not compiled in CPU-only builds. Fix deferred to CUDA session.
2. Infrastructure components (`md::ZeroForces`) hardcode `md::Phase::pre_force` in default_bindings — can't be reused by minimizer. Fix: create `cpu::minimizer::ZeroForces`/`ZeroEnergy` with bindings to `minimizer::Phase::force`.
3. CPU Driver subclass headers MUST include `cpu/backend_traits.h` for the `backend_traits<cpu::device_tag>` specialization. Without it, `this->add()` / `this->on()` fail at template instantiation.

## [2026-04-04] molrs integration complete

> **Superseded [2026-07-02]:** `cpu::MolRec` (and `io/molrec.{h,cpp}`) is
> deleted — the bridge is now raw-array FFI wrappers (`cpu/io/frame_writer.h`
> per-frame writers; `cpu::compute::*` batch computes + streaming accumulators
> in `cpu/observer/trajectory_metrics.h`, PIMPL). See the 2026-07-02 entry.

molrs CXX FFI integration is done. Key decisions that diverged from original spec:
- No `MolFrame` type in C++ -- all frame data lives in `molrs::MolRec` on the Rust side
- `cpu::MolRec` uses PIMPL (not direct `rust::Box`), keeping `rust/cxx.h` out of public headers
- Bridge is precision-templated (`bridge.rs` generated from template via `build.rs`)
- `ZarrWriter` delegates to `molrs::MolRec::write_zarr` (no separate `zarrs_ffi` dependency)
- Mulliken analysis stays in C++ (xtensor-based), not migrated to molrs
- `ATV_HAS_ZARR` / `ATV_BUILD_ZARR` flags removed -- Zarr is always available via molrs

## [2026-04-06] Test directory mirrors src/ layout

Tests are structured as `tests/{cpu,cuda,torch,facade}/...` mirroring the `src/` directory layout. Each test file uses GoogleTest fixtures (`TEST_F`) to organize tests by component/driver.

**Directory mapping:**
- `tests/cpu/md/` → CPU MD components (VelocityVerlet, LJ, CellList, bonded)
- `tests/cpu/scf/` → CPU SCF (integrals, solvers, gradient, grid, XC, evaluator)
- `tests/cpu/minimizer/` → CPU minimizer (L-BFGS)
- `tests/cuda/{md,potentials,minimizer,scf}/` → CUDA backend tests
- `tests/facade/` → Plan builder tests
- `tests/torch/` → LibTorch backend tests

**CMake naming convention:** `test_{backend}_{domain}_{component}`, e.g. `test_cpu_scf_gradient`, `test_cuda_md_thermostat`. Helper functions `atv_add_test(name, source_path)` and `atv_add_cpp_test(name, source_path)` take explicit source paths.

**Rule**: New test files go in the directory matching their source module. Use `TEST_F` fixtures to group tests by component. Each fixture class = one component or driver.
