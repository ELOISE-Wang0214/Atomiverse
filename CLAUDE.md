---
mol_project:
  name: atomiverse
  language: cpp
  build:
    install: "cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON"
    check: "cmake --build build -j$(nproc)"
    test: "ctest --test-dir build --output-on-failure"
    test_single: "ctest --test-dir build -R {path} --output-on-failure"
    coverage: "ctest --test-dir build --output-on-failure"
  arch:
    style: backend-pillars
    rules_section: "## Architecture"
  doc:
    style: doxygen
  science:
    required: true
  ci:
    config: .github/workflows/ci.yml
  notes_path: .claude/notes/notes.md
  specs_path: .claude/specs/
---

# CLAUDE.md

Atomiverse is a dual-backend (CPU + CUDA) molecular-simulation engine with one
shared Driver/Component API spanning classical molecular dynamics, Monte Carlo
sampling, geometry optimization, vibrational analysis, and electronic
structure (Gaussian-basis SCF, plane-wave DFT, Gaussian-and-plane-waves DFT,
and NDDO/AM1 semiempirics). C++17, `ATV::` namespace.

> **Branch note:** active development is on **`dev`**. `master` is a single
> squashed snapshot with the *old* pre-2026 layout (`cpu/frame.h`, `cpu/md/lj.h`,
> no `pw/`/`potentials/`); do not read it for current behaviour. New git
> worktrees default-base from `origin/master`, so branch a worktree off `dev`
> explicitly (`git worktree add -b <name> <path> dev`).

## Build & Test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
./build/tests/test_md             # single test binary
```

Requires: C++17, CMake ≥ 3.18, a Rust toolchain (for molrs), xtensor
(auto-fetched). CUDA toolkit only for the CUDA backend; FFTW3 only for
plane-wave / k-space.

### Toolchain via Lmod modules — never hardcode `/software/...` paths

Run everything in a **login shell** (`bash -lc '...'`) so Lmod's `module`
function is defined. The build env does **not** auto-load — load it explicitly:

- **CPU (x86_64 epyc login node):** `module load buildenv-gcc/2026.03-mpich`
  (GCC 14.3 + OpenBLAS/LAPACK/ScaLAPACK + FFTW + MPICH). The older
  `buildenv-gcc/2026.03` is **deprecated** — use the `-mpich` variant.
- **GPU (aarch64 GH200 node):** `module load GPU/buildenv-gcccuda/2026.03-cu13.0`
  (GCC 14.3, CUDA 13.0, GPU-aware MPICH, NCCL, cuFFT, FFTW, OpenBLAS).

A correct CPU configure+build+test from an automated agent, copy-pasteable:

```bash
bash -lc '
  module load buildenv-gcc/2026.03-mpich
  export OPENBLAS_ROOT="$OPENBLASDIR"                       # see BLAS note below
  export RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu   # see Rust note below
  export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"   # see runtime note below
  cmake -B build -DCMAKE_BUILD_TYPE=Release \
        -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON -DATV_BUILD_PW=ON \
        -DATV_BUILD_CUDA=OFF -DBLA_VENDOR=OpenBLAS \
        -DATV_MOLRS_PATH=/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/molrs
  cmake --build build -j$(nproc)
  ctest --test-dir build --output-on-failure
'
```

**BLAS:** the module exports `OPENBLASDIR` and adds OpenBLAS to
`CMAKE_PREFIX_PATH` — it does **not** set `OPENBLAS_ROOT`. Mirror it
(`export OPENBLAS_ROOT="$OPENBLASDIR"`) and pass `-DBLA_VENDOR=OpenBLAS`.

**Rust:** on the x86_64 login node the *default* rustup toolchain is wrongly
`aarch64`, so corrosion's `cargo metadata` aborts ("no such file or
directory"). Force `RUSTUP_TOOLCHAIN=stable-x86_64-unknown-linux-gnu`. On the
GPU node the default (aarch64) is correct — do **not** set the override there.

**molrs path:** molrs is a hard dependency built via corrosion. The CMake
default `../molrs` breaks from a worktree, and `$HOME/work/...` is wrong when
the tree lives under `/nobackup`. Always pass an absolute
`-DATV_MOLRS_PATH=<abs path to molrs>` (here:
`/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts/molrs`).

**Runtime gotcha (#1 cause of `libopenblas.so.0` / `GLIBCXX_3.4.32 not found`):**
the module sets `LIBRARY_PATH` (link time) but not `LD_LIBRARY_PATH` (run
time). Mirror it before `ctest` or before running any binary:
`export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"`.

### GPU node (CUDA build & test on Arrhenius)

The **login nodes are x86_64 with no CUDA toolkit** — the CUDA backend
(`ATV_BUILD_CUDA=ON`, any `.cu`) builds and runs only on the **aarch64
Grace/GH200 GPU node**. The CPU backend (incl. RHF and plane-wave) builds and
runs fine on the login node with `-DATV_BUILD_CUDA=OFF`.

- **Interactive:** `hpc gpu` (NAISS `interactive` srun wrapper; defaults from
  `~/.hpc/config.sh`: partition `gpu`, account `naiss2026-4-715-gpu`, 1 GPU, 3 h).
- **Batch** (use from an automated agent — `hpc gpu` needs a TTY):
  ```bash
  sbatch -p gpu -A naiss2026-4-715-gpu --gpus 1 -t 01:00:00 -c 8 -o out.log job.sh
  ```
- **Inside the job script** (`#!/bin/bash -l` so Lmod is defined; never pipe
  `module load` — a subshell loses the env):
  ```bash
  #!/bin/bash -l
  module load GPU/buildenv-gcccuda/2026.03-cu13.0
  cmake -B build-cuda -DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON \
        -DATV_MOLRS_PATH=$HOME/work/molcrafts/molrs   # adjust to the real path
  cmake --build build-cuda -j$(nproc)
  ctest --test-dir build-cuda --output-on-failure
  ```

### CMake options

| Option | Default | Enables / sets |
|---|---|---|
| `ATV_BUILD_CUDA` | ON | CUDA backend; defines `ATV_HAS_CUDA` |
| `ATV_BUILD_RHF` | OFF | Gaussian SCF (RHF/UHF/UKS); defines `ATV_HAS_RHF` |
| `ATV_BUILD_PW` | OFF | Plane-wave DFT (needs FFTW3) |
| `ATV_BUILD_GPW` | OFF | Gaussian-and-plane-waves DFT (needs `ATV_BUILD_RHF` + `ATV_BUILD_PW`) |
| `ATV_BUILD_KSPACE` | OFF | Reciprocal-space electrostatics PME (needs FFTW3) |
| `ATV_BUILD_TORCH` | OFF | LibTorch backend; defines `ATV_HAS_TORCH` |
| `ATV_BUILD_TESTS` / `ATV_BUILD_BENCHMARKS` / `ATV_BUILD_EXAMPLES` | OFF | test / benchmark / example targets |
| `ATV_USE_MATPLOT` | OFF | plot generation in selected examples (needs gnuplot) |
| `ATV_REAL_STORAGE` / `ATV_REAL_WORK` | `double` / `double` | LAMMPS-style real tiers (`float`\|`double`): per-atom storage / kernel-work precision. Accumulators (`AccReal`) are always `double`. |
| `ATV_ID_SIZE` | `bigbig` | integer-ID width: `bigbig` (64-bit `TagInt`, default — preserves legacy int64 tags) or `smallbig` (32-bit, LAMMPS/molrs-aligned memory opt-in). `BigInt` global counts always 64-bit. |

Performance: `-march=native` is required for the CPU SCF ERI kernels.

## Repository layout

```
include/atomiverse/ + src/   mirrored header/impl trees (architecture below)
tests/                       unit tests ONLY — each exercises one function or
                             method, and the WHOLE SUITE runs in seconds (it is
                             the fast-iteration loop). E2E IS STRICTLY
                             FORBIDDEN here: no full SCF/MD/vib solve, no MC
                             equilibration to convergence, no timing, no
                             external-engine comparison, no third-party package
                             read at build or test time. There is no opt-in
                             switch to run e2e from ctest: if it is registered,
                             it runs.
regressions/                 <domain>/<case>/ — the ONLY home of end-to-end.
                             One minimal textbook case per delivered feature,
                             written as an OUT-OF-LIBRARY consumer (public
                             headers + link atomiverse, nothing else), checked
                             against a cited LITERATURE value, <= 60 s. Built
                             behind ATV_BUILD_REGRESSIONS; run via the
                             `regressions` target, never by default ctest.
                             Every feature ships one. See regressions/README.md.
benchmarks/                  <domain>/<case>/{validation,performance} — see benchmarks/README.md
                             Parity against an external ENGINE (LAMMPS/PySCF/
                             CP2K) + speed/scaling. Run at acceptance time,
                             never by ctest.
examples/                    flat numbered tutorial sequence — see examples/README.md
docs/                        narrative docs (user-guide / tutorial / developer / api)
tools/                       data generators + debug harnesses
python/                      scikit-build-core + nanobind wheel (drives root CMake)
cmake/                       CMake helper modules
```

## Architecture

Three pillars: **Driver** (owns the loop), **Component** (one single-step
operation), **Frame** (device-resident state). Three execution backends
(`cpu/`, `cuda/`, `torch/`) sit behind one CUDA-free facade.

### Backend directory layering

Each backend root contains **directories only — no loose files**. Every
directory is a pillar layer, an infrastructure cluster, or a method domain:

```
<backend>/
  frame/       Frame pillar: device-resident state + residency ops
  driver/      Driver pillar: the loop node (CPUDriver), LoopHandler, OutcomePolicy,
               SubDriver (nesting), concrete Evaluate, single_point,
               step executor, backend_traits, step-loop instrumentation
  component/   Component pillar: ComponentBase, Launch/Context, generic + shared bases
  io/          external formats / molrs bridge (extxyz, frame_writer, lammps_data; CPU only)
  device/      CUDA/torch runtime utilities (device_manager, error, launch, timer)
  comm/        infra: domain decomposition, halo exchange, migration, reductions
  neighbor/    infra: cell/Verlet neighbor structures + the published PairList
  async/       infra: async analysis pipeline (snapshot ring, D2H bridge, executors)
  linalg/      infra: shared dense linear algebra (eigensolve, density build; CPU only)
  md/          domain: integrators, thermostats, barostats, deform, velocity init
  mc/          domain: Monte Carlo — the MC preset, movers,
               trial transaction (CPU only)
  minimizer/   domain: SD, CG, L-BFGS + opt recorders (LbfgsHistory: the curvature
               history + two-loop recursion, over a plain vector — shared with path/)
  replica/     infra: ReplicaSet + ReplicaExpansion (one Frame -> N coupled replicas)
  path/        domain: coupled-replica methods — CI-NEB (tangent/spring/climb as
               `couple`-slot components); PIMD/REMD would live here too
  scf/         domain: Gaussian-basis SCF (solver, integrals, gradients, evaluator)
  semiempirical/ domain: NDDO/AM1 semiempirical SCF (NDDOSolver, AM1Evaluator,
               AM1-BCC charge assigner; CPU only)
  pw/          domain: plane-wave DFT (PWSolver, PWEvaluator)
  gpw/         domain: Gaussian-and-plane-waves DFT (GPWSolver, GPWEvaluator; CPU only)
  xc/          domain: shared exchange-correlation pieces (Padé LDA parameterization)
  vib/         domain: finite-difference Hessian, thermochemistry (CPU only)
  kspace/      domain: reciprocal-space electrostatics (Ewald, PME)
  potentials/  domain: real-space force-field evaluators (pair + bonded styles)
  observer/    domain: derived quantities (RDF, virial stress, Mulliken, density)
  recorder/    domain: I/O (XYZ, cube, thermo, timing report)
```

**Rules:**
- **No loose files at a backend root.** Shared code goes in the matching pillar
  dir (`frame/`, `driver/`, `component/`, `io/`, `device/`); a new cohesive
  multi-file infra cluster gets its own root dir (like `comm/`, `neighbor/`,
  `async/`). Never park shared code inside a domain subdir.
- **Orthogonality: the same feature lives in the same directory on every
  backend.** Virial stress is `observer/`, neighbor structures `neighbor/`,
  decomposition `comm/` — on CPU and CUDA alike.
- Electrostatics splits LAMMPS-style: real-space pair styles in `potentials/`
  (`coul_long`), reciprocal-space solvers in `kspace/`.
- Namespaces are not yet aligned to directories (`cpu::md::mpi::Reduce` lives in
  `cpu/comm/`, `cpu::Frame` in `cpu/frame/`); aligning them is a deliberate
  follow-up, not a drive-by.

**Box** (`cpu::Box` in `cpu/frame/box.h`, `cuda::Box` in `cuda/frame/box.cuh`)
is the canonical cell type in each backend's `frame/` layer: owns the cell
matrix plus a per-step-immutable cache (`inv_cell`, `inv_L`, `half_L`, `ortho`)
and exposes `length(d)`, `volume()`, `min_image(...)`, `wrap_axis(...)`,
`scale_*`, `set_cell`. CPU adds an xsimd `min_image_batch<B>` in
`cpu/frame/box_simd.h`. LJ keeps an inlined branchy min-image in its xsimd
ket-loop because member-load-through-`this` defeats hoisting there (~7%/step).

### Facade (CUDA-free)

```
core/        Real, Int, Mat3, units, constants, device tags, backend_traits
facade/      Driver<PhaseEnum, DeviceTag>, Plan, Stage (phase enums)
potentials/  backend-agnostic pair force/energy kernels shared by CPU/CUDA
molrs.hpp    stable CXX FFI header
```

### Runtime / Driver / Handler — there is no Driver "category"

**A Driver is a staged loop. Termination is a Handler.** There are no category base
classes: `EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver` / `ResponseDriver` /
`SamplingDriver` are **deleted** (CPU and CUDA alike). They were six classes that
differed from one another by one `if`, every new method demanded a seventh, and no user
could write an eighth — "when do we stop" was a private `if` inside a driver class.

Stripped to the bone a Driver is three things: **stages** (a typed phase enum),
**iteration** (walk them repeatedly), **termination** (when to stop). Only the third was
welded into the hierarchy. Hand it to a policy object and no category remains:

| | what it is | how many |
|---|---|---|
| **Runtime** (`cpu::Backend`) | frame binding, node storage, timing, MPI, domain decomposition | **one per run**, shared by the whole driver tree |
| **Driver** (`cpu::CPUDriver<P>` / `cuda::CUDADriver<P>`) | **one loop node**: stages + plan + handlers | any number, **nestable** |
| **Handler** (`cpu::LoopHandler`) | asked once per iteration | any number, composable, **user-definable** |

```
Driver<PhaseEnum, DeviceTag>      plan / dispatch primitive
  └─ cpu::CPUDriver<PhaseEnum>    THE loop node
       └─ {MD, Minimize, MC, Evaluate, Hessian, PathDriver}    presets
```

Two levels, not four. A **preset** is *not* an architectural layer: it is a class whose
**constructor** attaches the default handler + infrastructure a given kind of run wants.
It contributes defaults, not semantics, and the user overrides any of them. (A class with
a constructor — never a `make_md()` factory; see the OOP rule below.)

| Preset | Phase enum | Default handler | Termination |
|---|---|---|---|
| `cpu::Evaluate` / `cuda::Evaluate` | `evaluate::Phase` | — | one pass |
| `cpu::md::MD` / `cuda::md::MD` | `md::Phase` | — | the step budget |
| `cpu::minimizer::Minimize` / `cuda::…::Minimize` | `minimizer::Phase` | `ForceConverged` | convergence, or the budget |
| `cpu::vib::Hessian` | `vib::Phase` | — | fixed 6N displacements + finalize |
| `cpu::mc::MC` | `sampling::Phase` | — (`Metropolis` **outcome**) | the trial budget |
| `cpu::PathDriver` | `path::Phase` | `AllReplicasConverged` | every replica converged |

**Two orthogonal policy slots.** Every loop has both:

- **`LoopHandler`** — *should the loop stop?* `add_handler(…)`, **plural**, composed by
  severity: **Stop > Converged > Continue**. Deliberately **not** a template: it sees a
  `LoopState` and *never a phase*, which is exactly what lets ONE handler serve MD,
  Minimize, MC and Path alike. Ships `ForceConverged`, `AllReplicasConverged`,
  `StopOnNaN`; a user writes `StopWhenTooHot{2000}` (LAMMPS `fix halt`) from **outside**
  the library. Its verdict is all-reduced across MPI ranks **by the driver**, once —
  handler authors never think about it.
- **`OutcomePolicy`** — *does this iteration count?* `set_outcome(…)`, **singular**,
  `Commit`/`Rollback`. MC's accept/reject was never a termination: rejecting a trial move
  discards one iteration's work and keeps going. `Metropolis` is one implementation. The
  policy **decides**; whoever recorded the changeset **undoes**.

**Algorithms never judge themselves.** `LBFGS`/`CG` **publish metrics** (`max_force`,
`rms_force`, displacement, ΔE) into `LoopMetrics`; a **handler judges** them.
`LBFGS::converged()` is deleted — `Driver::converged()` replaces it, because the *loop*
knows and the algorithm does not. Changing the convergence criteria does not touch the
optimizer. A component may **request** a stop; it **cannot declare convergence** — there
is no field left for it to write.

**Status belongs to the loop, not the Context.** Each loop node **owns** its
`IterationStatus` + `LoopMetrics`; `Context::loop` (a `LoopServices` slice) only *borrows*
the active loop's. `run_nested` re-points them at the child and restores the parent's on
exit, so **an inner loop cannot stop its enclosing loop** — a pointer that points at the
current loop cannot reach the loop above it. This bug shipped twice (LBFGS's unconditional
`converged` write; the first converged NEB image stopping the whole band) and is now
unwritable rather than merely forbidden. Pinned by `tests/cpu/driver/test_sub_driver.cpp`.

**A Driver is an executable node.** Same lifecycle as a component (build → launch × N →
destroy), so a Plan slot holds a component **or** a sub-driver (`cpu::SubDriver`). The
child runs its own Plan on the **parent's Runtime** (`ctx.exec.runtime`): no second
Backend, no rebuilt decomposition, no re-bound frame, no split communicator. **Type
erasure happens at the NODE boundary, never at the stage boundary** — the child keeps its
own phase enum, and stages stay typed enums all the way down (see *Stage typing*, below).

**Replica count was never an axis of a taxonomy — it is a nested loop.** So REMD/PIMD need
**no new category**: they are "a replica loop + a budget handler + a different `couple`
component". The replica machinery is infrastructure (`cpu/replica/`: `ReplicaSet`,
`ReplicaExpansion`), never inside a driver.

**Replicas are generated by an expansion strategy, never supplied by the user.**
`ReplicaSet` is built from `(Frame, ReplicaExpansion&)`, and `run()` takes ONE Frame and
writes ONE Frame back (CI-NEB: the climbing image = the saddle), so Frame stays the sole
I/O currency. Per-replica parameters (temperature ladder, bias centre, bead index, path
coordinate) belong to the *expansion* — putting them on the driver is what welded a
category to one method.

**Known ceiling: components cannot be cloned.** Running one inner Driver body over N
replicas would need N independent component states (N L-BFGS curvature histories), and
there is no `clone()` on `ComponentBase`. Nesting makes this **sharper**, not softer — a
`SubDriver` holds ONE inner driver with ONE component set. Undecided; see `.claude/NOTES.md`.

**Adding a new preset is cheap and adding a new "category" is meaningless.** If you catch
yourself wanting a new driver class, ask what *terminates* differently — and write a
`LoopHandler` instead.

### Stage typing — `slot_of` is the only cast

The six Phase enums all start at 0, so an `int` phase **cannot** be told apart across
workflows (`stage == 1` is simultaneously `md::pre_force`, `minimizer::line_search` and
`evaluate::evaluate`). `Backend::run_phase(int)` is gone. Everything a human writes goes
through a typed entry point:

| Use | Typed entry point |
|---|---|
| the ONE enum→slot conversion | `ATV::slot_of(phase)` (`constexpr`, zero cost) |
| the `COUNT` sentinel as a slot count | `ATV::slot_count<PhaseEnum>()` |
| replay one phase | `backend.run_stage(phase, step)` |
| "is this my stage?" (component) | `ctx.at(phase)` — cpu **and** cuda `LaunchContext` |
| declare a default binding | `ATV::make_binding(phase, order)` |

`test_stage_typing_scrub` (pure grep, every configuration) fails the build on any
hand-written `static_cast<int>` over a Phase enum anywhere in the tree.

`DefaultBinding` still carries `{PhaseKind, int phase, int order}` — that pair **is** the
erasure, not a leak: `ComponentBase` is not a template, so `default_bindings()` must
return one concrete type for every workflow. The `PhaseKind` beside the int is what makes
it safe, `make_binding` is the only way to build the pair, the cross-workflow misbind gate
throws at registration, and no caller ever types the int.

**Methods are never categories.** A *method* — how energy/forces are computed
(Gaussian SCF, plane-wave DFT, GPW, NDDO/AM1) — enters the architecture as an
**inner Solver** (a public low-level class that owns its iteration loop and
workspace; **not** a `Driver` subclass) wrapped by an **Evaluator component**
that any driver binds at its force/energy phase:

| Method | Inner solver | Evaluator component |
|---|---|---|
| Gaussian SCF (RHF/UHF/UKS) | `cpu::scf::SCFSolver` / `cuda::scf::SCFSolver` | `cpu::scf::SCFEvaluator` / `cuda::scf::SCFEvaluator` |
| Plane-wave DFT | `cpu::pw::PWSolver` | `cpu::pw::PWEvaluator` |
| Gaussian-and-plane-waves DFT | `cpu::gpw::GPWSolver` | `cpu::gpw::GPWEvaluator` |
| NDDO/AM1 semiempirics | `cpu::semiempirical::NDDOSolver` (AM1 parameterization) | `cpu::semiempirical::AM1Evaluator` |

*Properties* (charges, dipoles, Mulliken populations, AM1-BCC charges, …) are
component **Capabilities** — frame fields an Evaluator or Observer advertises
and writes — and never appear at the Driver level: there is no "charge
driver", only components with `Capability::Charge`.

Single-point energy in user code is `Evaluate + <method>Evaluator` (or the
`single_point` helper). There is **no `SinglePoint` class** and **no SCF phase
enum** — `ATV::scf::Phase` is deleted and grep-gated (`test_scf_phase_scrub`,
`test_cuda_scf_phase_scrub`). The `scf_solve()` / `pw_scf()` / `gpw_scf()`
shims remain as thin functional wrappers over a transient solver. Solver step
helpers (`*::scf::internal::*`, `*::pw::internal::*`) are owned by the solver,
never inherit `ComponentBase`, and are never registered into a Plan;
facade-public headers must not include them (tests may).

```
Preset::run(frame, budget)              # the preset composes the loop primitives
  setup_pass(frame)                     # binds ctx.loop at THIS loop's status+metrics
  for each iteration:
    status.reset()
    backend.step(idx)                   # components measure into LoopMetrics
    outcome = outcome_policy?           # Commit | Rollback  (MC only)
    verdict = poll_handlers(idx)        # Stop > Converged > Continue, MPI-allreduced
    if verdict != Continue: break
  teardown_pass()

Driver::run_nested(ctx, frame, budget)  # the SAME loop, on a BORROWED Runtime
  save ctx.loop; point it at MY status+metrics
  ... loop, via runtime->step_plan(my_plan, idx) ...
  restore ctx.loop                      # the child's verdict stays on the child
```

User code does not subclass `Driver<>` directly — it uses a **preset**, or subclasses the
loop node (`CPUDriver<PhaseEnum>`) when it needs a genuinely new loop shape. A preset adds
default handlers + infrastructure in its **constructor**; no runtime behaviour beyond the
loop node.

### Component tagging

Two orthogonal axes (`facade/component_role.h`):
- `ComponentRole role()` — structural slot: `Evaluator`, `Integrator`,
  `Optimizer`, `Thermostat`, `Barostat`, `ConvergenceCheck`, `Recorder`,
  `Observer`, `Mover`, `None`.
- `Capability capabilities()` — bitset of frame fields `launch()` may write:
  `Energy`, `Force`, `Stress`, `Charge`, `Dipole`, `DensityMatrix`,
  `Orbitals`. Bitwise `|`/`&`/`has_capability`.

**Algorithms are components, not driver subclasses.** `LBFGS`/`CG`/future
`FIRE` are Optimizer components plugged into `Minimize`;
`VelocityVerlet`/`NoseHoover`/`FixLangevin`/`Barostat` are
Integrator/Thermostat/Barostat components plugged into `MD`;
`TranslateMove` is a Mover component plugged into `MC`.

**PhaseKind-typed bindings** (`core/binding.h`, `facade/stage.h`): a
`DefaultBinding` is `{PhaseKind kind, int phase, int order}` constructed via
`make_binding(Phase, order)` — the `kind` discriminant tags which workflow's
phase enum the raw `int` came from. Driver registration validates every
binding's kind against the driver's own phase domain and throws naming both
workflows (no silent cross-workflow misbinds). Components used in more than
one workflow (NeighborList, pair/bonded styles, Thermo, …) must **not**
declare `default_bindings()` — they are registered explicitly via
`.on(phase).add(comp, order)` or a categorical setter, which binds them in
the calling driver's own phase domain.

**Opt-in protocols** — `ComponentBase` stays minimal; cross-cutting
infrastructure needs are declared by implementing a protocol the engine
probes for (no `dynamic_cast`, no dead default virtuals):

| Protocol | Declares | Engine probe |
|---|---|---|
| `HaloParticipant` (`component/halo_participant.h`) | declarative halo `CommField`s + endpoints; engine-side packing in `comm::HaloSend/HaloRecv` | `halo_participant()` |
| `NeighborConsumer` (`component/neighbor_protocol.h`) | `neighbor_cutoff()` the neighbor stage must cover | `as_neighbor_consumer()` |
| `NeighborProvider` (same header) | this component builds/publishes the pair list | `as_neighbor_provider()` |
| `PairEvaluator` (`component/pair_evaluator.h`) | shared non-bonded pair-style base (cutoff, special-bonds) | — |
| `AsyncComponentBase` (`async/async_component_base.h`; CUDA analogue `AsyncBridge`) | deferred snapshot-ring analysis off the hot path | `scheduling() == deferrable` |

**Iteration status channel** (`facade/iteration_status.h`):
`IterationStatus` is `{converged, stop_requested, reason}`, **owned by the loop node**
(`CPUDriver::status_`) and reset before each iteration. The Context only *borrows* the
active loop's, through the non-owning `ctx.loop` (`LoopServices{status, metrics}`) — so a
nested loop physically cannot reach its parent's. A component **measures** into
`ctx.loop.metrics` and may **request** a stop; only a `LoopHandler` declares convergence,
and only the driver writes it. SCF/PW internals never touch any of it.

**Rules:**
- **Driver owns the loop.** Evaluators, integrators, thermostats are
  Components — never standalone classes with their own loops.
- **Termination is a Handler, never a class.** Want a new stopping rule? Write a
  `LoopHandler`. A new driver class is almost never the answer.
- **Capability over role for evaluator validation.** `Evaluate` accepts any
  `role()==Evaluator` with non-empty capabilities; `MD` requires any component with
  `Capability::Force`; `Minimize` requires an Energy/Force evaluator + an Optimizer.
  Each **preset** owns its own `validate_components()`.
- **No duplication across methods.** Shared skeletons (gradient assembly, FD
  validation) are single parameterized functions.
- **Component lifecycle:** `build()` (one-time alloc) → `launch()` × N (per-step
  compute) → `destroy()`. No allocation in `launch()`.
- **Two registration paths:** `driver.add(comp)` uses `default_bindings()`;
  `driver.on(phase).add(comp, order)` for explicit binding.
- **Observer vs Recorder:** Observer computes derived quantities; Recorder
  writes I/O. Never both.
- **Hot-path data structures carry only fundamental physical quantities** —
  positions, velocities, integrals, shell descriptors; never precomputed
  indices, caches, or derived maps. Consumers derive their own views.

### Context tiers

`Context` (cpu `component/context.h`, cuda `component/context.cuh`) is a flat
**core** plus per-workflow **service slices**. Core fields serve every
workflow: threading/MPI topology (`n_threads`, `mpi_rank`, `mpi_size`), the
timing registry, and `outer_status` (owned by the category driver). Workflow
state lives in its slice — a component reaches only the slice of the workflow
it serves:

| Slice | Struct | Serves | Key fields (CPU) |
|---|---|---|---|
| `ctx.md` | `DynamicsServices` | Dynamics force pipeline | `dd` (domain decomposition), `pair_list`, halo/migration hooks |
| `ctx.mc` | `SamplingServices` | Sampling | `proposal` (`mc::Transaction` — Mover writes, driver rolls back on reject) |
| `ctx.opt` | `OptimizeServices` | Optimize | `reeval_forces` line-search callback, `line_search_trial` flag |
| `ctx.analysis` | `AnalysisServices` | async observers/recorders | `executor` (snapshot-ring analysis) |

Slices are plain aggregates and the Context stays non-polymorphic
(static-asserted). Per-backend: `cuda::Context` carries streams + the `md`
slice only (`cuda::DynamicsServices`: `pair_list`, `flat_pairs_requested`,
exclusion CSR, …); the cuda SCF solver deliberately does **not** take the
driver Context — it receives a narrow `SolverContext{stream, frame, step}`
(`cuda/scf/internal/solver_context.cuh`, isolation grep-gated by
`tests/cuda/scf_solver_isolation.sh`).

### Frame — multi-backend design

A Frame is the device-resident system state. The **semantic** field set is
shared across backends (atoms incl. `charge`/`tag`/`image`, bonded topology,
the `core/system_scalars.h` `SystemScalars` bundle, box) — but each backend
owns its **representation**, and the differences are load-bearing:

| Frame | Representation | Topology layout |
|---|---|---|
| `molrs.Frame` | flat Block columns — the exchange schema at the boundary | flat bonds/angles/dihedrals blocks |
| `cpu::Frame` | xtensor SoA on host; zero-copy in Driver | per-atom LAMMPS layout (`num_bond`/`bond_atom` with global-tag partners) |
| `cuda::Frame` | raw device pointers + host-cached counts | **dual**: single-GPU flat `BondList`/`AngleList`/`DihedralList` (kernels not yet ported to per-atom, so `cpu::to_cuda` re-flattens on every transfer); multi-GPU per-atom SoA with global-tag partners for migration |
| `torch::Frame` | `at::Tensor` blocks (autograd-capable) | flat bonds/angles/dihedrals sub-structs |

Conversion points, each owning exactly one boundary:

- `materialize()` / `write_back()` (`cpu/frame/materialize.h`) — molrs.Frame ⇄
  cpu::Frame, run once per `run()`, never in a step loop. `write_back`
  follows a pinned round-trip contract (golden-tested in
  `tests/cpu/test_materialize.cpp`): x/v/f, charge, image, and tag are
  written back (tag with checked i32 narrowing — overflow throws, no silent
  truncation); type, mass, topology, and box are per-run static inputs and
  **not** written; scalars/stress are a write-only channel.
- `cpu::to_cuda()` / `cuda::to_cpu()` — host ⇄ device; `cuda::write_back`
  inherits the cpu contract via `to_cpu`.
- `torch::transfer` — cpu ⇄ torch tensors.

Precision tiers map 1:1 across backends: `core/precision.h`
`StorageReal`/`WorkReal`/`AccReal` (accumulators always double) ⇄
`cuda/frame/types.cuh` `real_t`/`work_t`/`prec_t`; integer IDs by LAMMPS
lmptype role (`LocalInt`/`TagInt`/`ImageInt`/`BigInt`).

Per the Frame-purity rule the frame holds **no** method-specific scratch: the
SCF workspace is owned by `*::scf::SCFSolver`, the plane-wave workspace by
`cpu::pw::PWSolver`, the neighbor pair list by `*::NeighborList` (published as
a non-owning `ctx.md.pair_list`). No implicit transfer in low-level operators;
CUDA pipelines stay on device between explicit transfers.

### CUDA execution model

Components classify as **compute** (device-only launch), **CPU bridge** (D2H +
host callback via AsyncBridge), or **host-dependent compute** (documented
exception). See `.claude/NOTES.md` for the classification tables.

### molrs integration

molrs (Rust) is a hard dependency for I/O and trajectory analysis. The engine
*measures*; molrs *analyzes/stores*: bridge calls are raw-array FFI wrappers —
per-frame writers in `cpu/io/frame_writer.h`, batch computes + streaming
accumulators (`cpu::compute::*`) in `cpu/observer/trajectory_metrics.h`
(PIMPL) — with no `rust/cxx.h` in public headers. `molrs.hpp` is the only
include for CXX-generated code. **No analysis math lives in C++ on either
backend**: both the CPU and CUDA RDF observers stream positions + per-frame
cell into the molrs accumulator (binning + normalization are molrs-only,
grep-gated by `test_cuda_rdf_molrs_migration`); the one g(r) file emitter is
`cpu::rdf_write_curve` (`cpu/observer/rdf_curve.h`, CUDA-safe).

## Conventions

- **Namespace:** `ATV`. Sub: `ATV::cpu`, `ATV::cuda`, `ATV::md`,
  `ATV::minimizer`, `ATV::sampling`, `ATV::path`, `ATV::scf`, `ATV::pw`,
  `ATV::evaluate`.
  Inner solvers: `ATV::cpu::scf::SCFSolver`, `ATV::cuda::scf::SCFSolver`,
  `ATV::cpu::pw::PWSolver`, `ATV::cpu::gpw::GPWSolver`,
  `ATV::cpu::semiempirical::NDDOSolver`. SCF/PW internals:
  `ATV::cpu::scf::internal` etc.
- **File extensions:** CPU `.h`/`.cpp`, CUDA `.cuh`/`.cu`, facade/core `.h`.
  Documented exceptions: deliberately CUDA-free bridge headers under `cuda/`
  keep `.h` (`cuda/log_callback.h`, `cuda/timer.h`).
- **Driver template:** use a **preset** (`cpu::md::MD`, `cpu::minimizer::Minimize`,
  `cpu::mc::MC`, `cpu::Evaluate`, `cpu::vib::Hessian`, `cpu::PathDriver`), or subclass the
  loop node `cpu::CPUDriver<PhaseEnum>` / `cuda::CUDADriver<PhaseEnum>` for a genuinely
  new loop shape. Never `Driver<>` directly. There are no category bases.
  Driver subclass headers include the matching `backend_traits.h`.
- **OOP, not free-function soup.** Atomiverse is object-oriented: behaviour lives
  on classes, state is constructed, invariants are established in constructors.
  - **No factory functions.** Components, potentials, basis/PP data, config —
    all are classes with constructors. Never `make_X()` / `builtin_X()` /
    `resolve_X()` / `parse_X()` free functions, and never a new config-builder
    overload where a constructor would do. A lookup is a constructor:
    `GthPotential(Element, name)`, `EAMPotential(Element, name)`.
  - **Do not extract single-use code into a function.** If a block runs exactly
    once, from exactly one call site, it stays inline where it runs. A function
    exists to be called from more than one place, or to be tested on its own —
    not to give a paragraph a name. Premature extraction is what turns a solver
    into a maze of one-shot helpers that must each be read before the algorithm
    is visible. Extract on the *second* caller, not the first.
- **No fallback values** — lookup tables throw on unsupported input.
- **Three test tiers, and they never mix** (layout above):
  - `tests/` — unit only. One function or method per test. The whole suite runs
    in **seconds**; it is the fast-iteration loop. **E2E is strictly forbidden.**
  - `regressions/` — **every delivered feature ships exactly one** minimal
    textbook end-to-end case here, as an out-of-library consumer of the public
    API, checked against a **cited literature value** (not a golden snapshot of
    our own output), `<= 60 s`. Never in `tests/`.
  - `benchmarks/` — external-engine parity and performance, at acceptance time.
- **Phase naming:** Evaluate `prepare, evaluate, finalize`; MD LAMMPS-style
  (`initial_integrate, pre_force, force, post_force, final_integrate,
  end_of_step`); Minimizer `search_direction, line_search, force,
  convergence_check, end_of_step`; Vib `setup, displace, force, accumulate,
  diagonalize, thermochem, end_of_step`; Sampling `propose, energy,
  end_of_step`. **No SCF phase enum** anywhere.
- **Units:** eV, Å, fs, amu, eV/Å, K. Conversions in `core/units.h`. SCF/PW
  report internally in Hartree/Ry and convert eV onto `frame.scalars.energy`.
- **Precision:** LAMMPS-style tiered, in `core/precision.h`. Reals:
  `StorageReal` (per-atom state) / `WorkReal` (kernel math) — lowerable to
  `float` via `-DATV_REAL_STORAGE` / `-DATV_REAL_WORK`; `AccReal` (accumulators,
  energy, virial, SCF) is **always** `double`. Integers by role (LAMMPS
  `lmptype.h`): `LocalInt` (32-bit counts), `TagInt`/`ImageInt` (IDs/image flags,
  32-bit under `-DATV_ID_SIZE=smallbig`, 64-bit under `bigbig`), `BigInt`
  (global totals, always 64-bit). Legacy `Real`=`AccReal`, `Int`=`LocalInt`.
  CUDA mirror: `real_t`/`work_t`/`prec_t` in `cuda/frame/types.cuh`.
- **libxc:** when enabled, native XC must match libxc to 1e-12.
- **ERI performance:** `-march=native` required. xsimd 11.2.0 linked directly to
  `atomiverse` for `native_integrals.cpp`. Do **not** set `XTENSOR_USE_XSIMD`
  globally (it activates SIMD-backed xtensor storage which breaks
  `reshape_view` in `xc_grid_engine.h`). `eri_block()` dispatches by Ltot to a
  compile-time-Ltot template; ss|ss uses an xsimd ket-loop.
- **Doxygen:** every public symbol needs `@brief`, `@param` (with units),
  `@return`.

## Key user-facing APIs (current, verified)

- **Frame:** `cpu::Frame f(n);` then `f.atoms.x[i]`/`.y`/`.z`/`.vx..`/`.mass[i]`/
  `.type[i]` (type = atomic number Z for SCF); `f.box = Mat3::cubic(L)`;
  scalars on `f.scalars.energy` / `f.scalars.kinetic`. Header
  `atomiverse/cpu/frame/frame.h`.
- **Gaussian integrals:** `cpu::scf::build_integrals(frame, BasisSet, charge=0,
  compute_eri=false)`. **`compute_eri` defaults to false** — `rhf_solve` builds
  Fock directly; pass `true` only when you read `integrals.ERI` or feed the CUDA
  RHFFock kernel. `BasisSet`: `STO3G, Pople631G, Pople631Gs, Def2SVP, Def2TZVP`.
- **SCF solves:** `cpu::scf::rhf_solve(ints, tol, max_iter)`,
  `uhf_solve(ints, n_alpha, n_beta, ...)`, `uks_solve(ints, n_alpha, n_beta,
  XCType, n_radial, n_angular, ...)`. `XCType`: `LSDA, PBE, B3LYP, PBE0`.
  Config builders: `rhf_config`, `uhf_config`, `uks_config`, `default_config`
  (Restricted KS / PBE / MINAO).
- **Evaluators (no categorical setter — bind explicitly):**
  `eval.on(evaluate::Phase::evaluate).add(cpu::pw::PWEvaluator{cfg})`;
  `aimd.on(md::Phase::force).add(cpu::scf::SCFEvaluator{cfg})`. Or the
  `cpu::single_point(frame, evaluator)` helper.
- **Plane wave:** `cpu::pw::PWConfig{upf_path, ecut_ry, k_mesh, tol, max_iter,
  mix_beta, pulay_depth}` → `PWSolver{cfg}.build(frame); run(frame)` →
  `PWResult{energy [Ry], energy_ha, converged, iterations}`. Shim:
  `cpu::pw::pw_scf(frame, upf_path, ...)`. Empty `upf_path` = toy `V_loc`
  (runs, unphysical). A real Si UPF ships at
  `benchmarks/pw/si/validation/reference/Si_lda_nc.upf`.
- **Categorical setters** (LAMMPS-shaped, on `cpu::md::MD` /
  `cpu::minimizer::Minimize` / `cpu::mc::MC`): `set_pair, set_bond, set_angle,
  set_dihedral, set_neighbor, set_observer, set_record, set_thermo` (shared);
  `set_integrator, set_thermostat, set_barostat` (MD); `set_optimizer`
  (Minimize); `set_mover` (MC). NeighborList is user-set; halo exchange auto-injects under
  `ATV_BUILD_MPI` when `set_pair` is called; `VirialStress` auto-injects when an
  NPT/Berendsen barostat declares `reads()==Stress`. `add()` /
  `on(phase).add()` remain public for hooks the categorical API lacks.

```cpp
// Single-point RHF (CPU)
auto ints = cpu::scf::build_integrals(frame, cpu::scf::BasisSet::STO3G);
auto res  = cpu::scf::rhf_solve(ints, 1e-10, 100);   // res.energy [Hartree]

// Classical MD (Dynamics) — categorical setters
cpu::md::MD sim;
sim.set_pair(cpu::PairLJCut{0.0104, 3.4, 6.8});
sim.set_neighbor(cpu::md::NeighborList<>{6.8, 0.3});
sim.set_integrator(cpu::md::VelocityVerlet{0.01});
sim.set_thermo(100);
sim.build();  sim.run(frame, 1000);

// AIMD (Dynamics + SCFEvaluator force provider) — explicit bind (no default_bindings)
cpu::md::MD aimd;
aimd.set_integrator(cpu::md::VelocityVerlet{0.5});
aimd.on(md::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
aimd.build();  aimd.run(frame, 100);

// Geometry optimization
cpu::minimizer::Minimize opt;
opt.set_optimizer(cpu::minimizer::LBFGS{5, 0.01, 10});   // memory, alpha0, max_ls
opt.add_handler(cpu::ForceConverged{0.05, 0.03, 0.1, 0.05, 1e-3});  // the tolerances
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{cfg});
opt.build();  opt.run(frame, 100);

// Monte Carlo sampling (NVT Metropolis)
cpu::mc::MC mc{/*T [K]=*/300.0, /*seed=*/12345u};
mc.set_pair(cpu::PairLJCut{eps, sigma, cutoff});
mc.set_neighbor(cpu::md::NeighborList<>{cutoff, skin, 1, 0, true});
mc.set_mover(cpu::mc::TranslateMove{/*delta [Å]=*/0.2});
mc.set_thermo(100);
mc.build();  mc.run(frame, 100000);   // 100k trial moves
double e = mc.current_energy();       // [eV]; mc.acceptance_ratio() in [0,1]

// Plane-wave DFT single point
cpu::pw::PWEvaluator pwe{cpu::pw::PWConfig{.ecut_ry = 30.0, .k_mesh = 4}};
cpu::single_point(frame, pwe);   // frame.scalars.energy [eV]; pwe.last_result() [Ry/Ha]

// Termination is a policy — compose handlers, or write your own (LAMMPS `fix halt`)
class StopWhenTooHot : public cpu::LoopHandler {
  public:
    explicit StopWhenTooHot(double t_max) : t_max_(t_max) {}
    const char* name() const override { return "stop_when_too_hot"; }
    ATV::Verdict after_iteration(const cpu::LoopState& s) override {
        return (s.frame.scalars.temperature > t_max_) ? ATV::Verdict::Stop
                                                      : ATV::Verdict::Continue;
    }
  private:
    double t_max_;
};
sim.add_handler(cpu::StopOnNaN{});            // plural — Stop > Converged > Continue
sim.add_handler(StopWhenTooHot{2000.0});      // written OUTSIDE the library
sim.run(frame, 1000);
if (sim.converged()) { /* ... */ }            // the LOOP knows; the algorithm does not
sim.stop_reason();                            // which handler ended it

// A Driver is a node: nest one in another, on ONE shared Runtime
cpu::minimizer::Minimize inner;                     // its own phases + handlers
inner.set_optimizer(cpu::minimizer::LBFGS{5, 0.01, 10});
inner.on(minimizer::Phase::force).add(forces);
cpu::md::MD outer;                                  // a DIFFERENT phase enum
outer.set_integrator(cpu::md::VelocityVerlet{0.5});
outer.on(md::Phase::end_of_step).add(cpu::SubDriver{inner, /*budget=*/50});
outer.build();                                      // builds BOTH; inner allocates once
outer.run(frame, 1000);                             // inner's convergence CANNOT stop outer
```

## Implemented components

### Drivers — all **presets** over one loop node; there is no category

| Preset | Backend | Phase enum | Default policy |
|--------|---------|------------|----------------|
| `cpu::Evaluate` / `cuda::Evaluate` | CPU/CUDA | `evaluate::Phase` | — (one pass) |
| `cpu::md::MD` / `cuda::md::MD` | CPU/CUDA | `md::Phase` | — (the step budget) |
| `cpu::minimizer::Minimize` / `cuda::minimizer::Minimize` | CPU/CUDA | `minimizer::Phase` | handler `ForceConverged` |
| `cpu::vib::Hessian` | CPU | `vib::Phase` | — (fixed 6N + finalize) |
| `cpu::mc::MC` | CPU | `sampling::Phase` | **outcome** `Metropolis` |
| `cpu::PathDriver` | CPU | `path::Phase` | handler `AllReplicasConverged` |

Each derives the loop node (`cpu::CPUDriver<P>` / `cuda::CUDADriver<P>`) **directly** and
configures its defaults in its **constructor**. `cpu::SubDriver<P>` wraps any of them as a
node in another driver's Plan (nesting, on one shared Runtime).

SCF/PW/GPW/NDDO are inner solvers (`*::scf::SCFSolver`, `cpu::pw::PWSolver`,
`cpu::gpw::GPWSolver`, `cpu::semiempirical::NDDOSolver`), not drivers.

### MD components
VelocityVerlet, NVELimit, InitVelocity, NoseHoover (Thermostat), NoseHooverNPT
(Integrator+thermostat+barostat), FixLangevin (Thermostat + Force), Barostat
(Berendsen), Deform, NeighborList, PairLJCut (Energy|Force|Stress),
BondHarmonic / AngleHarmonic / DihedralMultiHarmonic (Energy|Force),
VirialStress. **LAMMPS-style names are canonical** — the old
`LJ`/`HarmonicBond`/`HarmonicAngle`/`MultiHarmonicDihedral` CPU aliases are
deleted (the torch backend keeps its own `LJ`/`HarmonicBond`).

### MC components (CPU)
Three presets over one shared loop (`mc::McBase`): `mc::MC` (NVT), `mc::NptMc`
(NPT), `mc::Gcmc` (μVT) — each attaches its `OutcomePolicy` (`Metropolis`,
`NptMetropolis`, `GcmcMetropolis`) in its constructor. Movers: `TranslateMove`
(single-particle displacement δ [Å]), `VolumeMove` (isotropic, uniform-in-V),
`InsertDeleteMove` (single species; tags from `Atom::acquire_tag`), composed
one-trial-per-sweep by `MoveSelector{}.add(mover, weight)…` (ref-qualified
chaining — a temporary chain stays an rvalue so the driver owns the copy). The
trial changeset is `mc::Transaction`, a **Frame patch** (touched atoms'
original rows + box + pre-trial `nlocal`; rows travel via
`comm::pack/unpack_exchange`), published on `ctx.mc.proposal`: the Mover
records before mutating, the policy decides, `Transaction::rollback` splices
the patch back bitwise (topology included). Per-verb attempt/accept tallies on
`McBase::attempts/accepts(verb)`. Reuses the MD force pipeline for trial
energies (`set_pair`/`set_neighbor`); `NeighborList` runs in always-rebuild
mode; trial energy is a deliberate O(N) full recompute. Single-rank v1.

### Path components (CPU)
`cpu::path::NebCoupling` (Evaluator|Force, bound at `couple`) — Henkelman improved
tangent + spring + climbing image; overwrites each replica's true force with the NEB
force. `cpu::path::BandLBFGS` (Optimizer) — L-BFGS over the concatenated band
(`n_images x 3N`), reusing `minimizer::LbfgsHistory`; trust-radius step, no line search
(NEB forces are not the gradient of anything, so there is no phi(alpha) to search).
Expansions: `cpu::LinearInterpolation` (NEB; owns the product state as its own config)
and `cpu::IdentityReplication` (the PIMD/REMD shape; no in-tree user yet).

### Minimizer components
LBFGS, CG (Polak-Ribière+), OptLog / TrajectoryWriter / FinalGeometry
(Recorders). Optimizer is hot-swappable on the same `Minimize` driver. MPI:
`install_defaults()` folds Migrate/Halo/Reduce into the `force` phase; LBFGS/CG
are DD-aware; recorders gather owners to rank 0.

### Gaussian SCF (`-DATV_BUILD_RHF=ON`)
Inner solvers `cpu::scf::SCFSolver` / `cuda::scf::SCFSolver` (Hamiltonian →
diagonalize → density → DIIS → convergence; loop internal). Evaluators
`cpu::scf::SCFEvaluator` (Energy|DensityMatrix|Orbitals, +Force with
`with_gradient(true)`; plugs into Evaluate/MD/Minimize) and
`cuda::scf::SCFEvaluator` (energy/density only — no CUDA SCF gradient).
Gradients: `scf_gradient(frame, integrals, result, config)` (HF + KS).
Utilities: `basis_eval.h`, `grid.h`, `xc_functional.h`, `aspc_predictor.h`.

### GPW data (in-tree static `.h` tables, no parsing)
GTH pseudopotentials, MOLOPT/PADE basis sets and the DFT-D3 reference tables
are **static C++ data headers** (`cpu/gpw/basis_data/gpw_tables.h`,
`basis_data/dftd3_tables.h`), transcribed once from the published tables
(Krack 2005; VandeVondele-Hutter 2007; Grimme 2010; generator preserved at
`tools/gen_gpw_tables.cpp`). Access is **OOP lookup constructors only**:
`GthPotential(Element, name)` / `GaussianBasis(Element, name)`
(`cpu/gpw/gth_potential.h`, `gaussian_basis.h`) and `D3Params{}`
(`internal/dftd3_params.h`). `Element` is molrs's periodic-table enum
(`atomiverse/molrs.hpp` aliases it as `ATV::Element`) — the single
element-identity type everywhere; a GPW `Kind{element, basis, pot}`
(CP2K's `&KIND`) pairs it with its data. **There are no parsers, no
`builtin_*`/`make_*` factory functions, and nothing reads any data file at
build, test or run time.** Unsupported element/name throws (no fallback).

`GpwConfig.local_pp_on_grid` (default true) selects CP2K's route of collocating
the GTH local short-range PP on the grid. **Energy and analytic force both read
it** — the gradient must be taken of the functional the SCF converged, or the two
disagree by the grid discretization error of V_loc^SR (this was a real bug: ZnO
analytic-vs-FD 0.53 Ha/Bohr). See `cpu/gpw/internal/local_pp.h`.

### Plane-wave DFT (`-DATV_BUILD_PW=ON`, FFTW3)
Inner solver `cpu::pw::PWSolver` (Poisson → XC → Veff → per-k Davidson →
density + Pulay mix → convergence; `build()` parses UPF, geometry setup rebuilt
per `run()`). Evaluator `cpu::pw::PWEvaluator` (Energy; writes
`frame.scalars.energy` [eV], Ry/Ha on `last_result()`). Shim `pw_scf()`. MPI:
k-point round-robin over `MPI_COMM_WORLD`. CUDA: `cuda::pw::pw_scf` runs the
Davidson hot path on device (cuFFT/cuBLAS) and reuses the CPU `internal`
helpers for parity (`tests/cuda/pw/test_pw_si.cu`).

### Vibrational (CPU)
`cpu::vib::Hessian` (Response): finite-difference Hessian + thermochemistry.

### Output (cross-driver)
Thermo, Recorder (XYZ), RDF, VirialStress, Mulliken, density.

## Skills & agents

General workflow skills (impl / fix / debug / refactor / review / docs / test /
spec / litrev / note) come from the **molcrafts `mol` plugin** (`/mol:impl`,
`/mol:fix`, …), auto-adapted via the `mol_project:` frontmatter above.

Atomiverse-specific skills layer on top for domain audits: `/av-arch`,
`/av-cuda`, `/av-hpc`, `/av-perf` — each delegating to a paired agent
(`av-architect`, `av-cuda`, `av-hpc`, `av-optimizer`) that encodes
Atomiverse-only knowledge (CUDA execution-model classification, xsimd 11.2.0 /
xtensor 0.25.0 constraints, ERI shell-block dispatch + `-march=native`, Frame
purity, dual-backend hot-path anti-patterns).

Knowledge flow: `.claude/NOTES.md` (evolving) → CLAUDE.md (stable) ← agents
(read-only). Specs live under `.claude/specs/` while in flight and are deleted
at `done` (completed-spec history is the git log).
