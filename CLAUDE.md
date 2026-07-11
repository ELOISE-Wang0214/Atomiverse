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
                             method at ms-scale. No end-to-end runs (a full
                             SCF/MD/vib solve or MC equilibration to
                             convergence), no timing, no external-engine
                             comparison, and no third-party package read at
                             build or test time. There is no opt-in switch to
                             run e2e from ctest: if it is registered, it runs.
benchmarks/                  <domain>/<case>/{validation,performance} — see benchmarks/README.md
                             Home of every e2e / parity check (val_*), run at
                             acceptance time, never by ctest.
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
  driver/      Driver pillar: category bases, concrete Evaluate, single_point,
               step executor, backend_traits, step-loop instrumentation
  component/   Component pillar: ComponentBase, Launch/Context, generic + shared bases
  io/          external formats / molrs bridge (extxyz, frame_writer, lammps_data; CPU only)
  device/      CUDA/torch runtime utilities (device_manager, error, launch, timer)
  comm/        infra: domain decomposition, halo exchange, migration, reductions
  neighbor/    infra: cell/Verlet neighbor structures + the published PairList
  async/       infra: async analysis pipeline (snapshot ring, D2H bridge, executors)
  linalg/      infra: shared dense linear algebra (eigensolve, density build; CPU only)
  md/          domain: integrators, thermostats, barostats, deform, velocity init
  mc/          domain: Monte Carlo — SamplingDriver's concrete MC, movers,
               trial transaction (CPU only)
  minimizer/   domain: SD, CG, L-BFGS + opt recorders
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

### Driver taxonomy

A Driver category is a **task type** — the analogue of a Gaussian job type
(SP, Opt, Freq, …): what kind of calculation the user asked for. Category is
decided by **termination semantics**, never by algorithm or method. CPU has
five category bases, each subclassing `CPUDriver<PhaseEnum>` (which wraps
`Driver<PhaseEnum, DeviceTag>`):

| Category base | Gaussian analogue | Phase enum | Termination | Concrete CPU driver |
|---|---|---|---|---|
| `cpu::EvaluateDriver` | SP | `evaluate::Phase` | exactly one pass | `cpu::Evaluate` |
| `cpu::DynamicsDriver` | MD | `md::Phase` | full step budget; honours `outer_status.stop_requested` | `cpu::md::MD` |
| `cpu::OptimizeDriver` | Opt | `minimizer::Phase` | early stop on `outer_status.converged` | `cpu::minimizer::Minimize` |
| `cpu::ResponseDriver` | Freq | `vib::Phase` | fixed 6N-displacement budget + finalize; no early stop | `cpu::vib::Hessian` |
| `cpu::SamplingDriver` | MC | `sampling::Phase` | trial-move budget with inline Metropolis accept/reject; honours `stop_requested`, no convergence path | `cpu::mc::MC` |

CUDA mirrors three of these: `cuda::Evaluate`, `cuda::md::MD`,
`cuda::minimizer::Minimize`, each subclassing `cuda::CUDADriver<PhaseEnum>`.

Adding a sixth user-facing category requires a genuinely distinct outer-loop
semantic (e.g. Path for NEB — terminate when every image of a coupled chain
has converged) — never an algorithm wrapper.

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
CPUDriver<PhaseEnum>::run(frame, ...)
  setup_pass(frame)
  for each iteration (semantic depends on category):
    reset outer_status
    backend.step(idx)
    consult outer_status (Optimize: break on converged; Dynamics: on stop_requested)
  teardown_pass()
```

User code does not subclass `Driver<>` directly — use a category base. Concrete
drivers add only infrastructure components in their constructor; no runtime
behaviour beyond the base.

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
`Context::outer_status` is `{converged, stop_requested, reason}`, owned by the
category driver and reset before each iteration. Optimizer components write it
during their `convergence_check`-phase binding. SCF/PW internals never touch it.

**Rules:**
- **Driver owns the loop.** Evaluators, integrators, thermostats are
  Components — never standalone classes with their own loops.
- **Capability over role for evaluator validation.** `EvaluateDriver` accepts
  any `role()==Evaluator` with non-empty capabilities; `DynamicsDriver`
  requires any component with `Capability::Force`; `OptimizeDriver` requires an
  Energy/Force evaluator + an Optimizer component.
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
  `ATV::minimizer`, `ATV::sampling`, `ATV::scf`, `ATV::pw`, `ATV::evaluate`.
  Inner solvers: `ATV::cpu::scf::SCFSolver`, `ATV::cuda::scf::SCFSolver`,
  `ATV::cpu::pw::PWSolver`, `ATV::cpu::gpw::GPWSolver`,
  `ATV::cpu::semiempirical::NDDOSolver`. SCF/PW internals:
  `ATV::cpu::scf::internal` etc.
- **File extensions:** CPU `.h`/`.cpp`, CUDA `.cuh`/`.cu`, facade/core `.h`.
  Documented exceptions: deliberately CUDA-free bridge headers under `cuda/`
  keep `.h` (`cuda/log_callback.h`, `cuda/timer.h`).
- **Driver template:** subclass a category base
  (`cpu::EvaluateDriver`/`DynamicsDriver`/`OptimizeDriver`/`ResponseDriver`/
  `SamplingDriver` or the CUDA equivalents), never `Driver<>` directly.
  Driver subclass headers include the matching `backend_traits.h`.
- **No factory functions** — components are classes with constructors.
- **No fallback values** — lookup tables throw on unsupported input.
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
opt.set_optimizer(cpu::minimizer::LBFGS{5, 0.01, 0.05});
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
```

## Implemented components

### Drivers

| Driver | Backend | Category | Phase enum |
|--------|---------|----------|------------|
| `cpu::Evaluate` / `cuda::Evaluate` | CPU/CUDA | Evaluate | `evaluate::Phase` |
| `cpu::md::MD` / `cuda::md::MD` | CPU/CUDA | Dynamics | `md::Phase` |
| `cpu::minimizer::Minimize` / `cuda::minimizer::Minimize` | CPU/CUDA | Optimize | `minimizer::Phase` |
| `cpu::vib::Hessian` | CPU | Response | `vib::Phase` |
| `cpu::mc::MC` | CPU | Sampling | `sampling::Phase` |

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
`TranslateMove` (Mover — single-particle random displacement, max step δ [Å]);
the trial changeset is `mc::Transaction`, published on `ctx.mc.proposal`
(Mover writes the proposal, the driver commits or rolls back per the inline
Metropolis test). Reuses the MD force pipeline for trial energies
(`set_pair`/`set_neighbor`); `NeighborList` runs in always-rebuild mode.

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

### GPW data (in-tree, no external package)
GTH-PBE pseudopotentials, DZVP-MOLOPT-SR-GTH basis sets and the DFT-D3 reference
C6 table are **compiled in** (`cpu/gpw/basis_data/gpw_data.h` →
`builtin_gth_potential` / `builtin_gth_basis`; `internal/dftd3_params.h` →
`builtin_dftd3_params`), transcribed once from the published tables (Krack 2005;
VandeVondele-Hutter 2007; Grimme 2010). **Nothing in the tree reads a third-party
data file at build or test time.** Unsupported element/name throws (no fallback).
The `cp2k_io.h` / `cp2k_basis.h` file parsers remain as an optional ingest path
for data outside the supported set, but no test or example depends on them.

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
