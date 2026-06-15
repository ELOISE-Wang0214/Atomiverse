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
  notes_path: .claude/NOTES.md
  specs_path: .claude/specs/
---

# CLAUDE.md

Atomiverse is a dual-backend (CPU + CUDA) molecular-simulation engine with one
shared Driver/Component API spanning classical molecular dynamics, geometry
optimization, vibrational analysis, and electronic structure (Gaussian-basis
SCF and plane-wave DFT). C++17, `ATV::` namespace.

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
| `ATV_BUILD_KSPACE` | OFF | Reciprocal-space electrostatics PME (needs FFTW3) |
| `ATV_BUILD_TORCH` | OFF | LibTorch backend; defines `ATV_HAS_TORCH` |
| `ATV_BUILD_TESTS` / `ATV_BUILD_BENCHMARKS` / `ATV_BUILD_EXAMPLES` | OFF | test / benchmark / example targets |
| `ATV_USE_MATPLOT` | OFF | plot generation in selected examples (needs gnuplot) |
| `ATV_REAL` / `ATV_INT` | `double` / `int` | compile-time `Real` / `Int` precision |

Performance: `-march=native` is required for the CPU SCF ERI kernels.

## Repository layout

```
include/atomiverse/ + src/   mirrored header/impl trees (architecture below)
tests/                       unit tests — fast correctness only; never timing or
                             external-engine comparison (those are benchmarks)
benchmarks/                  <domain>/<case>/{validation,performance} — see benchmarks/README.md
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
  io/          external formats / molrs bridge (molrec, lammps_data; CPU only)
  device/      CUDA/torch runtime utilities (device_manager, error, launch, timer)
  comm/        infra: domain decomposition, halo exchange, migration, reductions
  neighbor/    infra: cell/Verlet neighbor structures + the published PairList
  async/       infra: async analysis pipeline (snapshot ring, D2H bridge, executors)
  md/          domain: integrators, thermostats, barostats, deform, velocity init
  minimizer/   domain: SD, CG, L-BFGS + opt recorders
  scf/         domain: Gaussian-basis SCF (solver, integrals, gradients, evaluator)
  pw/          domain: plane-wave DFT (PWSolver, PWEvaluator)
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

**Box** (`cpu::Box` in `cpu/box.h`, `cuda::Box` in `cuda/box.cuh`) is the
canonical cell type at each backend root: owns the cell matrix plus a
per-step-immutable cache (`inv_cell`, `inv_L`, `half_L`, `ortho`) and exposes
`length(d)`, `volume()`, `min_image(...)`, `wrap_axis(...)`, `scale_*`,
`set_cell`. CPU adds an xsimd `min_image_batch<B>` in `cpu/box_simd.h`. LJ keeps
an inlined branchy min-image in its xsimd ket-loop because member-load-through-
`this` defeats hoisting there (~7%/step).

### Facade (CUDA-free)

```
core/        Real, Int, Mat3, units, constants, device tags, backend_traits
facade/      Driver<PhaseEnum, DeviceTag>, Plan, Stage (phase enums)
potentials/  backend-agnostic pair force/energy kernels shared by CPU/CUDA
molrs.hpp    stable CXX FFI header
```

### Driver taxonomy

CPU has four **outer task-category** drivers, each subclassing
`CPUDriver<PhaseEnum>` (which wraps `Driver<PhaseEnum, DeviceTag>`). Category is
decided by **termination semantics**, not algorithm:

| Category base | Phase enum | Termination | Concrete CPU driver |
|---|---|---|---|
| `cpu::EvaluateDriver` | `evaluate::Phase` | exactly one pass | `cpu::Evaluate` |
| `cpu::DynamicsDriver` | `md::Phase` | full step budget; honours `outer_status.stop_requested` | `cpu::md::MD` |
| `cpu::OptimizeDriver` | `minimizer::Phase` | early stop on `outer_status.converged` | `cpu::minimizer::Minimize` |
| `cpu::ResponseDriver` | `vib::Phase` | fixed 6N-displacement budget + finalize; no early stop | `cpu::vib::Hessian` |

CUDA mirrors three of these: `cuda::Evaluate`, `cuda::md::MD`,
`cuda::minimizer::Minimize`, each subclassing `cuda::CUDADriver<PhaseEnum>`.

Adding a fifth user-facing category requires a genuinely distinct outer-loop
semantic (e.g. Sampling for MC/HMC/REMD, Path for NEB) — never an algorithm
wrapper.

**SCF is *not* an outer category.** Both Gaussian SCF and plane-wave DFT are
*inner solvers* — public low-level classes, **not** `Driver` subclasses —
driven by an Evaluator component:

- Gaussian: `cpu::scf::SCFSolver` / `cuda::scf::SCFSolver`, behind
  `cpu::scf::SCFEvaluator` / `cuda::scf::SCFEvaluator`.
- Plane-wave: `cpu::pw::PWSolver`, behind `cpu::pw::PWEvaluator`.

Single-point energy in user code is `Evaluate + <method>Evaluator` (or the
`single_point` helper). There is **no `SinglePoint` class** and **no SCF phase
enum** — `ATV::scf::Phase` is deleted and grep-gated (`test_scf_phase_scrub`,
`test_cuda_scf_phase_scrub`). The `scf_solve()` / `pw_scf()` shims remain as
thin functional wrappers over a transient solver. SCF/PW step helpers
(`*::scf::internal::*`, `*::pw::internal::*`) are owned by the solver, never
inherit `ComponentBase`, and are never registered into a Plan; facade-public
headers must not include them (tests may).

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
  `Observer`, `None`.
- `Capability capabilities()` — bitset of frame fields `launch()` may write:
  `Energy`, `Force`, `Stress`, `Charge`, `Dipole`, `DensityMatrix`,
  `Orbitals`. Bitwise `|`/`&`/`has_capability`.

**Algorithms are components, not driver subclasses.** `LBFGS`/`CG`/future
`FIRE` are Optimizer components plugged into `Minimize`;
`VelocityVerlet`/`NoseHoover`/`FixLangevin`/`Barostat` are
Integrator/Thermostat/Barostat components plugged into `MD`.

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

### Data residency

The three frames carry an **identical field set** (atoms incl.
`charge`/`tag`/`image`, bonds/angles/dihedrals, the `core/system_scalars.h`
`SystemScalars` bundle, box) and differ only in device representation:

- `cpu::Frame` = xtensor on host. Zero-copy in Driver (no transfer).
- `cuda::Frame` = device pointers. `cpu::to_cuda()` / `cuda::to_cpu()`.
- `torch::Frame` = `at::Tensor` (autograd-capable). `torch::transfer` bridges.
- `molrs.Frame` ⇄ a device frame via `materialize()` (H2D) / `write_back()`
  (D2H), run once per `run()`.

Per the Frame-purity rule the frame holds **no** method-specific scratch: the
SCF workspace is owned by `*::scf::SCFSolver`, the plane-wave workspace by
`cpu::pw::PWSolver`, the neighbor pair list by `*::NeighborList` (published as a
non-owning `Context::pair_list`). No implicit transfer in low-level operators;
CUDA pipelines stay on device between explicit transfers.

### CUDA execution model

Components classify as **compute** (device-only launch), **CPU bridge** (D2H +
host callback via AsyncBridge), or **host-dependent compute** (documented
exception). See `.claude/NOTES.md` for the classification tables.

### molrs integration

molrs (Rust) is a hard dependency for I/O. `cpu::MolRec` uses PIMPL — no
`rust/cxx.h` in public headers. `molrs.hpp` is the only include for
CXX-generated code.

## Conventions

- **Namespace:** `ATV`. Sub: `ATV::cpu`, `ATV::cuda`, `ATV::md`,
  `ATV::minimizer`, `ATV::scf`, `ATV::pw`, `ATV::evaluate`. Inner solvers:
  `ATV::cpu::scf::SCFSolver`, `ATV::cuda::scf::SCFSolver`,
  `ATV::cpu::pw::PWSolver`. SCF/PW internals: `ATV::cpu::scf::internal` etc.
- **File extensions:** CPU `.h`/`.cpp`, CUDA `.cuh`/`.cu`, facade/core `.h`.
  Documented exceptions: deliberately CUDA-free bridge headers under `cuda/`
  keep `.h` (`cuda/log_callback.h`, `cuda/timer.h`).
- **Driver template:** subclass a category base
  (`cpu::EvaluateDriver`/`DynamicsDriver`/`OptimizeDriver`/`ResponseDriver` or
  the CUDA equivalents), never `Driver<>` directly. Driver subclass headers
  include the matching `backend_traits.h`.
- **No factory functions** — components are classes with constructors.
- **No fallback values** — lookup tables throw on unsupported input.
- **Phase naming:** Evaluate `prepare, evaluate, finalize`; MD LAMMPS-style
  (`initial_integrate, pre_force, force, post_force, final_integrate,
  end_of_step`); Minimizer `search_direction, line_search, force,
  convergence_check, end_of_step`; Vib `setup, displace, force, accumulate,
  diagonalize, thermochem, end_of_step`. **No SCF phase enum** anywhere.
- **Units:** eV, Å, fs, amu, eV/Å, K. Conversions in `core/units.h`. SCF/PW
  report internally in Hartree/Ry and convert eV onto `frame.scalars.energy`.
- **Precision:** compile-time `Real` (`-DATV_REAL`), `Int` (`-DATV_INT`).
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
  `cpu::minimizer::Minimize`): `set_pair, set_bond, set_angle, set_dihedral,
  set_neighbor, set_observer, set_record, set_thermo` (shared);
  `set_integrator, set_thermostat, set_barostat` (MD); `set_optimizer`
  (Minimize). NeighborList is user-set; halo exchange auto-injects under
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

SCF/PW are inner solvers (`*::scf::SCFSolver`, `cpu::pw::PWSolver`), not drivers.

### MD components
VelocityVerlet, NVELimit, InitVelocity, NoseHoover (Thermostat), NoseHooverNPT
(Integrator+thermostat+barostat), FixLangevin (Thermostat + Force), Barostat
(Berendsen), Deform, NeighborList, PairLJCut (Energy|Force|Stress),
BondHarmonic / AngleHarmonic / DihedralMultiHarmonic (Energy|Force),
VirialStress. **LAMMPS-style names are canonical** — the old
`LJ`/`HarmonicBond`/`HarmonicAngle`/`MultiHarmonicDihedral` CPU aliases are
deleted (the torch backend keeps its own `LJ`/`HarmonicBond`).

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
