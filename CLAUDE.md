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

## Build & Test

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
./build/tests/test_md          # single test
./build/benchmarks/bench_md_step  # requires -DATV_BUILD_BENCHMARKS=ON
```

Requires: C++17, CMake ≥ 3.18, CUDA toolkit, Rust toolchain (for molrs), xtensor (auto-fetched).

### Toolchain via modules — never hardcode library paths

Use **Lmod modules**, not hardcoded `/software/...` paths. Run everything in a **login shell** (`bash -l`) so Lmod is defined and the build env auto-loads:

- **CPU (x86_64/epyc login node):** `buildenv-gcc/2026.03` auto-loads in a login shell (GCC 14.3 + OpenBLAS/LAPACK/ScaLAPACK).
- **GPU (aarch64/GH200 node):** `module load GPU/buildenv-gcccuda/2026.03-cu13.0` (GCC 14.3, CUDA 13.0, GPU-aware MPICH, NCCL, cuFFT, FFTW, OpenBLAS).

**Runtime gotcha (the #1 source of `libopenblas.so.0` / `GLIBCXX_3.4.32 not found` at test time):** the buildenv module sets `LIBRARY_PATH` + `*_ROOT` (link time) but **not** `LD_LIBRARY_PATH` (runtime). Before `ctest` / running any binary, mirror it:

```bash
bash -lc 'export LD_LIBRARY_PATH="$LIBRARY_PATH:$LD_LIBRARY_PATH"; ctest --test-dir build --output-on-failure'
```

This resolves OpenBLAS + the GCC-14 libstdc++ with **zero hardcoded paths**, on both the login node and the GPU node. CMake BLAS discovery uses `$OPENBLAS_ROOT` (module-provided); pass `-DBLA_VENDOR=OpenBLAS`.

### GPU node (CUDA build & test on Arrhenius)

The **login nodes are x86_64 with no CUDA toolkit** — the CUDA backend (`ATV_BUILD_CUDA=ON`, any `.cu`) can only be built and tested on the **aarch64 Grace/GH200 GPU node**. The CPU backend builds fine on the login node with `-DATV_BUILD_CUDA=OFF`.

- **Interactive GPU node:** `hpc gpu` (NAISS `interactive` srun wrapper; defaults from `~/.hpc/config.sh`: partition `gpu`, account `naiss2026-4-715-gpu`, 1 GPU, 3 h).
- **Batch / non-interactive** (use this from an automated agent — `hpc gpu` needs a TTY):
  ```bash
  sbatch -p gpu -A naiss2026-4-715-gpu --gpus 1 -t 01:00:00 -c 8 -o out.log job.sh
  ```
- **Inside the job script** (`#!/bin/bash -l` so Lmod's `module` is defined):
  ```bash
  #!/bin/bash -l
  module load GPU/buildenv-gcccuda/2026.03-cu13.0   # NB: never pipe `module load` (subshell loses env)
  # GCC 14.3, CUDA 13.0 (nvcc on PATH, GH200/sbsa), GPU-aware MPICH, NCCL,
  # cuFFT/cuFFTMp, FFTW, OpenBLAS — all via CPATH / LIBRARY_PATH.
  cmake -B build-cuda -DATV_BUILD_CUDA=ON -DATV_BUILD_TESTS=ON \
        -DATV_MOLRS_PATH=$HOME/work/molcrafts/molrs   # ../molrs breaks from a worktree
  cmake --build build-cuda -j$(nproc)
  ctest --test-dir build-cuda --output-on-failure
  ```
- On the GPU node the **default rustup toolchain (aarch64) is correct** — do NOT set `RUSTUP_TOOLCHAIN` (that override is only needed on the x86_64 login node, whose default is wrongly aarch64).
- cuFFT-backed PME and other `ATV_BUILD_KSPACE` CUDA code link `CUDA::cufft` automatically; FFTW/cuFFT are on `CPATH`/`LIBRARY_PATH` after the module load.

## Repository layout

```
include/atomiverse/ + src/   mirrored header/impl trees (layout below)
tests/                       unit tests — fast correctness only; never timing or
                             external-engine comparison (those are benchmarks)
benchmarks/                  <domain>/<case>/{validation,performance} — see benchmarks/README.md
examples/                    flat numbered tutorial sequence — see examples/README.md
tools/                       data generators (gen_basis_data.py, gen_minao_data.py)
                             and debug harnesses (debug_amorphous/)
python/                      scikit-build-core + nanobind wheel (drives root CMake)
docs/, cmake/                narrative docs; CMake helper modules
```

## Architecture

Dual-backend (CPU + CUDA) molecular simulation platform. Three pillars: **Driver** (owns loop), **Component** (single-step operation), **Frame** (device-resident state).

### Backend directory layering

Each backend (`cpu/`, `cuda/`, `torch/`) contains **directories only — no
loose files at the backend root**. Every directory is one of: a pillar layer
(Frame / Driver / Component), an infrastructure cluster, or a method domain
(one class of components):

```
<backend>/
  frame/       ← Frame pillar: device-resident state + residency ops
                 (frame, atom, box(+box_simd), transfer, materialize)
  driver/      ← Driver pillar: category bases (evaluate/dynamics/optimize/
                 response_driver), concrete Evaluate + single_point, backend
                 step executor, backend_traits, step-loop instrumentation
  component/   ← Component pillar: ComponentBase, LaunchContext/Context,
                 generic components (Zero*), shared bases (PairEvaluator)
                 and helpers (owned_range)
  io/          ← external formats / molrs bridge (molrec, lammps_data; CPU only)
  device/      ← CUDA/torch runtime utilities (device_manager, error, launch,
                 timer, log_callback / stream_guard)
  comm/        ← infra: domain decomposition, halo exchange, migration,
                 reductions (MPI on CPU; multi-GPU + NCCL on CUDA)
  neighbor/    ← infra: cell/Verlet neighbor structures + the published PairList
  async/       ← infra: async analysis pipeline (snapshot ring, D2H bridge,
                 executors)
  md/          ← domain: time integration (integrators, thermostats, barostats,
                 deform, velocity init)
  minimizer/   ← domain: SD, CG, L-BFGS + opt recorders
  scf/         ← domain: Gaussian-basis SCF (solver, integrals, gradients,
                 evaluator)
  pw/          ← domain: plane-wave DFT (PWSolver, PWEvaluator)
  vib/         ← domain: finite-difference Hessian, thermochemistry (CPU only)
  kspace/      ← domain: reciprocal-space electrostatics (Ewald, PME)
  potentials/  ← domain: real-space force-field evaluators (pair + bonded styles)
  observer/    ← domain: derived quantities (RDF, virial stress, Mulliken, density)
  recorder/    ← domain: I/O (XYZ, cube, thermo, timing report)
```

**Rules:**
- **No loose files at a backend root.** New shared code goes into the
  matching pillar directory (`frame/`, `driver/`, `component/`, `io/`,
  `device/`); a new cohesive multi-file infra cluster gets its own root-level
  directory (like `comm/`, `neighbor/`, `async/`). Never park shared code
  inside a domain subdirectory.
- **Orthogonality: the same feature lives in the same directory on every
  backend.** Virial stress is `observer/`, neighbor structures are
  `neighbor/`, decomposition is `comm/` — on CPU and CUDA alike.
- Electrostatics splits LAMMPS-style: real-space pair styles in
  `potentials/` (`coul_long`), reciprocal-space solvers in `kspace/`.
- Namespaces are unchanged by the 2026-06 layout passes (`cpu::md::mpi::Reduce`
  lives in `cpu/comm/`, `cpu::md::NeighborList` in `cpu/neighbor/`,
  `cpu::Frame` in `cpu/frame/`); aligning namespaces to directories is a
  deliberate follow-up, not a drive-by.

**Box** (`cpu::Box` in `cpu/box.h`, `cuda::Box` in `cuda/box.cuh`) is the canonical simulation-cell type at each backend root. It owns the cell matrix plus a per-step-immutable cache (`inv_cell`, per-axis `inv_L`, `half_L`, `bool ortho`) and exposes `length(d)`, `volume()`, `min_image(dx, dy, dz)`, `wrap_axis(q, d)`, `scale_diagonal/isotropic/axis`, `set_cell`. The CPU side adds an xsimd `min_image_batch<B>` overload in `cpu/box_simd.h` for SIMD callers. All scattered `Mat3 box(d, d)` reads and the duplicated `static inline minimum_image` / `wrap_to_box` helpers from potentials/, md/, observer/, scf/, and comm/ have been folded into Box; LJ keeps an inlined branchy min-image in its xsimd ket-loop because member-load-through-`this` defeats compiler hoisting in that hot path (~7% per-step cost).

### Facade (CUDA-free)

```
core/        ← Real, Int, Mat3, units, constants, device tags, backend_traits
facade/      ← Driver<PhaseEnum, DeviceTag>, Plan, Stage (phase enums)
potentials/  ← backend-agnostic pair force/energy kernels (LJKernel, CoulKernel,
               SumKernel + width-agnostic math) shared by CPU/CUDA evaluators
molrs.hpp    ← stable CXX FFI header
```

### Driver taxonomy

CPU has four **outer task category** drivers, each subclassing
`CPUDriver<PhaseEnum>` which itself wraps the underlying
`Driver<PhaseEnum, DeviceTag>` template. The category is decided by
**termination semantics**, not by algorithm:

| Category base | Phase enum | Termination | Concrete CPU driver |
|---|---|---|---|
| `cpu::EvaluateDriver` | `evaluate::Phase` | exactly one pass | `cpu::Evaluate` |
| `cpu::DynamicsDriver` | `md::Phase` | full step budget; honours `outer_status.stop_requested` | `cpu::md::MD` |
| `cpu::OptimizeDriver` | `minimizer::Phase` | early stop on `outer_status.converged` | `cpu::minimizer::Minimize` |
| `cpu::ResponseDriver` | `vib::Phase` | fixed 6N displacement budget + one finalize pass; no early stop, no `outer_status` | `cpu::vib::Hessian` |

Adding a fifth user-facing category requires a genuinely distinct
outer-loop semantic (e.g. Sampling for MC/HMC/REMD, Path for NEB).
Never to wrap an algorithm under a new class name.

**SCF is *not* an outer category.** It is an inner solver represented
by `cpu::scf::SCFSolver` (a public low-level class, *not* a Driver
subclass), used internally by force-providing evaluator components
like `cpu::scf::SCFEvaluator`. Single-point energy in user code is
`cpu::Evaluate + cpu::scf::SCFEvaluator`. There is no `SinglePoint`
class. The CPU **plane-wave** SCF follows the identical shape:
`cpu::pw::PWSolver` (inner solver owning `internal::PWWorkspace`)
behind the `cpu::pw::PWEvaluator` Evaluator component
(`Capability::Energy`, writes `frame.energy` in eV; Ry/Ha on
`last_result()`); its step helpers (`PoissonSolver`, `XCBuilder`,
`VeffBuilder`, `DavidsonSolver`, `DensityBuilder`, `ConvergenceCheck`)
live in `cpu::pw::internal::` and are never Plan components. The
historical `ATV::scf::Phase` enum is deleted — no SCF anywhere carries
an outer phase enum (grep-gated by `test_scf_phase_scrub`). `pw_scf()`
survives as a thin shim over a transient `PWSolver`. Under
`ATV_HAS_MPI`, `PWSolver::run()` distributes k-points round-robin
across `MPI_COMM_WORLD` and Allreduces the density / band-energy
accumulations (parity gated by `test_mpi_pw_kpoint_n{1,2,4}`).

CUDA mirrors the same three-category taxonomy. `cuda::Evaluate`,
`cuda::md::MD`, and `cuda::minimizer::Minimize` are the three concrete
CUDA category drivers, each subclassing `cuda::CUDADriver<PhaseEnum>`
(the CUDA backend layer, parallel to `cpu::CPUDriver<PhaseEnum>`).
CUDA SCF is `cuda::scf::SCFSolver` — an inner solver, not a driver —
used by the `cuda::scf::SCFEvaluator` Evaluator component. Single-point
CUDA energy is `cuda::Evaluate + cuda::scf::SCFEvaluator` (or the
`cuda::single_point` helper); there is no `cuda::scf::SinglePoint`
class. The optimizer algorithm is a component:
`cuda::minimizer::SteepestDescent` plugs into `cuda::minimizer::Minimize`.

```
CPUDriver<PhaseEnum>::run(frame, ...)
  setup_pass(frame)
  for each iteration (semantic depends on category):
    reset outer_status
    backend.step(idx)
    consult outer_status (Optimize: break on converged; Dynamics: break on stop_requested)
  teardown_pass()
```

**Canonical facade.** `ATV::Driver<PhaseEnum, DeviceTag>` (in
`facade/driver.h`) is the underlying plan/dispatch primitive. User
code does not subclass it directly — use one of the three category
bases. Concrete drivers add only infrastructure components in their
constructor; they add no runtime behaviour beyond the base.

**Component tagging.** Two orthogonal axes (`facade/component_role.h`):
- `ComponentRole role()` — structural slot in the loop:
  `Evaluator`, `Integrator`, `Optimizer`, `Thermostat`, `Barostat`,
  `ConvergenceCheck`, `Recorder`, `Observer`, `None` (structural).
- `Capability capabilities()` — bitset of frame fields that
  `launch()` may write: `Energy`, `Force`, `Stress`, `Charge`,
  `Dipole`, `DensityMatrix`, `Orbitals`. Bitwise `|` / `&` /
  `has_capability(set, bit)` provided.

Algorithm classes are **components**, not driver subclasses:
`cpu::minimizer::LBFGS` / `CG` / future `FIRE` are Optimizer
components plugged into `Minimize`. `cpu::md::VelocityVerlet` /
`NoseHoover` / `FixLangevin` / `Barostat` are Integrator /
Thermostat / Barostat components plugged into `MD`.

**Iteration status channel** (`facade/iteration_status.h`).
`cpu::Context::outer_status` is a `{converged, stop_requested,
reason}` struct, owned by the category driver, reset before each
iteration. The Optimizer components write it during their
`convergence_check`-phase binding (LBFGS / CG own their stopping
criteria — there is no separate ConvergenceCheck-tagged component;
the `ComponentRole::ConvergenceCheck` tag is reserved for future
standalone checks such as a NaN/blow-up detector). SCF internals
never touch it (they manage their own internal convergence inside
`SCFSolver`).

**Rules:**
- **Driver owns the loop.** Force evaluators, integrators, thermostats are Components — never standalone classes with their own loops.
- **Four current outer driver categories** (Evaluate / Dynamics / Optimize / Response). Adding a fifth requires a distinct outer-loop semantic, never an algorithm wrapper.
- **Algorithm = component.** LBFGS/CG/FIRE are Optimizer components, never driver subclasses. Same for Thermostat / Barostat policies.
- **Capability over role for evaluator validation.** `EvaluateDriver` accepts any component with `role()==Evaluator` and non-empty `capabilities()`. `DynamicsDriver` requires any component (any role) with `Capability::Force`. `OptimizeDriver` requires any Energy/Force-capable evaluator + an Optimizer component.
- **No duplication across methods.** Shared algorithmic skeletons (gradient assembly, FD validation) are single parameterized functions.
- **Component lifecycle:** `build()` (one-time alloc) → `launch()` × N (per-step compute) → `destroy()` (cleanup). No allocation in `launch()`.
- **Two registration paths:** `driver.add(comp)` uses `default_bindings()`; `driver.on(phase).add(comp, order)` for explicit binding.
- **Observer vs Recorder:** single responsibility. Observer computes derived quantities. Recorder writes I/O. Never both.
- **Hot-path data structures carry only fundamental physical quantities.** `cpu::Frame`, `cpu::scf::IntegralPack`, `cpu::scf::ShellList`, etc. store positions, velocities, integrals, shell descriptors — never precomputed indices, caches, or derived maps. Consumers that need a derived view compute it themselves from the fundamentals.
- **SCF internals are not Plan components.** `cpu::scf::internal::*` (HamiltonianBuilder / EigenSolver / DensityBuilder / DIISMixer / DualConvergenceCheck) are owned by `SCFSolver` and never inherit `cpu::ComponentBase`. They live in `include/atomiverse/cpu/scf/internal/` and are not included from any facade-public header (tests may include them).

### Data residency

The three execution frames — `cpu::Frame` (xtensor on host), `cuda::Frame`
(raw device pointers), `torch::Frame` (`at::Tensor`) — carry an **identical
field set** (atoms incl. `charge`/`tag`/`image`, bonds/angles/dihedrals, the
`core/system_scalars.h` `SystemScalars` bundle, box); they differ only in
device representation. Per the Frame purity rule the frame holds **no**
method-specific scratch: the SCF workspace is owned by `*::scf::SCFSolver`
(`cuda::scf::SCFWorkspace`), the neighbor pair list by `*::NeighborList`
(`cuda::PairList`), published to consumers as a non-owning pointer on
`Context::pair_list`.

- `cpu::Frame` = xtensor arrays on host. Zero-copy in Driver (no transfer).
- `cuda::Frame` = device pointers. `cpu::to_cuda()` / `cuda::to_cpu()` for transfer.
- `torch::Frame` = `at::Tensor` (autograd-capable). `torch::transfer` bridges to cpu/cuda.
- `molrs.Frame` ⇄ a device frame via `cpu`/`cuda`/`torch::materialize()` (explicit H2D) and `write_back()` (explicit D2H) — backend roots, sibling of `transfer`, run once per `run()`.
- No implicit transfer in low-level operators. CUDA pipelines stay on device between explicit transfers.

### CUDA execution model

Components classify as: **compute** (device-only launch), **CPU bridge** (D2H + host callback via AsyncBridge), or **host-dependent compute** (documented exception). See `.claude/NOTES.md` for detailed classification tables.

### molrs integration

molrs (Rust) is a hard dependency for I/O. `cpu::MolRec` uses PIMPL — no `rust/cxx.h` in public headers. `molrs.hpp` is the only include for CXX-generated code.

## Conventions

- **Namespace:** `ATV`. Sub: `ATV::cpu`, `ATV::cuda`, `ATV::md`, `ATV::minimizer`, `ATV::scf`, `ATV::evaluate`. SCF internals: `ATV::cpu::scf::internal` and `ATV::cuda::scf::internal`. Inner SCF solvers: `ATV::cpu::scf::SCFSolver`, `ATV::cuda::scf::SCFSolver`.
- **File extensions:** CPU `.h`/`.cpp`, CUDA `.cuh`/`.cu`, facade/core `.h`. Documented exceptions: deliberately CUDA-free bridge headers under `cuda/` keep `.h` (`cuda/log_callback.h`, `cuda/timer.h`).
- **Driver template:** `Driver<PhaseEnum, DeviceTag>` is the underlying primitive. User-facing CPU code subclasses one of `cpu::EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver`, and CUDA code one of `cuda::EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver`, not `Driver<>` directly. CPU/CUDA driver subclass headers must include the matching `backend_traits.h`.
- **Component tagging:** `role()` returns `ComponentRole` (structural slot); `capabilities()` returns `Capability` bitset (frame fields written). Both default to `None`.
- **No factory functions** — components are classes with constructors
- **No fallback values** — lookup tables throw on unsupported input
- **Phase naming:** Evaluate: `prepare, evaluate, finalize`. MD: LAMMPS-style. Minimizer: `search_direction, line_search, force, convergence_check, end_of_step`. Vib: `setup, displace, force, accumulate, diagonalize, thermochem, end_of_step`. There is **no SCF phase enum** — gaussian (`cpu::scf::SCFSolver`, `cuda::scf::SCFSolver`) and plane-wave (`cpu::pw::PWSolver`) SCF are inner solvers with no outer phases (`ATV::scf::Phase` deleted; grep-gated).
- **Units:** eV, Å, fs, amu, eV/Å, K. Conversion constants in `core/units.h`.
- **Precision:** compile-time `Real` (`-DATV_REAL=float|double`), `Int` (`-DATV_INT=int|int64_t`)
- **libxc:** when enabled, native XC must match libxc to 1e-12
- **ERI performance:** `-march=native` is required for CPU SCF. xsimd 11.2.0 is linked directly to `atomiverse` for use in `native_integrals.cpp`. Do NOT set `XTENSOR_USE_XSIMD` as a compile definition globally — it activates SIMD-backed xtensor storage which breaks `reshape_view` in `xc_grid_engine.h`. `eri_block()` dispatches by **Ltot** (A.l+B.l+C.l+D.l) to a compile-time-Ltot template; ss|ss case uses xsimd ket-loop vectorization.
- **Doxygen:** every public symbol needs `@brief`, `@param` with units, `@return`

## Implemented components

### Drivers

| Driver | Backend | Category | Phase enum |
|--------|---------|----------|------------|
| `cpu::Evaluate` | CPU | `EvaluateDriver` | `evaluate::Phase` |
| `cpu::md::MD` | CPU | `DynamicsDriver` | `md::Phase` |
| `cpu::minimizer::Minimize` | CPU | `OptimizeDriver` | `minimizer::Phase` |
| `cpu::vib::Hessian` | CPU | `ResponseDriver` | `vib::Phase` |
| `cuda::Evaluate` | CUDA | `EvaluateDriver` | `evaluate::Phase` |
| `cuda::md::MD` | CUDA | `DynamicsDriver` | `md::Phase` |
| `cuda::minimizer::Minimize` | CUDA | `OptimizeDriver` | `minimizer::Phase` |

SCF is **not** a driver category on either backend. The inner SCF
solver is `cpu::scf::SCFSolver` / `cuda::scf::SCFSolver`, driven by the
`SCFEvaluator` Evaluator component. Single-point energy is
`Evaluate + scf::SCFEvaluator` — there is no `SinglePoint` class.

Single-point energy in user code:
```cpp
cpu::scf::SCFEvaluator scf{cfg};
scf.with_gradient(false);  // energy-only
cpu::single_point(frame, scf);
// frame.energy populated; scf.last_result().energy in Hartree.
```

### MD components
VelocityVerlet (Integrator), NVELimit (Integrator), InitVelocity, NoseHoover (Thermostat), NoseHooverNPT (Integrator — also acts as thermostat + barostat), FixLangevin (Thermostat, also Force-capable), Barostat (Barostat — Berendsen), Deform, NeighborList, PairLJCut (Evaluator: Energy|Force|Stress), BondHarmonic / AngleHarmonic / DihedralMultiHarmonic (Evaluator: Energy|Force), VirialStress. LAMMPS-style names are canonical — the 2026-05 `LJ` / `HarmonicBond` / `HarmonicAngle` / `MultiHarmonicDihedral` alias headers are deleted (torch backend keeps its own `LJ` / `HarmonicBond` names).

### Minimizer components
LBFGS (Optimizer — was `LBFGSStep`), CG (Optimizer — Polak-Ribière+),
OptLog / TrajectoryWriter / FinalGeometry (Recorder).
Optimizer choice is hot-swappable: the same `cpu::minimizer::Minimize`
driver accepts either by adding the component.

**MPI:** `Minimize::install_defaults()` mirrors MD's auto-injection,
folded into the `force` phase (Migrate @ -50, HaloSend @ -30 /
HaloRecv @ +10 when a pair evaluator is set, `md::mpi::Reduce` @ +20
so the optimizer reads globally-summed energy before
`convergence_check`). LBFGS / CG are DD-aware: owned-row loops,
globally-reduced inner products / step caps / convergence metrics, and
curvature-history reset on every `dd->neighbor_generation` bump
(migration reshuffles local rows). Recorders are DD-correct:
TrajectoryWriter / FinalGeometry gather owners to rank 0
(`md::mpi::gather_to_root`) before writing; OptLog reduces its metrics
across ranks.

### CUDA SCF (requires `-DATV_BUILD_RHF=ON`)
- **Inner solver:** `cuda::scf::SCFSolver` — public low-level class,
  NOT a `Driver<>` subclass. Owns the device-resident SCF cycle
  (Fock/Hamiltonian build → diagonalize → density build → convergence
  check), looping internally to convergence. Three run paths: synthetic,
  RHF from host integrals, RHF from GPU-native integral kernels.
- **Internal step helpers** (`cuda::scf::internal::`,
  `include/atomiverse/cuda/scf/internal/`): `Hamiltonian`, `RHFFock`,
  `RHFFockDevice`, `EigenSolver`, `DensityBuilder`, `RHFDensity`,
  `ConvergenceCheck`. Owned by `SCFSolver`, driven directly in a fixed
  order, carry no phase enum, NOT registered into any Plan, NOT
  `ComponentBase` subclasses.
- **Force evaluator:** `cuda::scf::SCFEvaluator` — `cuda::ComponentBase`
  with `role()==Evaluator` and `capabilities()=Energy|DensityMatrix|
  Orbitals`. Holds an `SCFSolver` member built once. Energy/density
  only — there is no CUDA SCF gradient kernel. Plugs into
  `cuda::Evaluate`.

### CPU SCF (requires `-DATV_BUILD_RHF=ON`)
- **Inner solver:** `cpu::scf::SCFSolver` — public low-level class.
  Owns the SCF cycle (Hamiltonian build → diagonalize → density →
  DIIS → convergence check). NOT a `Driver<>` subclass. Used by
  `SCFEvaluator` and available directly for low-level / debugging
  workflows.
- **Internal step helpers** (`cpu::scf::internal::`,
  `include/atomiverse/cpu/scf/internal/`): `HamiltonianBuilder`,
  `EigenSolver`, `DensityBuilder`, `DIISMixer`, `DualConvergenceCheck`.
  Owned by `SCFSolver`, NOT registered into any Plan, NOT
  ComponentBase subclasses. Tests may include the headers directly;
  facade-public production headers must not.
- **Force evaluator:** `cpu::scf::SCFEvaluator` — `cpu::ComponentBase`
  with `role()==Evaluator` and `capabilities()=Energy|DensityMatrix|
  Orbitals` (plus `Force` when `with_gradient(true)`). Holds an
  `SCFSolver` member built once. Plugs into `cpu::Evaluate`,
  `cpu::md::MD` (AIMD), and `cpu::minimizer::Minimize` (geometry
  optimization).
- **`scf_solve(integrals, config)`** — kept as a ≤30-line shim that
  constructs a transient `SCFSolver` and runs it. Used internally
  by `rhf_solve` / `uhf_solve` / `uks_solve`.
- **Gradients:** `scf_gradient(frame, integrals, result, config)` — unified analytic gradient (HF + KS)
- **Utilities:** `basis_eval.h`, `grid.h`, `xc_functional.h`, `aspc_predictor.h`

### CPU PW (requires `-DATV_BUILD_PW=ON`, FFTW3)
- **Inner solver:** `cpu::pw::PWSolver` — public low-level class, NOT a
  `Driver<>` subclass. Configured by `PWConfig` (upf_path, ecut_ry,
  k_mesh, tol, max_iter, mix_beta, pulay_depth); owns the
  facade-private `internal::PWWorkspace` and the SCF loop
  (Poisson → XC → Veff → per-k Davidson → density build + Pulay mix →
  convergence). `build()` parses the UPF; geometry-dependent setup
  (Ewald, G-vectors, FFT plans, V_loc, KB projectors) is rebuilt per
  `run()`, reusing grid-sized scratch when FFT dims are unchanged.
- **Internal step helpers** (`cpu::pw::internal::`,
  `include/atomiverse/cpu/pw/internal/`): `PoissonSolver`, `XCBuilder`,
  `VeffBuilder`, `DavidsonSolver`, `DensityBuilder`,
  `ConvergenceCheck` + numerics (`ewald`, `fft`, `g_vector`, `kpoint`,
  `upf_reader`, `pw_workspace`). Owned by `PWSolver`, NOT
  ComponentBase subclasses, NOT registered into any Plan; public
  headers must not include them (tests may).
- **Evaluator:** `cpu::pw::PWEvaluator` — `cpu::ComponentBase`,
  `role()==Evaluator`, `capabilities()==Energy` (no PW forces yet).
  Writes `frame.energy` [eV]; `last_result().energy` [Ry] /
  `.energy_ha` [Ha]. Plugs into `cpu::Evaluate` via
  `eval.on(evaluate::Phase::evaluate).add(pwe)`.
- **`pw_scf(frame, upf, ecut, k, tol, max_iter)`** — thin shim over a
  transient `PWSolver` (back-compat + scripting).
- **MPI:** k-point round-robin over `MPI_COMM_WORLD` inside `run()`
  (gated `ATV_HAS_MPI` + `MPI_Initialized`); density + band-energy
  Allreduced per iteration; Davidson initial-guess noise seeded per
  subspace index so results are rank-layout independent.
- **Shared problem assembly** (`internal/pw_setup.h`,
  `internal/pw_subspace.h`): `setup_pw_problem()` (Ewald, G-vectors,
  FFT grids, V_loc, per-k KB projectors / atomic wfc on q = k+G,
  initial density), `pulay_mix()`, and the Davidson subspace utilities
  (Y_lm, Jacobi, `init_subspace`) are reused verbatim by the CUDA
  backend for numerical parity.

### CUDA PW (requires `-DATV_BUILD_CUDA=ON` + `-DATV_BUILD_PW=ON`)
- **CUDA:** `cuda::pw::pw_scf(...)` (`cuda/pw/pw_scf.cuh`) — same
  contract as the CPU `pw_scf()` shim. The Davidson hot path (H·ψ via
  cuFFT Z2Z, KB projector products and subspace algebra via cuBLAS,
  per-band density accumulation) runs on device in
  `cuda::pw::internal::DeviceBandSolver`; ALL other steps (problem
  assembly, Poisson, LDA, V_eff, Pulay mixing, convergence, energy
  report) reuse the CPU code via the shared `cpu::pw::internal`
  helpers (`pw_setup.h`, `pw_subspace.h`) for numerical parity
  (`tests/cuda/pw/test_pw_si.cu` asserts CPU/CUDA energy agreement).
  Inner-solver contract: no Driver, no phase enum in the CUDA tree
  (`test_cuda_scf_phase_scrub` enforces).

### Output (cross-driver)
Thermo, Recorder (XYZ), RDF

## Example

User-facing CPU drivers expose a LAMMPS-shaped categorical setter API
(`set_pair`, `set_bond`, `set_neighbor`, `set_integrator`,
`set_thermostat`, `set_barostat`, `set_observer`, `set_record`,
`set_thermo`, `set_optimizer`). NeighborList is user-set (LAMMPS
`neighbor`); halo exchange auto-injects under `#ifdef ATV_BUILD_MPI`
when `set_pair` is called; `VirialStress` auto-injects when an
NPT / Berendsen barostat declares `reads()==Stress`. The
underlying `add()` and `on(phase).add()` are still public — use
them to register user-defined hooks the categorical API does not
yet have a slot for.

```cpp
// Single-point SCF energy + density (no forces)
cpu::scf::SCFEvaluator scf{cfg};
scf.with_gradient(false);
cpu::single_point(frame, scf);
double E_eV = frame.energy;
double E_Ha = scf.last_result().energy;

// Classical MD (Dynamics) — categorical setters
cpu::md::MD sim;
sim.set_pair(cpu::PairLJCut{0.0104, 3.4, 6.8});
sim.set_neighbor(cpu::md::NeighborList<>{6.8, /*skin=*/0.3});
sim.set_integrator(cpu::md::VelocityVerlet{0.01});
sim.set_thermo(100);
sim.build();
sim.run(frame, 1000);

// NPT MD — VirialStress auto-injected via reads()==Stress
cpu::md::MD npt;
npt.set_pair(cpu::PairLJCut{0.0104, 3.4, 6.8});
npt.set_neighbor(cpu::md::NeighborList<>{6.8, 0.3});
npt.set_integrator(cpu::md::NoseHooverNPT{
    {.T_target=300.0, .tau_T=100.0, .P_target=0.0, .tau_P=1000.0, .dt=0.5}});
npt.build();        // VirialStress injected at post_force @ -50 automatically
npt.run(frame, 1000);

// AIMD (Dynamics + SCFEvaluator as ForceProvider)
// cpu::scf::SCFEvaluator has no categorical setter and no
// default_bindings() (it serves three driver categories, so the phase
// cannot be a component-level default) — bind it explicitly with
// on(<driver's force phase>).add(...). A bare add() throws.
cpu::md::MD aimd;
aimd.set_integrator(cpu::md::VelocityVerlet{0.5});
aimd.on(md::Phase::force).add(cpu::scf::SCFEvaluator{rhf_config});
aimd.build();
aimd.run(frame, 100);

// Geometry optimization (same explicit binding rule as AIMD)
cpu::minimizer::Minimize opt;
opt.set_optimizer(cpu::minimizer::LBFGS{5, 0.01, 0.05});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
opt.build();
opt.run(frame, /*max_iter=*/100);

// CUDA single-point SCF energy — cuda::Evaluate + SCFEvaluator
cuda::Evaluate ceval;
ceval.add(cuda::scf::SCFEvaluator{integrals, tol, /*max_iter=*/100});
ceval.build();
ceval.run(frame);            // frame.scalars.energy populated
// or the free-function helper:
cuda::scf::SCFEvaluator cscf{integrals, tol, 100};
cuda::single_point(frame, cscf);

// CUDA classical MD (Dynamics)
cuda::md::MD csim;
csim.on(md::Phase::force)
    .add(cuda::NeighborList{6.8}, -10)
    .add(cuda::PairLJCut{0.0104, 3.4, 6.8});
csim.add(cuda::md::VelocityVerlet{0.01});
csim.build();
csim.run(frame, 1000);
```

## Skills & agents

General workflow skills (impl / fix / debug / refactor / review / docs / test / spec / litrev / note) come from the **molcrafts `mol` plugin** (`/mol:impl`, `/mol:fix`, …). They auto-adapt to this project via the `mol_project:` frontmatter at the top of this file.

**Atomiverse-specific** skills layer on top for domain audits the generic mol skills can't perform:

`/av-arch`, `/av-cuda`, `/av-hpc`, `/av-perf`

Each delegates to a paired Atomiverse-specific agent (`av-architect`, `av-cuda`, `av-hpc`, `av-optimizer`) that encodes Atomiverse-only knowledge: CUDA execution-model classification, xsimd 11.2.0 / xtensor 0.25.0 compatibility constraints, ERI shell-block dispatch + `-march=native` requirement, Frame purity rules, and hot-path anti-patterns specific to the dual-backend layout.

Knowledge flow: `.claude/NOTES.md` (evolving) → CLAUDE.md (stable) ← agents (read-only).
Specs live under `.claude/specs/` while in flight and are deleted at `done`
(completed-spec history is the git log).
