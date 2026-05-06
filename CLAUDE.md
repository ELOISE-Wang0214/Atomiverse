# CLAUDE.md

## Build & Test

```bash
cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_RHF=ON
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
./build/tests/test_md          # single test
./build/benchmarks/bench_md_step  # requires -DATV_BUILD_BENCHMARKS=ON
```

Requires: C++17, CMake ≥ 3.18, CUDA toolkit, Rust toolchain (for molrs), xtensor (auto-fetched).

## Architecture

Dual-backend (CPU + CUDA) molecular simulation platform. Three pillars: **Driver** (owns loop), **Component** (single-step operation), **Frame** (device-resident state).

### Backend directory layering

Each backend (`cpu/`, `cuda/`, `torch/`) follows the same two-tier layout:

```
<backend>/
  frame.h, context.h, component.h, backend.h, ...  ← shared infrastructure (root = common)
  molrec.h, transfer.h                              ← cross-domain utilities
  md/          ← domain: integrators, thermostats, neighbor lists
  scf/         ← domain: solvers, integrals, gradients, evaluators
  minimizer/   ← domain: SD, CG, L-BFGS
  observer/    ← domain: analysis (RDF, Mulliken, density)
  recorder/    ← domain: I/O (XYZ, Zarr, thermo output)
  potentials/  ← domain: force-field evaluators (LJ, bonds, angles)
  bonded/      ← domain: bonded interactions (CPU only)
```

**Rule:** Files at the backend root are **shared infrastructure** used by multiple domains. Domain subdirectories hold domain-specific code. When adding new shared code (data models, utilities, base classes), place it at the root — never inside a domain subdirectory.

### Facade (CUDA-free)

```
core/        ← Real, Int, Mat3, units, constants, device tags, backend_traits
facade/      ← Driver<PhaseEnum, DeviceTag>, Plan, Stage (phase enums)
molrs.hpp    ← stable CXX FFI header
```

### Driver taxonomy

CPU has three **outer task category** drivers, each subclassing
`CPUDriver<PhaseEnum>` which itself wraps the underlying
`Driver<PhaseEnum, DeviceTag>` template. The category is decided by
**termination semantics**, not by algorithm:

| Category base | Phase enum | Termination | Concrete CPU driver |
|---|---|---|---|
| `cpu::EvaluateDriver` | `evaluate::Phase` | exactly one pass | `cpu::Evaluate` |
| `cpu::DynamicsDriver` | `md::Phase` | full step budget; honours `outer_status.stop_requested` | `cpu::md::MD` |
| `cpu::OptimizeDriver` | `minimizer::Phase` | early stop on `outer_status.converged` | `cpu::minimizer::Minimize` |

Adding a fourth user-facing category requires a genuinely distinct
outer-loop semantic (e.g. Sampling for MC/HMC/REMD, Path for NEB,
Response for finite-difference Hessian). Never to wrap an algorithm
under a new class name.

**SCF is *not* an outer category.** It is an inner solver represented
by `cpu::scf::SCFSolver` (a public low-level class, *not* a Driver
subclass), used internally by force-providing evaluator components
like `cpu::scf::SCFEvaluator`. Single-point energy in user code is
`cpu::Evaluate + cpu::scf::SCFEvaluator`. There is no `SinglePoint`
class.

CUDA still uses the per-algorithm shape (`cuda::md::MD`,
`cuda::minimizer::CG`/`SD`, `cuda::scf::SCF`) — migration tracked in
`.claude/specs/cuda-driver-taxonomy.md`.

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
iteration. `ConvergenceCheck`-tagged components write it; SCF
internals never touch it (they manage their own internal convergence
inside `SCFSolver`).

**Rules:**
- **Driver owns the loop.** Force evaluators, integrators, thermostats are Components — never standalone classes with their own loops.
- **Three current outer driver categories.** Adding a fourth requires a distinct outer-loop semantic, never an algorithm wrapper.
- **Algorithm = component.** LBFGS/CG/FIRE are Optimizer components, never driver subclasses. Same for Thermostat / Barostat policies.
- **Capability over role for evaluator validation.** `EvaluateDriver` accepts any component with `role()==Evaluator` and non-empty `capabilities()`. `DynamicsDriver` requires any component (any role) with `Capability::Force`. `OptimizeDriver` requires any Energy/Force-capable evaluator + an Optimizer component.
- **No duplication across methods.** Shared algorithmic skeletons (gradient assembly, FD validation) are single parameterized functions.
- **Component lifecycle:** `build()` (one-time alloc) → `launch()` × N (per-step compute) → `destroy()` (cleanup). No allocation in `launch()`.
- **Two registration paths:** `driver.add(comp)` uses `default_bindings()`; `driver.on(phase).add(comp, order)` for explicit binding.
- **Observer vs Recorder:** single responsibility. Observer computes derived quantities. Recorder writes I/O. Never both.
- **Hot-path data structures carry only fundamental physical quantities.** `cpu::Frame`, `cpu::scf::IntegralPack`, `cpu::scf::ShellList`, etc. store positions, velocities, integrals, shell descriptors — never precomputed indices, caches, or derived maps. Consumers that need a derived view compute it themselves from the fundamentals.
- **SCF internals are not Plan components.** `cpu::scf::internal::*` (HamiltonianBuilder / EigenSolver / DensityBuilder / DIISMixer / DualConvergenceCheck) are owned by `SCFSolver` and never inherit `cpu::ComponentBase`. They live in `include/atomiverse/cpu/scf/internal/` and are not included from any facade-public header (tests may include them).

### Data residency

- `cpu::Frame` = xtensor arrays on host. Zero-copy in Driver (no transfer).
- `cuda::Frame` = device pointers. `cpu::to_cuda()` / `cuda::to_cpu()` for transfer.
- No implicit transfer in low-level operators. CUDA pipelines stay on device between explicit transfers.

### CUDA execution model

Components classify as: **compute** (device-only launch), **CPU bridge** (D2H + host callback via AsyncBridge), or **host-dependent compute** (documented exception). See `.claude/NOTES.md` for detailed classification tables.

### molrs integration

molrs (Rust) is a hard dependency for I/O. `cpu::MolRec` uses PIMPL — no `rust/cxx.h` in public headers. `molrs.hpp` is the only include for CXX-generated code.

## Conventions

- **Namespace:** `ATV`. Sub: `ATV::cpu`, `ATV::cuda`, `ATV::md`, `ATV::minimizer`, `ATV::scf`, `ATV::evaluate`. SCF internals: `ATV::cpu::scf::internal`.
- **File extensions:** CPU `.h`/`.cpp`, CUDA `.cuh`/`.cu`, facade/core `.h`
- **Driver template:** `Driver<PhaseEnum, DeviceTag>` is the underlying primitive. User-facing CPU code subclasses one of `cpu::EvaluateDriver` / `DynamicsDriver` / `OptimizeDriver`, not `Driver<>` directly. CPU driver subclass headers must include `cpu/backend_traits.h`.
- **Component tagging:** `role()` returns `ComponentRole` (structural slot); `capabilities()` returns `Capability` bitset (frame fields written). Both default to `None`.
- **No factory functions** — components are classes with constructors
- **No fallback values** — lookup tables throw on unsupported input
- **Phase naming:** Evaluate: `prepare, evaluate, finalize`. MD: LAMMPS-style. Minimizer: `search_direction, line_search, force, convergence_check, end_of_step`. SCF: `build_hamiltonian, diagonalize, build_density, convergence_check, end_of_step` (CUDA-only legacy; CPU SCF is an inner solver and uses no outer phase enum).
- **Units:** eV, Å, fs, amu, eV/Å, K. Conversion constants in `core/units.h`.
- **Precision:** compile-time `Real` (`-DATV_REAL=float|double`), `Int` (`-DATV_INT=int|int64_t`)
- **libxc:** when enabled, native XC must match libxc to 1e-12
- **ERI performance:** `-march=native` is required for CPU SCF. xsimd 11.2.0 is linked directly to `atomiverse` for use in `native_integrals.cpp`. Do NOT set `XTENSOR_USE_XSIMD` as a compile definition globally — it activates SIMD-backed xtensor storage which breaks `reshape_view` in `xc_grid_engine.h`. `eri_block()` dispatches by **Ltot** (A.l+B.l+C.l+D.l) to a compile-time-Ltot template; ss|ss case uses xsimd ket-loop vectorization. See `.claude/specs/eri-hpc-opt.md`.
- **Doxygen:** every public symbol needs `@brief`, `@param` with units, `@return`

## Implemented components

### Drivers

| Driver | Backend | Category | Phase enum |
|--------|---------|----------|------------|
| `cpu::Evaluate` | CPU | `EvaluateDriver` | `evaluate::Phase` |
| `cpu::md::MD` | CPU | `DynamicsDriver` | `md::Phase` |
| `cpu::minimizer::Minimize` | CPU | `OptimizeDriver` | `minimizer::Phase` |
| `md::MD` (CUDA) | CUDA | (legacy — pending migration) | `md::Phase` |
| `cuda::minimizer::SD` / `CG` | CUDA | (legacy — algorithm-as-driver) | `minimizer::Phase` |
| `cuda::scf::SCF` | CUDA | (legacy — pending migration) | `scf::Phase` |

CUDA migration to the same three-category taxonomy is tracked in
`.claude/specs/cuda-driver-taxonomy.md`.

Single-point energy in user code:
```cpp
cpu::scf::SCFEvaluator scf{cfg};
scf.with_gradient(false);  // energy-only
cpu::single_point(frame, scf);
// frame.energy populated; scf.last_result().energy in Hartree.
```

### MD components
VelocityVerlet (Integrator), NVELimit (Integrator), InitVelocity, NoseHoover (Thermostat), NoseHooverNPT (Integrator — also acts as thermostat + barostat), FixLangevin (Thermostat, also Force-capable), Barostat (Barostat — Berendsen), Deform, NeighborList, LJ (Evaluator: Energy|Force|Stress), HarmonicBond / HarmonicAngle / MultiHarmonicDihedral (Evaluator: Energy|Force), VirialStress

### Minimizer components
LBFGS (Optimizer — was `LBFGSStep`), CG (Optimizer — Polak-Ribière+).
Optimizer choice is hot-swappable: the same `cpu::minimizer::Minimize`
driver accepts either by adding the component.

### SCF components (CUDA)
RHFFock, RHFFockDevice, EigenSolver, RHFDensity, DensityBuilder, ConvergenceCheck, Hamiltonian

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

### Output (cross-driver)
Thermo, Recorder (XYZ), RDF

## Example

```cpp
// Single-point SCF energy + density (no forces)
cpu::scf::SCFEvaluator scf{cfg};
scf.with_gradient(false);
cpu::single_point(frame, scf);
double E_eV = frame.energy;
double E_Ha = scf.last_result().energy;

// Classical MD (Dynamics)
cpu::md::MD sim;
sim.add(cpu::md::VelocityVerlet{0.01});
sim.on(md::Phase::force).add(cpu::md::LJ{0.0104, 3.4, 6.8});
sim.build();
sim.run(frame, 1000);

// AIMD (Dynamics + SCFEvaluator as ForceProvider)
cpu::md::MD aimd;
aimd.add(cpu::md::VelocityVerlet{0.5});
aimd.on(md::Phase::force).add(cpu::scf::SCFEvaluator{rhf_config});
aimd.build();
aimd.run(frame, 100);

// Geometry optimization (Optimize + LBFGS or CG as Optimizer component)
cpu::minimizer::Minimize opt;
opt.add(cpu::minimizer::LBFGS{5, 0.01, 0.05});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
opt.build();
opt.run(frame, /*max_iter=*/100);
```

## Skills & agents

Skills delegate to expert agents. All reference CLAUDE.md as single source of truth.

`/av-impl`, `/av-fix`, `/av-debug`, `/av-refactor`, `/av-review`, `/av-arch`, `/av-perf`, `/av-docs`, `/av-test`, `/av-cuda`, `/av-hpc`, `/av-spec`, `/av-litrev`, `/av-note`

Knowledge flow: `.claude/NOTES.md` (evolving) → CLAUDE.md (stable) ← agents (read-only).
Specs: `.claude/specs/INDEX.md` with lifecycle tracking.
