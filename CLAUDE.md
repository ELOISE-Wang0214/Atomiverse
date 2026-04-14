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

### Driver / Component model

```
Driver<PhaseEnum, DeviceTag>::run(frame, n_steps)
  for each step:
    for each phase in PhaseEnum:
      for each Component bound to this phase (sorted by order):
        component->launch(ctx)
```

**Canonical facade.** `ATV::Driver<PhaseEnum, DeviceTag>` (in `facade/driver.h`) is the **single** user-facing class. It owns the component set, compiles the `Plan` in `build()`, and drives the per-step/per-phase dispatch in `run()`. Concrete subclasses (`cpu::md::MD`, `cuda::scf::SCF`, `cpu::minimizer::LBFGS`, …) only exist to inject infrastructure components via their constructor — they add no runtime behaviour beyond the template. There is **no** `Spec` / `Plan` / `Driver` three-way split in the current tree; proposals in `.claude/specs/rust-frame-lowering.md` are forward-looking and not in force.

**Rules:**
- **Driver owns the loop.** Force evaluators, integrators, thermostats are Components — never standalone classes with their own loops.
- **No duplication across methods.** Shared algorithmic skeletons (gradient assembly, FD validation, AIMD integration) are single parameterized functions. Method-specific logic is injected via config/enum dispatch.
- **Component lifecycle:** `build()` (one-time alloc) → `launch()` × N (per-step compute) → `destroy()` (cleanup). No allocation in `launch()`.
- **Two registration paths:** `driver.add(comp)` uses `default_bindings()`; `driver.on(phase).add(comp, order)` for explicit binding.
- **Infrastructure components** (ZeroForces, ZeroEnergy) are driver-specific — each driver type (MD, minimizer, SCF) has its own set with `default_bindings()` pointing to the correct phase.
- **Observer vs Recorder:** single responsibility. Observer computes derived quantities. Recorder writes I/O. Never both.
- **Hot-path data structures carry only fundamental physical quantities.** `cpu::Frame`, `cpu::scf::IntegralPack`, `cpu::scf::ShellList`, etc. store positions, velocities, integrals, shell descriptors — never precomputed indices, caches, or derived maps (e.g. `basis_to_atom`). Consumers that need a derived view (an observer computing a per-atom summary, a Fock builder that wants pair-list bounds) compute it themselves from the fundamentals. Keeps data layout auditable and prevents silent staleness.

### Data residency

- `cpu::Frame` = xtensor arrays on host. Zero-copy in Driver (no transfer).
- `cuda::Frame` = device pointers. `cpu::to_cuda()` / `cuda::to_cpu()` for transfer.
- No implicit transfer in low-level operators. CUDA pipelines stay on device between explicit transfers.

### CUDA execution model

Components classify as: **compute** (device-only launch), **CPU bridge** (D2H + host callback via AsyncBridge), or **host-dependent compute** (documented exception). See `.claude/NOTES.md` for detailed classification tables.

### molrs integration

molrs (Rust) is a hard dependency for I/O. `cpu::MolRec` uses PIMPL — no `rust/cxx.h` in public headers. `molrs.hpp` is the only include for CXX-generated code.

## Conventions

- **Namespace:** `ATV`. Sub: `ATV::cpu`, `ATV::cuda`, `ATV::md`, `ATV::minimizer`, `ATV::scf`
- **File extensions:** CPU `.h`/`.cpp`, CUDA `.cuh`/`.cu`, facade/core `.h`
- **Driver template:** `Driver<PhaseEnum, DeviceTag>` — always 2 explicit template args. CPU driver subclass headers must include `cpu/backend_traits.h`.
- **No factory functions** — components are classes with constructors
- **No fallback values** — lookup tables throw on unsupported input
- **Phase naming:** MD: LAMMPS-style. Minimizer: `search_direction, line_search, force, convergence_check, end_of_step`. SCF: `build_hamiltonian, diagonalize, build_density, convergence_check, end_of_step`.
- **Units:** eV, Å, fs, amu, eV/Å, K. Conversion constants in `core/units.h`.
- **Precision:** compile-time `Real` (`-DATV_REAL=float|double`), `Int` (`-DATV_INT=int|int64_t`)
- **libxc:** when enabled, native XC must match libxc to 1e-12
- **ERI performance:** `-march=native` is required for CPU SCF. xsimd 11.2.0 is linked directly to `atomiverse` for use in `native_integrals.cpp`. Do NOT set `XTENSOR_USE_XSIMD` as a compile definition globally — it activates SIMD-backed xtensor storage which breaks `reshape_view` in `xc_grid_engine.h`. `eri_block()` dispatches by **Ltot** (A.l+B.l+C.l+D.l) to a compile-time-Ltot template; ss|ss case uses xsimd ket-loop vectorization. See `.claude/specs/eri-hpc-opt.md`.
- **Doxygen:** every public symbol needs `@brief`, `@param` with units, `@return`

## Implemented components

### Drivers

| Driver | Backend | Phase enum | Auto-adds |
|--------|---------|------------|-----------|
| `md::MD` / `cpu::md::MD` | CUDA / CPU | `md::Phase` | ZeroForces, ZeroEnergy, ZeroKinetic, KineticEnergy |
| `minimizer::SD`, `minimizer::CG` | CUDA | `minimizer::Phase` | ZeroForces, ZeroEnergy |
| `cpu::minimizer::LBFGS` | CPU | `minimizer::Phase` | ZeroForces, ZeroEnergy |
| `scf::SCF` | CUDA | `scf::Phase` | ZeroDensity |

### MD components
VelocityVerlet, InitVelocity, NoseHoover, Barostat, Deform, NeighborList, LJ, HarmonicBond, HarmonicAngle, MultiHarmonicDihedral, VirialStress

### SCF components (CUDA)
RHFFock, RHFFockDevice, EigenSolver, RHFDensity, DensityBuilder, ConvergenceCheck, Hamiltonian

### CPU SCF (requires `-DATV_BUILD_RHF=ON`)
- **Solvers:** `scf_solve(integrals, config)` — unified RHF/UHF/UKS engine
- **Gradients:** `scf_gradient(frame, integrals, result, config)` — unified analytic gradient (HF + KS)
- **SCFEvaluator** — `cpu::ComponentBase` at `md::Phase::force` for AIMD. Single component for all methods. Optional ASPC predictor.
- **LBFGSStep** — `cpu::ComponentBase` for geometry optimization. Two-loop L-BFGS with step capping.
- **Utilities:** `basis_eval.h`, `grid.h`, `xc_functional.h`, `aspc_predictor.h`

### Output (cross-driver)
Thermo, Recorder (XYZ), RDF

## Example

```cpp
// Classical MD
cpu::md::MD sim;
sim.add(cpu::md::VelocityVerlet{0.01});
sim.on(md::Phase::force).add(cpu::md::LJ{0.0104, 3.4, 6.8});
sim.build();
sim.run(frame, 1000);

// AIMD (ab initio MD)
cpu::md::MD aimd;
aimd.add(cpu::md::VelocityVerlet{0.5});
aimd.on(md::Phase::force).add(cpu::scf::SCFEvaluator{rhf_config});
aimd.build();
aimd.run(frame, 100);

// Geometry optimization
cpu::minimizer::LBFGS opt;
opt.add(cpu::minimizer::LBFGSStep{5, 0.01, 0.05});
opt.on(minimizer::Phase::force).add(cpu::scf::SCFEvaluator{config});
opt.build();
opt.run(frame, 100);
```

## Skills & agents

Skills delegate to expert agents. All reference CLAUDE.md as single source of truth.

`/av-impl`, `/av-fix`, `/av-debug`, `/av-refactor`, `/av-review`, `/av-arch`, `/av-perf`, `/av-docs`, `/av-test`, `/av-cuda`, `/av-hpc`, `/av-spec`, `/av-litrev`, `/av-note`

Knowledge flow: `.claude/NOTES.md` (evolving) → CLAUDE.md (stable) ← agents (read-only).
Specs: `.claude/specs/INDEX.md` with lifecycle tracking.
