# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test

```bash
# Configure (from repo root)
cmake -B build -DATV_BUILD_TESTS=ON -DATV_BUILD_BENCHMARKS=ON -DATV_BUILD_EXAMPLES=ON

# Build
cmake --build build -j$(nproc)

# Run all tests
ctest --test-dir build --output-on-failure

# Run a single test
./build/tests/test_md

# Run benchmarks
./build/benchmarks/bench_md_step

# Run example
./build/examples/lj

# Generate docs
cmake -B build/docs -DATV_DOCS_ONLY=ON
cmake --build build/docs --target docs
cmake --build build/docs --target docs-open

# Low-level fallback
doxygen Doxyfile   # output: docs/html/
```

Requires: CUDA toolkit, C++17, CMake >= 3.18, Rust toolchain (for molrs), xtensor (fetched automatically). GPU architecture auto-detected (`CMAKE_CUDA_ARCHITECTURES=native`). molrs is a hard dependency; Corrosion fetches automatically but requires `cargo` on PATH.

## Architecture

Atomiverse is a C++/CUDA molecular simulation platform with a **dual-backend architecture**: a unified facade API dispatches to fully decoupled CPU (xtensor) and CUDA (native device arrays) backends. MD, DFT, and MLIP are interchangeable modules sharing this common data model.

> **Backend characterization:** The CUDA backend is **CUDA-primary, not CUDA-only**. MD per-step physics (integrate, force, neighbor list) runs device-resident end-to-end. SCF iteration currently has host scalar dependencies (convergence norms, energy accumulation). Observer/recorder components are explicit **CPU bridge** sidecars that perform D2H + host analysis. See [CUDA execution model](#cuda-execution-model) for the full classification.

### Dual-backend design

```
Facade (Driver, Plan, Stage)    -- orchestration, CUDA-free headers
  |
  +-- cpu:: device             -- xtensor storage, .h/.cpp files
  |     cpu::Frame             -- host-side data (sub-structs: atoms, bonds, angles, dihedrals)
  |     cpu::to_cuda()         -- cpu::Frame -> cuda::Frame (allocate + H2D copy)
  |
  +-- cuda:: device            -- device pointers, .cuh/.cu files
        cuda::Frame            -- GPU-side data (sub-structs: atoms, bonds, angles, dihedrals)
        cuda::Context          -- CUDA streams (void* to keep headers CUDA-free)
        cuda::to_cpu()         -- cuda::Frame -> cpu::Frame (D2H copy)
        cuda::destroy()        -- free all device memory
```

**Core rules:**
- CPU and CUDA are both **devices** with symmetric naming (`cpu::to_cuda()`, `cuda::to_cpu()`)
- `cpu::Frame` and `cuda::Frame` are **separate, independent device types** -- no shared template, no shared base
- Only `to_*` transfer functions (no `from_*`) -- each direction has one function
- Data residency is a **type-level property**, not a runtime flag
- **No implicit transfer** in low-level operators -- transfer functions are always explicit
- CUDA pipelines stay **fully on device** between explicit transfer calls
- CPU backend has **no dependency on CUDA**; CUDA backend has **no dependency on xtensor**

### Directory layout

```
include/atomiverse/
  core/          -- shared types: Real, Int, Mat3, units, constants, Scheduling, PairList contract, GTO
  facade/        -- orchestration, CUDA-free: Driver, Plan, Stage (phase enums)
  molrs.hpp      -- stable CXX FFI header (includes generated bridge.h)
  cpu/           -- CPU device: Frame (xtensor), Context, transfer (to_cuda), MolRec (PIMPL)
    scf/         -- CPU SCF: integrals, RHF/UHF solvers
    observer/    -- CPU observers: RDF, Mulliken, density grid
    recorder/    -- CPU recorders: XYZ, ZarrWriter, trajectory, thermo
  cuda/          -- CUDA backend: Frame, Context, transfer, ComponentBase, Backend, Timer
    neighbor_list.cuh -- cross-driver NeighborList (cluster-based pair builder)
    potentials/  -- energy/force evaluators: LJ, HarmonicBond, HarmonicAngle, etc.
    md/          -- MD-specific: integrators, thermostats, barostat, deform, init
    minimizer/   -- minimizer components
    scf/         -- SCF components (.cuh headers)
    observer/    -- observer components
    recorder/    -- recorder components
    scf/         -- SCF components
  torch/         -- LibTorch backend (optional, gated by ATV_BUILD_TORCH)
    neighbor_list.cuh -- torch NeighborList (holds cuda::NeighborList)
    potentials/  -- differentiable operators: LJ, HarmonicBond
    transfer.cuh -- to_torch, to_cuda, to_cpu, from_cpu, destroy
    frame.cuh    -- torch::Frame (Tensor-based system state)
src/
  cpu/           -- CPU backend implementations (.cpp)
    scf/         -- SCF sources (integrals, solvers)
    observer/    -- observer sources (RDF, Mulliken, density)
    recorder/    -- recorder sources (XYZ, trajectory, thermo)
    molrec.cpp   -- MolRec PIMPL forwarding to molrs
  cuda/          -- CUDA backend implementations (.cu)
    neighbor_list.cu -- NeighborList + md/neighborlist/ sub-kernels
    potentials/  -- LJ, HarmonicBond, HarmonicAngle, MultiHarmonicDihedral sources
    md/          -- MD component sources (integrators, thermostats, init)
    minimizer/   -- minimizer sources
    scf/         -- SCF component sources
    observer/    -- observer sources
    recorder/    -- recorder sources
    scf/         -- SCF sources
    transfer.cu  -- cpu::to_cuda / cuda::to_cpu / cuda::destroy
    driver.cu    -- Backend + Driver template instantiation
  torch/         -- LibTorch backend sources
```

### molrs integration (CXX FFI)

molrs (Rust, sibling repo `../molrs`) is a **hard dependency** for I/O and analysis. Linked via Corrosion + CXX bridge.

```
Atomiverse (C++)
  |
  +-- molrs_cxx (CXX bridge target, auto-generated)
  |     |
  |     +-- molrs-cxxapi (Rust staticlib)
  |           +-- molrs-core (Frame, Block, XYZ writer, RDF)
  |
  +-- include/atomiverse/molrs.hpp   -- stable C++ header (includes generated bridge.h)
```

**Key rules:**
- `cpu/molrec.h` uses **PIMPL** — no `rust/cxx.h` in the header, safe for CUDA includers
- `molrs.hpp` is the only include for CXX-generated code; business code never includes bridge headers directly
- Precision matching: CMake sets Cargo feature `f64` when `ATV_REAL=double`; `build.rs` generates `bridge.rs` from `CXX_BRIDGE_SCHEMA` template with matching float type
- Bridge surface: only coordinate slices (`&[{F}]`) and per-atom fields use the precision-dependent type; scalars cross as `f64` always
- RDF: pure C++ histogram (`cpu::rdf_accumulate` / `cpu::rdf_write`), not migrated to molrs
- `ZarrWriter` calls `molrs::molrec_write_zarr()` via FFI
- All frame data lives in `molrs::MolRec` on the Rust side (no C++ mirror)
- Mulliken analysis: stays in C++ (depends on electronic structure context — basis sets, overlap matrix — not raw frame data)
- Thermo output: pure C++, never goes through molrs
- Density grid evaluation (`cpu::eval_density_grid`): stays in C++, Zarr output not wired

### Data flow

The user always works with `cpu::Frame`. `Driver::run(cpu::Frame&, int)` handles the full backend lifecycle internally.

**CUDA path:**
```
cpu::Frame (xtensor)
  --> Driver::run(frame, n_steps)
      |-- Backend::init()         -- create CUDA streams (cuda::Context)
      |-- Backend::materialize()  -- cpu::to_cuda(frame) --> cuda::Frame
      |-- Backend::build_components() -- one-time resource allocation
      |-- Backend::run(n_steps)   -- simulation loop (launch per step)
      |-- Backend::destroy_components() -- release build-time resources
      |-- Backend::retrieve()     -- cuda::to_cpu(cuda_frame, frame)
      |-- Backend::shutdown()     -- cuda::destroy + free streams
```

**CPU path:**
```
cpu::Frame (xtensor)
  --> Driver::run(frame, n_steps)
      |-- Backend::init()         -- create cpu::Context
      |-- Backend::build_components() -- one-time resource allocation (xtensor arrays)
      |-- Backend::run(n_steps)   -- simulation loop (launch per step, xtensor ops)
      |-- Backend::destroy_components() -- release build-time resources
```

### Component lifecycle

Components follow a strict **build / launch / destroy** lifecycle on both backends:

```
build(BuildContext)     -- one-time resource preparation
launch(LaunchContext)   -- per-step computation
destroy(DestroyContext) -- resource cleanup
```

- **build()**: Called once after Context and Frame exist. CUDA: may perform cudaMalloc, cudaMemcpy, cudaStreamSynchronize. CPU: may allocate xtensor arrays and precompute lookup tables.
- **launch()**: Called every simulation step. For **physics phases** (initial_integrate through post_force): compute-only (no allocation, no I/O). For **end_of_step**: async D2H (CUDA) or output (CPU) allowed for output components.
- **destroy()**: Called once after the simulation loop. Releases all resources allocated in `build()`. Executed in reverse Plan order. CUDA: cudaFree. CPU: xtensor deallocation.

### Data residency (hard rules)

**CUDA backend:**
- **cuda::Frame is all-device pointers + host-cached metadata.** Host-cached counts (`n_atoms()`, `n_bonds()`, etc.) are plain `Int` fields set at materialize time.
- **Zero D2H in physics phases.** Build-time D2H is allowed.
- **Async D2H allowed at end_of_step.** Output components (Thermo, Recorder, RDF) only.
- **DMA over compute for zeroing.** Use `cudaMemsetAsync` instead of custom zeroing kernels.

**CPU backend:**
- **cpu::Frame is xtensor arrays on host.** All data is host-resident, no transfer needed.
- **No allocation in launch().** Workspace arrays should be allocated in build() and reused.
- **xtensor-native operations for zeroing.** Use `xt::zeros_like` or `std::fill`, not manual loops.

### Component model

Components are OOP classes inheriting from `ComponentBase` (in `cuda/component.cuh`).

**Single-responsibility rule:** A component is either an **observer** (computes a derived quantity) or a **recorder** (writes data to I/O). Never both.

- **Observer**: computes and returns results (e.g. Mulliken charges → `vector<Real>`). Does not know about MolRec, file formats, or output destinations.
- **Recorder**: consumes data and writes it somewhere (e.g. XYZ file, Zarr store, stdout). Calls molrs or CPU I/O functions.
- **No `populate_*` glue functions** that mix observation with recording. The caller decides what to do with observer output.

**`.cu` / `.cpp` boundary rule:** `.cu` files handle CUDA operations only (kernel launch, D2H issue, stream management). CPU-side logic (xtensor wraps, analysis, I/O) lives in `.cpp` files. AsyncBridge's `on_data_ready()` should call a `.cpp` helper — the `.cu` file passes raw pointers, the `.cpp` file does xtensor wrapping internally.

**molrs delegation rule:** An observer/recorder is delegated to molrs **only if** its computation depends solely on raw per-frame data (positions, velocities, pair distances, etc.) — not on simulation-specific context (basis sets, overlap matrices, force field parameters). If it needs context from the computational method, it stays in C++.

Two registration paths:

```cpp
sim.add(md::VelocityVerlet{0.001});   // uses default_bindings()
sim.on(md::Phase::force)              // explicit phase override
   .add(NeighborList{3.0}, -10)
   .add(LJ{1.0, 1.0, 3.0});
```

- Components with natural phases override `default_bindings()`.
- Cross-driver components return empty `default_bindings()` and must be bound explicitly.
- **Infrastructure components** in `cuda/md/init.cuh` -- ZeroForces, ZeroEnergy, ZeroKinetic, InitVelocity. Each Driver auto-adds its own.
- **Output components** (Thermo, Recorder, RDF) are cross-driver and live in `cuda/observer/` or `cuda/recorder/`.

### No contract system

Operations within a stage slot execute **sequentially**. No data contracts, access validation, or auto-zero inference.

### CUDA execution model

Components in the CUDA backend fall into three categories. New components must declare which category they belong to.

#### Compute components (device-resident)

Per-step `launch()` is fully device-resident: kernel launches and device-to-device DMA only. No D2H, no host scalar reads, no stream synchronization, no allocation.

**Clean examples:** VelocityVerlet, LJ, NeighborList, HarmonicBond, NoseHoover, ZeroForces.

#### CPU bridge components (explicit host participation)

`launch()` performs D2H on the defer stream, synchronizes, then runs host-side analysis or I/O. These are **not** CUDA-only and must not be described as such. They are opt-in sidecars that do not affect compute-path claims.

**Examples:** Thermo, Recorder, RDF, MullikenObs, DensityRecorder.

**Pattern:** All bridges use `AsyncBridge` which is a deferred synchronous snapshot (not fire-and-forget):
1. `cudaEventRecord` on main stream
2. `cudaStreamWaitEvent` on defer stream
3. D2H memcpy on defer stream
4. `cudaStreamSynchronize(defer)` -- blocks host
5. Host callback (`on_data_ready`) -- calls a `.cpp` helper with raw pointers

**Hard rule:** The `.cu` file passes `Real*` + dimensions to the CPU helper. The CPU helper (`.cpp`) does xtensor wrapping, analysis, and I/O internally. **No xtensor includes in `.cu` files.**

#### Host-dependent compute (known violations)

These components participate in the main compute path but require host scalars or host-side math. They are documented exceptions, not design targets.

| Component | Host dependency | Location |
|-----------|----------------|----------|
| `ConvergenceCheck` | Host scalar norm + convergence decision | `src/cuda/scf/convergence_check.cu` |
| `RHFDensity` | Host scalar energy accumulation + writeback | `src/cuda/scf/rhf_density.cu` |
| `RHFFock` | CPU `S^{-1/2}` construction via xtensor-blas | `src/cuda/scf/rhf_fock.cu` (build) |
| `RHFFockDevice` | Host-side basis construction + eigenvalue D2H | `src/cuda/scf/rhf_fock.cu` (build) |
| `InitVelocity` | Host-side COM velocity correction | `src/cuda/md/init.cu` (build) |

**Target state:** SCF host scalar dependencies should move to device-side reductions. Build-time host work (InitVelocity, RHFFock integrals) is acceptable.

#### Synchronization discipline

| Sync point | Location | Justification |
|------------|----------|---------------|
| Per-slot dual-stream barrier | `Backend::step()` | Ensures phase ordering; over-synchronizes (collapses compute/defer overlap) |
| Defer-stream sync in AsyncBridge | `AsyncBridge::launch()` | Required for host callback correctness; blocks host per bridge invocation |
| Build-time stream sync | Various `build()` | Acceptable: one-time setup, not in hot path |
| SCF convergence sync | `ConvergenceCheck::launch()` | Required: host reads norm to decide loop continuation |

**Target state:** Replace per-slot barrier with dependency-driven sync. Allow bridge buffering to reduce defer-stream waits.

## Conventions

- **Namespace:** `ATV` (not `atomiverse`, not `AV`). Sub-namespaces: `ATV::cpu`, `ATV::cuda`, `ATV::md`, `ATV::minimizer`, `ATV::scf`
- **Backend as top-level, functionality mirrored underneath:** `cuda/potentials/`, `cuda/md/`, `cuda/minimizer/`, `cuda/scf/`, `cuda/observer/`, `cuda/recorder/` (same for `cpu/` and `torch/`). Cross-driver components (NeighborList) live at backend root (`cuda/neighbor_list.cuh`), not under a driver subdirectory.
- **File extensions:** CPU backend: `.h` + `.cpp`. CUDA backend: `.cuh` + `.cu`. Core/facade: `.h` (CUDA-free)
- **CMake registration:** new library sources go in root `CMakeLists.txt`; new tests use `atv_add_test(test_name)` in `tests/CMakeLists.txt`
- **No factory functions** -- components are classes with constructors
- **No `kernel_` prefix** on kernel functions
- **No fallback values** -- lookup tables (atomic radii, basis sets, element data) must cover all supported elements explicitly. Unsupported input **must** `throw std::invalid_argument` with a descriptive message, never silently return a default. This applies to Bragg-Slater radii, basis set parameters, and any element-indexed data.
- **libxc as XC ground truth** -- when `ATV_USE_LIBXC` is enabled, native XC functional output (exc, vrho, vsigma) must match libxc to 1e-12 at every test point. Any discrepancy is a bug in the native code, not a "numerical difference".
- **Acronyms in prose/comments:** always full uppercase (CPU, CUDA, XYZ, MD, DFT, MLIP, NVE, NVT, NPT, SoA, D2H, etc.), never title-case or lower-case
- **Phase naming:** LAMMPS-style for MD (`initial_integrate`, `pre_force`, `force`, `post_force`, `final_integrate`, `end_of_step`); minimizer phases (`search_direction`, `line_search`, `force`, `convergence_check`, `end_of_step`); SCF phases (`build_hamiltonian`, `diagonalize`, `build_density`, `convergence_check`, `end_of_step`)
- **Box:** `Mat3` (3x3 matrix), not separate scalars
- **Precision:** compile-time `Real` typedef set by CMake `-DATV_REAL=float|double`; compile-time `Int` typedef set by `-DATV_INT=int|int64_t`
- **Unit conversion:** `core/units.h` provides `kFtm2v` (force-to-acceleration, 1/103.6427), `kMvv2e` (mass-velocity^2-to-energy, 103.6427), `kBoltzmann` (8.617e-5 eV/K)
- **Doxygen:** every public symbol needs `@brief`, `@param` with units, `@return`, `@tparam`
- **Units:** eV (energy), Angstrom (distance), fs (time), amu (mass), eV/Angstrom (force), K (temperature)

## Key types

| Type | Layer | Header | Purpose |
|------|-------|--------|---------|
| `cpu::Frame` | CPU device | `cpu/frame.h` | System state: sub-structs (atoms, bonds, angles, dihedrals), Mat3 box, scalars |
| `cuda::Frame` | CUDA backend | `cuda/frame.cuh` | SoA device pointers in sub-structs (atoms, bonds, angles, dihedrals) + PairList + host-cached metadata |
| `cuda::Context` | CUDA backend | `cuda/context.cuh` | CUDA stream handles as `void*` (main + defer) |
| `cpu::Context` | CPU backend | `cpu/context.h` | Minimal CPU context (empty for now) |
| `BuildContext` | Lifecycle | `cuda/component.cuh` | Context for build(): Context + Frame (backend-specific) |
| `LaunchContext` | Lifecycle | `cuda/component.cuh` | Context for launch(): Context + Frame + phase + step |
| `DestroyContext` | Lifecycle | `cuda/component.cuh` | Context for destroy(): Context + Frame |
| `ComponentBase` | Backend | `cuda/component.cuh` | Abstract base: `name()`, `build()`, `launch()`, `destroy()`, `default_bindings()` |
| `Driver<PhaseEnum>` | Facade | `facade/driver.h` | `add()`, `on(phase).add()`, `build()`, `run(cpu::Frame&, n_steps)` |
| `Plan` | Planning | `facade/plan.h` | Compiled operation sequence (StageSlots sorted by order) |
| `Backend` | Execution | `cuda/backend.cuh` | Lifecycle manager + GPU dispatcher (internal) |
| `cpu::MolRec` | CPU I/O | `cpu/molrec.h` | Trajectory container (PIMPL over molrs::MolRec, no C++ mirror) |
| `cpu::XYZ` | CPU I/O | `cpu/recorder/xyz.h` | extxyz writer delegating to molrs |
| `cpu::ZarrWriter` | CPU I/O | `cpu/recorder/xyz.h` | Zarr V3 writer delegating to molrs::MolRec::write_zarr |
| `molrs::AtvMolRec` | molrs FFI | `molrs.hpp` | Opaque trajectory handle (Rust-side Frame/Block storage) |

## Implemented drivers and components

### MD (`md::MD`)
- Auto-adds: ZeroForces (pre_force, -1000), ZeroEnergy (pre_force, -999), ZeroKinetic (end_of_step, -1000), KineticEnergy (end_of_step, -100)
- **VelocityVerlet** -- symplectic integrator, auto-binds initial_integrate + final_integrate
- **InitVelocity** -- Maxwell-Boltzmann velocity generation + COM removal in build(), no-op launch()
- **NoseHoover** -- single-chain Nose-Hoover thermostat, auto-binds initial_integrate + final_integrate (order +100)
- **Barostat** -- Berendsen isotropic pressure coupling, auto-binds final_integrate (order +200)
- **Deform** -- uniaxial constant-strain-rate deformation, auto-binds initial_integrate (order -100)

### Minimizers
- **`minimizer::SD`** -- steepest descent, auto-adds ZeroForces/ZeroEnergy at force phase
- **`minimizer::CG`** -- conjugate gradient skeleton, auto-adds ZeroForces/ZeroEnergy at force phase

### Cross-driver components (empty default_bindings, bind via `.on().add()`)
- **NeighborList** (`cuda/neighbor_list.cuh`, `ATV::NeighborList`) -- cluster-based pair builder (morton sort, AABB, grid CSR, fused pair search, export)

### Potentials (`cuda/potentials/`, `ATV::` namespace)
- **LJ** -- Lennard-Jones 12-6 with shifted potential
- **HarmonicBond** -- harmonic bond stretching V=K(r-r0)^2
- **HarmonicAngle** -- harmonic angle bending V=K(theta-theta0)^2
- **MultiHarmonicDihedral** -- multi-harmonic dihedral V=sum(Cn cos^n(phi))

### Output/observer components (cross-driver, CPU bridge sidecars)
- **Thermo** -- thermodynamic output at end_of_step. CPU bridge: async D2H + host print. Supports `frequency()`.
- **Recorder** -- XYZ trajectory writer at end_of_step. CPU bridge: async D2H + `cpu::trajectory_append`. Supports `frequency()`.
- **RDF** -- radial distribution function at end_of_step. CPU bridge: D2H pair distances + `cpu::rdf_accumulate`.
- **VirialStress** -- virial stress tensor at post_force. Compute component (device-resident).

### SCF (`scf::SCF`)
- Auto-adds: ZeroDensity (build_hamiltonian, -1000)
- **Hamiltonian** -- toy SCF: H(P) = H0 + alpha*P, allocates workspace in build(), cuBLAS axpy
- **EigenSolver** -- dense symmetric eigendecomposition via cuSOLVER syevd (shared by toy + RHF)
- **DensityBuilder** -- toy SCF: P = C_occ * C_occ^T via cuBLAS gemm
- **ConvergenceCheck** -- ||P - P_old||_F via cuBLAS nrm2. Optional energy = dot(H, P) (flag `compute_energy`)
- **RHFFock** -- RHF: F = Hcore + G(P), F' = X^T F X. Uploads integrals in build() (requires `-DATV_BUILD_RHF=ON`)
- **RHFDensity** -- RHF: C = X*C', P = 2*C_occ*C_occ^T, E = 0.5*Tr[(Hcore+F)*P]+V_nuc

### Integral computation (RHF)
- `cpu::scf::build_integrals_sto3g(frame)` -- STO-3G integrals via **native C++** (S, T, V, Hcore, ERI, V_nuc)
- `cpu::scf::rhf_solve(integrals, tol, max_iter)` -- CPU RHF reference via **xtensor-blas**

Complete example:

```cpp
cpu::Frame frame(500);
// ... fill FCC positions, masses, box ...

// Phase 1: minimize
minimizer::SD min;
min.on(minimizer::Phase::force).add(NeighborList{6.8}, -10).add(LJ{0.0104, 3.4, 6.8});
min.on(minimizer::Phase::end_of_step).add(Thermo{100});
min.build();
min.run(frame, 500);

// Phase 2: NVE MD
md::MD sim;
sim.add(md::VelocityVerlet{0.01});
sim.add(md::InitVelocity{5.0, 42});
sim.on(md::Phase::force).add(NeighborList{6.8}, -10).add(LJ{0.0104, 3.4, 6.8});
sim.on(md::Phase::end_of_step).add(Thermo{100}).add(Recorder{"traj.xyz", 100}).add(RDF{10.0, 50, "rdf.dat"});
sim.build();
sim.run(frame, 1000);
sim.summary();

// Phase 3: SCF
scf::SCF solver;
solver.add(scf::Hamiltonian{8, 0.1});   // n_occ=8, alpha=0.1
solver.add(scf::EigenSolver{});
solver.add(scf::DensityBuilder{});
solver.add(scf::ConvergenceCheck{1e-5});
solver.on(scf::Phase::end_of_step).add(Thermo{5});
solver.build();
solver.run(frame, 50);

// Phase 4: RHF (requires -DATV_BUILD_RHF=ON)
cpu::Frame h2(2);
h2.atoms.x[0] = 0; h2.atoms.x[1] = 0.74; h2.atoms.type[0] = 1; h2.atoms.type[1] = 1; // H2
h2.box = Mat3::cubic(100.0);
auto integrals = cpu::scf::build_integrals_sto3g(h2);  // native C++
scf::SCF rhf;
rhf.add(scf::RHFFock{integrals});
rhf.add(scf::EigenSolver{});
rhf.add(scf::RHFDensity{integrals.n_occ()});
rhf.add(scf::ConvergenceCheck{1e-10, false});  // energy computed by RHFDensity
rhf.build();
rhf.run(h2, 100);  // E_total = -1.1168 Hartree
```

## Driver types

- **DFT is a Driver** (owns SCF loop/state/convergence), never an Evaluator component.
- Each driver type provides its own `Phase` enum ending with `COUNT`.
- Each driver type auto-adds its own infrastructure components at the appropriate phases.

## Skills & agents

Skills are thin user entry points that delegate to expert agents. All skills and agents reference CLAUDE.md as the single source of truth -- architecture rules are not duplicated.

**Implementation:** `/av-impl` (scope-aware: SMALL/MEDIUM/LARGE), `/av-fix` (quick fixes)
**Diagnosis:** `/av-debug` (read-only troubleshooting of build/test/runtime failures)
**Refactoring:** `/av-refactor` (restructure while maintaining invariants)
**Review:** `/av-review` (full review), `/av-arch`, `/av-perf`, `/av-docs`, `/av-test`
**CUDA compliance:** `/av-cuda` (CUDA execution model audit -- classifies components, flags host violations, validates sync discipline)
**CPU HPC:** `/av-hpc` (CPU backend parallelism and vectorization -- OpenMP, MPI, SIMD, cache optimization)
**Planning:** `/av-spec` (specification), `/av-litrev` (literature review)
**Notes:** `/av-note` (capture architectural decisions into `.claude/NOTES.md`)

### Architecture knowledge flow

```
.claude/NOTES.md   -- evolving: new decisions and corrections from development
       | (promote when stable)
CLAUDE.md          -- stable: proven architecture rules and conventions
       ^ (read by)
agents + skills    -- never duplicate rules, always reference the sources
```

### Specs

`/av-spec` persists structured specifications to `.claude/specs/{slug}.md` with lifecycle tracking (draft -> approved -> in-progress -> done). See `.claude/specs/INDEX.md` for the current spec list.
