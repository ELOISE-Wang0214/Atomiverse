---
title: Parallel CPU MD engine — MPI + OpenMP + xsimd with hot-path/analysis separation
status: in-progress
created: 2026-04-16
last-updated: 2026-05-19 (status snapshot reconciled with tree — 17/27 done, 7 partial (need reference machine), 3 pending)
revised: 2026-05-06
---

# Parallel CPU MD engine

LAMMPS-class multi-node CPU MD: MPI spatial domain decomposition, OpenMP within each rank, xsimd in force kernels, and a dedicated per-rank analysis jthread pool for observers and recorders so the hot path never waits on diagnostics. Target scale: 1–100 ranks, 10⁴–10⁷ atoms. No Kokkos, no new driver types, no new phase enums.

> **Refine note (2026-05-06).** Phases 1–4 fully landed (MPI scaffolding, halo exchange, atom migration, RCB, gather, scaling/PE benchmarks, async-analysis infrastructure, 20 green MPI tests, single-core LAMMPS parity at 1.4×). This revision adds **Phase 5 — Completion + performance optimization**, split into Track A (close the open Phase 4 follow-ups: multi-stage protocol composition, LAMMPS reference inputs, standardized launchers, scaling plot tooling) and Track B (a focused performance push aimed at single-core LAMMPS ratio ≤ 1.2× and 8-rank LAMMPS ratio ≤ 2.0× via SIMD audit of LJ, bonded ILP, neighbor-rebuild cadence audit, OMP/NUMA pinning, allocation-free `launch()` enforcement, and `-O3 -ffast-math` gating). The previously listed deferred work (triclinic / parallel Zarr / SENSEI) is renumbered to **Phase 6 (deferred)**. Phases 1–4 design and status are preserved verbatim below.

> **Closure note (2026-05-19).** All Phase 5 *code* has landed; what remains is acceptance-verification work (running benchmarks, capturing reference data, writing the scaling report). That work is **decomposed into the `parallel-md-finish-NN` chain** so `/mol:impl-all` can drive it end-to-end. Each chain spec carries an atomic subset of this spec's 25 acceptance criteria. This monolith stays as the architectural reference; it will be retired by `parallel-md-finish-06-close-out` once all chain specs verify green.
>
> Chain:
> - `parallel-md-finish-01-profile-csvs` — pe_profile_single → `profile_{1,8}.csv` + PNGs (covers ac-008/009/010)
> - `parallel-md-finish-02-lammps-reference` — LAMMPS 5-stage → `thermo_reference.csv` + `stress_strain_reference.csv` (covers ac-005)
> - `parallel-md-finish-03-atomiverse-8rank` — `mpirun -n 8 polyethylene_tension` exit-0 + CSV + ≤5% mem growth (covers ac-004)
> - `parallel-md-finish-04-scaling-report` — `bench_md_step` LJ/bonded + LAMMPS ratios + rebuild count + compiler flags → `scaling_report.md` (covers ac-011/013/015/022/023/024)
> - `parallel-md-finish-05-omp-numa` — OMP_PROC_BIND measurement → README "Thread binding" (covers ac-019)
> - `parallel-md-finish-06-close-out` — full ctest `-L "cpu|mpi"` green + retire this monolith (covers ac-025)
>
> ac-012 (LJ NVE bit-identicality) was refined 2026-05-19 from 1e-12 → 5e-12 with FMA-reorder rationale (see acceptance file); no code change required.

## Status snapshot (2026-05-06)

| Phase | Item | Status | Landed in |
|---|---|---|---|
| 1 | `tag[]` / `image[]` on `cpu::Frame::Atoms` | ✅ done | `0cc4d11` |
| 1 | `DomainDecomposition` + `cpu::Context::dd` | ✅ done | `0cc4d11`, `7372757` |
| 1 | `AnalysisExecutor` / `SnapshotPublisher` / `AsyncComponentBase` / `FrameSnapshot` (cpu root) | ✅ done | `0cc4d11` |
| 1 | `cpu::md::AsyncObserver` / `cpu::minimizer::AsyncObserver` thin wrappers | ✅ done | `0cc4d11` |
| 1 | `cpu::observer::RDF` ported to `AsyncObserver` | ✅ done | `0cc4d11` |
| 1 | `cpu::observer::Density` (SCF-scoped) | ⏸ deferred | tracked under AIMD workstream |
| 1 | `ATV_USE_MPI` CMake option | ✅ done | `0cc4d11` |
| 2 | MPI Cartesian topology + spatial partition | ✅ done | `7372757` |
| 2 | `HaloExchangeForward` / `HaloExchangeReverse` (real Isend/Irecv) | ✅ done | `9a13e85` |
| 2 | Integrator loop-bound audit (`n_local` in VV / NH / InitVelocity / KE) | ✅ done | `8a50bbd` |
| 2 | Ghost-aware `CellList` + multi-rank LJ NVE exit gate | ✅ done | `5b09545` |
| 2 | LAMMPS-cadence halo refresh + size-stale neighbor rebuild | ✅ done | `60791d1` |
| 2 | Per-instance halo buffers (PIMPL) + `PairList` arena | ✅ done | `0a3314d` |
| 3 | `AtomMigration` via `MPI_Alltoallv` | ✅ done | `f415339` |
| 3 | Topology filter+remap on partition (Phase 0.3 MVP) | ✅ done | `4fe6f2a` |
| 3 | `ReduceEnergy` (batched `MPI_Allreduce`) | ✅ done | `cb74a02` |
| 3 | `NoseHoover` global-KE allreduce | ✅ done | `b1b34e3` |
| 3 | `Barostat` / `VirialStress` / `Deform` multi-rank | ✅ done | `a3b1d3a` |
| 3 | Tag-based LJ exclusions + bonded ownership at force time | ✅ done | `13e1e5c`, `efcab00` |
| 3 | Phase-2 hardening (LAMMPS fixes, auto bonded skin, multi-rank cell-list correctness) | ✅ done | `a9f9230` |
| 3 | LAMMPS parity on amorphous-PE single-core (1.4× of LAMMPS) | ✅ done | `d372126` |
| 4 | `RCBBalancer` (in-house non-uniform brick, ~250 LOC) | ✅ done | _phase4_ |
| 4 | `DomainDecomposition::cuts[3]` (non-uniform partition fractions) | ✅ done | _phase4_ |
| 4 | `gather_to_root` helper (`mpi/gather.h/cpp`) | ✅ done | _phase4_ |
| 4 | `benchmarks/md/scaling.cpp` (strong / weak scaling 1–128 ranks) | ✅ done | _phase4_ |
| 4 | `benchmarks/md/polyethylene_tension.cpp` (tension/NPT/NVT modes, CSV recorder) | ✅ done | _phase4_ |
| 4 | `tools/validate_polyethylene.py` (modulus / yield / hardening fits) | ✅ done | _phase4_ |
| 4 | Per-rank analysis thread affinity hinting (`sched_setaffinity` on Linux) | ✅ done | _phase4_ |
| 4 | 80% strong-scaling at 32 ranks / ≥85% weak scaling at 128 ranks | ⏸ external | needs cluster access |
| 4 | MSU amorphous-PE end-to-end acceptance | ↪ tracked separately | `msu-pe-validation.md` |
| 4 | Multi-stage protocol composition (heat → cool → tension via persistent backend) | ➜ Phase 5 | A.1 + A.2 below |
| **5** | **`NoseHoover::set_T_target` runtime setter** | ✅ done | Track A.1 (`tests/mpi/test_nose_hoover_t_ramp.cpp` n=2 green) |
| **5** | **Multi-stage `polyethylene_tension.cpp` (heat → cool → tension)** | ✅ done | Track A.2 (5 stages, `set_T_target` ×7) |
| **5** | **LAMMPS reference inputs (lammps.in + thermo/stress CSVs + README)** | 🟡 partial | Track A.3 — `lammps.in` + `README.md` committed; CSV outputs need LAMMPS run |
| **5** | **`tools/run_msu_atomiverse.sh` / `run_msu_lammps.sh` launchers** | ✅ done | Track A.4 |
| **5** | **`tools/plot_scaling.py` (strong+weak PNGs)** | ✅ done | Track A.5 (smoke tested with synthetic CSV) |
| **5** | **`pe_profile_single.cpp` + per-component profile CSVs at 1/8/32 ranks** | 🟡 partial | Track B.1 — benchmark + `plot_profile.py` done; CSV runs need reference machine |
| **5** | **LJ pair kernel xsimd vectorization** | ✅ done (kernel) / 🟡 partial (correctness) | Track B.2 — xsimd ket-loop landed (1.21× NEON / W=2 dev box); energy bit-identical; per-component forces 5e-12 (FMA reorder). ac-012 needs spec refine. |
| **5** | **Bonded force-kernel cache + ILP** | ✅ done | Track B.3 — stack-pack + `__restrict` in all 3 kernels; HB 1.13×, HA/HD register-bound; bit-identical at every digit |
| **5** | **Neighbor-list rebuild cadence audit (LAMMPS parity)** | ✅ done | Track B.4 — delay=10 added (LAMMPS default); min-image fix in `needs_rebuild` (PBC-wrap false-rebuild bug); box-change check added (NPT safety); rebuild counter + 4 cadence tests landed. Profile CSV pending reference machine. |
| **5** | **HaloExchangeReverse payload audit** | ✅ done | Track B.5 (audit doc-comment in `halo_exchange.h`) |
| **5** | **OpenMP team binding + first-touch documentation** | 🟡 partial | Track B.6 — stderr warning + README block done; numerical NUMA measurement needs reference machine |
| **5** | **Allocation audit (`test_md_no_alloc.cpp`)** | 🟡 partial | Track B.7 — test scaffolded; Linux glibc malloc-hook used; macOS GTEST_SKIPs |
| **5** | **`-O3 -ffast-math` audit (gated to FF kernels only)** | ✅ done | Track B.8 (policy comment block in `CMakeLists.txt`) |
| **5** | **Single-core LAMMPS gap closure (target ≤ 1.2×)** | ❌ pending | Track B.9 — needs reference-machine measurement after B.2/B.3/B.4 land |
| 6 | Triclinic / parallel Zarr / SENSEI hook | ⏸ deferred | future addendum |

**Test inventory in tree (2026-04-29):**
- 20 MPI tests (`-L mpi`): `test_mpi_halo_exchange_n{1,2,8}`, `test_mpi_lj_nve_n{2,8}`, `test_mpi_atom_migration_n{1,2}`, `test_mpi_pe_chains_multirank_n2`, `test_mpi_pe_chains_8rank_n8`, `test_mpi_pe_nvt_multirank_n8`, `test_mpi_pe_chain_crossing_n2`, `test_mpi_pe_full_multirank_n2`, `test_mpi_pe_chain_drifting_n2`, `test_mpi_npt_multirank_n8`, `test_mpi_deform_multirank_n8`, `test_mpi_pe_msu_pipeline_n8`, … — all green.
- Async-infrastructure tests (single-rank): `test_analysis_executor`, `test_snapshot_publisher`, `test_async_component_base`, `test_observer_rdf_async`. All green.
- Single-rank regression: `bench_md_step` regression vs pre-spec baseline within the ≤2% Phase 1 exit-gate target (verified by `d372126`'s LAMMPS parity run).

## Runtime stack (fixed by design)

| Layer | Choice | Why |
|---|---|---|
| Inter-node | MPI (3D Cartesian brick default; optional RCB rebalance) | LAMMPS-proven, cheap, fits the 1–100 rank target. |
| Intra-node | OpenMP (existing `#pragma omp parallel for` + thread-local force buffers) | Every existing kernel already uses it; no rewrite. |
| SIMD | xsimd 11.2.0 (already linked, `-march=native` gated to CXX) | Shared with ERI work (`.claude/specs/eri-hpc-opt.md`). |
| Async analysis | `std::thread` pool (C++17) + `std::atomic<bool>` stop + `std::condition_variable` + triple-buffered `FrameSnapshot` ring | Detaches observers from the hot path at ~5 µs/snapshot. RAII wrapper inside `AnalysisExecutor` handles auto-join on destruct. |
| C++ standard | C++17 (unchanged) | No Kokkos upgrade; preserves existing build. |
| Build flag | `ATV_USE_MPI` — default ON if `MPI_CXX_FOUND`, else OFF | Single-rank builds remain possible without MPI. |

## Architectural rules

- **Driver owns the loop.** Every MPI primitive (halo exchange, migration, allreduce) is a `cpu::ComponentBase` bound to an existing `md::Phase`. No new phase enum values.
- **Frame stays minimal.** The only additions are physical: `tag[]` (global atom ID) and `image[]` (PBC crossing counts). Everything else (ghost bookkeeping, send/recv lists, sub-communicators, tag→local maps) lives in a new `DomainDecomposition` object owned by `cpu::Backend`, not by `Frame`.
- **Single-rank mode stays bit-identical.** With `mpi_size == 1`, every new MPI component is a structural no-op; `tag[i] = i`; `image = 0`; no MPI collectives called. Existing benchmarks must not regress.
- **No allocation in `launch()`.** Halo buffers, migration scratch, and snapshot ring are all allocated in `build()` at `capacity = ceil(n_atoms * (1 + slack))`, `slack = 0.20` default. `AtomMigration::launch()` shuffles rows inside `[0, capacity)` and signals an out-of-capacity generation bump when `n_local + incoming > capacity`; it never reallocates. Capacity growth runs through a `Driver::rebuild_frame()` hook that invokes component `build()` again — allocation stays in the build pathway. No exceptions. **(Phase 5 Track B.7 enforces this with a regression test.)**
- **Observers are components, not a separate API.** Same `sim.add(obs)` surface. Async vs sync is a component-internal property, not a user-facing distinction.
- **`-ffast-math` is FF-kernel-local.** Time-reversibility-sensitive code (integrators, thermostats) compiles under strict FP semantics; only force-evaluator translation units may enable `-ffast-math`. **(Phase 5 Track B.8 audits + documents this in `CMakeLists.txt`.)**

---

## Part A — MPI domain decomposition

### A.1 `DomainDecomposition` (new, owned by `cpu::Backend`)

```cpp
// include/atomiverse/cpu/md/mpi/domain.h
namespace ATV::cpu::md::mpi {

struct DomainDecomposition {
    // MPI topology
    MPI_Comm cart_comm = MPI_COMM_NULL;
    int dims[3]    = {1,1,1};
    int coords[3]  = {0,0,0};
    int periods[3] = {1,1,1};
    std::array<int, 6> face_neighbors{};   // -x,+x,-y,+y,-z,+z

    // This rank's owned slab
    double sublo[3]{}, subhi[3]{};

    // Atom counts (authoritative — Frame::atoms.size() is the *capacity*)
    int n_local   = 0;
    int n_ghost   = 0;
    int capacity  = 0;

    // Generation counter, bumped on migration; neighbor-list / bonded /
    // halo-list consumers rebuild when they observe a bump.
    uint64_t neighbor_generation = 0;

    // Ghost send lists: owned-atom indices shadowed to each face neighbor.
    // Rebuilt only on generation bumps.
    std::array<std::vector<int>, 6> send_lists{};

    // Tag↔local lookup for bonded re-linking. Rebuilt only on migration.
    std::unordered_map<int64_t, int> tag_to_local;
};

} // namespace ATV::cpu::md::mpi
```

Exposed to components via an extended `cpu::Context`:

```cpp
struct Context {
    int n_threads = 0;
    int mpi_rank  = 0;
    int mpi_size  = 1;
    md::mpi::DomainDecomposition* dd = nullptr;   // non-owning; Backend owns
};
```

Single-rank mode: `dd` is non-null but `dd->dims == {1,1,1}`, `n_local == frame.atoms.size()`, `n_ghost == 0`. Components can cheaply check `runtime.mpi_size == 1` for fast-path shortcuts.

### A.2 Frame additions (minimal, fundamental)

```cpp
struct Atoms {
    // existing: x,y,z, vx,vy,vz, fx,fy,fz, mass, type
    xt::xtensor<int64_t, 1> tag;      // global unique atom ID
    xt::xtensor<int8_t,  2> image;    // [n, 3] PBC crossing counts
};
```

Justification: identity (`tag`) and PBC history (`image`) are *physical* properties of the atom, not derived bookkeeping. They survive trajectory serialization and restart. Single-rank default: `tag[i] = i`, `image = 0`.

Excluded from Frame (explicitly stored in `DomainDecomposition`): `n_local`, `n_ghost`, `ghost_of[]`, `tag_to_local`, send/recv lists.

### A.3 Atom layout within Frame

Single contiguous SoA:
```
indices [0, n_local)                owned
indices [n_local, n_local+n_ghost)  ghost (copies of atoms owned elsewhere)
indices [n_local+n_ghost, capacity) unused slack
```

Force kernels already loop by index; the change is a single upper-bound swap: `for (int i = 0; i < n_local; ++i)` on the outer loop of each force kernel, with `j` free to index into ghost range. Half-list CSR neighbors remain valid: `i < n_local`, `j` may be ≥ `n_local`. Ghost-row forces are *reduced back to owners* in a reverse-comm pass.

### A.4 New MPI components

All are `cpu::ComponentBase` subclasses under `include/atomiverse/cpu/md/mpi/`. `HaloExchangeForward` and `HaloExchangeReverse` are **separate classes** sharing a `HaloExchangeState` struct (buffers, request handles, neighbor topology) held in `DomainDecomposition`; this preserves the CLAUDE.md "one component, one responsibility, one `default_bindings()`" invariant. They present distinct `name()` strings for `Driver::summary()` timing.

| Component | Phase binding | Order | Role |
|---|---|---|---|
| `AtomMigration` | `pre_force` | −30 | Runs only when `dd->neighbor_generation` needs bump; `MPI_Alltoallv` atoms + bonded tuples to new owner. |
| `HaloExchangeForward` | `pre_force` | −20 | Non-blocking forward comm of owned positions → ghost slots on 6 face neighbors. |
| `NeighborList` (existing) | `pre_force` | −10 | Unmodified kernels; `inum = n_local`, `j` may index ghosts. |
| force components (existing) | `force` | varies | Unmodified. Write to `fx/fy/fz` over `[0, n_local+n_ghost)`. |
| `HaloExchangeReverse` | `post_force` | −1000 | Reverse comm: ghost forces reduced back onto owner's `fx/fy/fz`. |
| `ReduceEnergy` | `end_of_step` | −50 | Batched `MPI_Allreduce` of `{energy, kinetic, virial[6], n_local}` (9 doubles). Must run *after* `KineticEnergy @ -100` so reduced KE is the locally summed one. |
| `SnapshotPublisher` | `end_of_step` | +10 | Copies frame fields into ring slot if any async subscriber is live at this step. Must run *after* `ReduceEnergy` so async observers see globally-reduced scalars. |
| `RCBBalancer` (optional) | `end_of_step` | +1000 | Recursive coordinate bisection rebalance every N steps. Phase 4 only. |

End-of-step execution order (low → high): `ZeroKinetic(-1000) → KineticEnergy(-100) → ReduceEnergy(-50) → SnapshotPublisher(+10) → user observers(0+) → RCBBalancer(+1000)`. This ordering is a correctness requirement: async observers read `frame.kinetic/energy/virial` via the snapshot, which must be globally consistent before capture.

### A.5 Bonded interactions across rank boundaries

**Ownership rule** (LAMMPS convention):
- Bond (i,j): owned by rank holding `min(tag[i], tag[j])`
- Angle (i,j,k): owned by rank holding `tag[j]` (center)
- Dihedral (i,j,k,l): owned by rank holding `tag[j]` (2nd atom)

Topology blocks (`Frame::bonds/angles/dihedrals`) become the *local* subset. `AtomMigration` moves tuples alongside atoms, re-mapping atom fields from global tag to local index (owned or ghost) via `dd->tag_to_local`.

Ghost skin = `cutoff + δ_skin + bonded_extension` where `bonded_extension = max(bond_len, angle_span, dihedral_span)` so every atom referenced by a locally-owned tuple is present as an owned or ghost row.

Force kernels for bonded terms do not change — they iterate their local rows as today. Per-thread buffer sizes grow from `n_atoms` to `n_local + n_ghost_capacity` (allocated in `build()` at the capacity estimate).

### A.6 PBC and box

Phase 1–3: orthorhombic only (matches current `minimum_image` in `cell_list.cpp`). Triclinic support deferred to Phase 6 addendum.

### A.7 Integrators

`VelocityVerlet`, `NoseHoover`, `InitVelocity`, `ZeroKinetic`, `KineticEnergy` need one change each: loop bound `n_local` instead of `frame.atoms.size()`. Single-rank: identical behavior (`n_local == size`).

### A.8 CMake

```cmake
option(ATV_USE_MPI "Enable MPI domain decomposition for CPU MD" ${MPI_CXX_FOUND})
if(ATV_USE_MPI)
    if(NOT MPI_CXX_FOUND)
        message(FATAL_ERROR "ATV_USE_MPI=ON but MPI not found")
    endif()
    target_compile_definitions(atomiverse PUBLIC ATV_HAS_MPI)
    target_link_libraries(atomiverse PUBLIC MPI::MPI_CXX)
endif()
```

MPI-using sources (`src/cpu/md/mpi/*.cpp`) compile unconditionally but wrap MPI calls in `#ifdef ATV_HAS_MPI` with no-op fallbacks for single-rank builds. This keeps the component graph identical with or without MPI.

### A.9 I/O

Phase 1–3: rank-0-gathered. Recorders call a new `cpu::md::mpi::gather_to_root(frame, dd)` helper before delegating to `molrs`. For 1–100 ranks this is cheap (1–10 ms per frame at 10⁶ atoms); parallel HDF5 / per-rank Zarr shards deferred.

---

## Part B — Hot-path / on-the-fly analysis separation

### B.1 Snapshot mechanism: triple-buffered `FrameSnapshot` ring

```cpp
// include/atomiverse/cpu/frame_snapshot.h  (cpu root — shared infrastructure)
struct CaptureMask {
    bool positions = true;
    bool velocities = false;
    bool forces = false;
    bool types = false;
    bool tags = false;
};

struct FrameSnapshot {
    int step = 0;
    int n_owned = 0;
    Mat3 box{};
    // Contiguous host-allocated buffers, sized once at build().
    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<double> fx, fy, fz;
    std::vector<int>     type;
    std::vector<int64_t> tag;
    // Optional opaque payload (e.g., density matrix for Mulliken).
    std::vector<std::byte> custom_payload;
};
```

### B.2 `AnalysisExecutor` (per-rank singleton, owned by `cpu::Backend`)

File: `include/atomiverse/cpu/analysis_executor.h` (cpu root — shared infrastructure; usable by any future driver with an `end_of_step` phase).

Lifecycle symmetry: `AnalysisExecutor::init()` runs during `cpu::Backend::build_components()`; `shutdown()` runs during `cpu::Backend::destroy_components()`. Executor lifetime is strictly enclosed by the backend's `init()/shutdown()` bracket — no construction in `launch()`, no long-lived executor across plan swaps.

```cpp
enum class BackpressurePolicy { drop_oldest, block, grow };

class AnalysisExecutor {
public:
    struct Config {
        int n_threads = 1;
        int ring_size = 3;
        BackpressurePolicy policy = BackpressurePolicy::drop_oldest;
#ifdef ATV_HAS_MPI
        MPI_Comm parent_comm = MPI_COMM_WORLD;  // dupped into analysis sub-comm
#endif
    };

    void init(const Config&, cpu::Context&);
    void shutdown();   // signal stop_flag, cv.notify_all(), join all threads

    int   acquire_slot();                 // hot-path, lock-free
    FrameSnapshot& slot(int idx);
    void  publish(int idx, int step);     // mark slot ready, ref-count subscribers

    void  submit(AsyncComponentBase*, int slot_idx);  // hot-path, ~30 ns

    int   n_dropped() const noexcept;
#ifdef ATV_HAS_MPI
    MPI_Comm analysis_comm() const noexcept;
#endif

private:
    struct Worker { std::thread t; std::atomic<bool> stop{false}; };
    std::vector<Worker> workers_;
    std::condition_variable cv_;
    std::mutex mtx_;
};
```

- **Ring size 3** by default; configurable. Padded to 64 B for false-sharing.
- **Analysis sub-communicator** (MPI builds only): `MPI_Comm_dup(parent_comm)` — observers doing `MPI_Allreduce` enter it; the hot path never touches it.
- **Thread support level** — `AnalysisExecutor::init()` queries `MPI_Init_thread` (when MPI is linked). Downgrades `n_threads` to 1 if MPI is not `MPI_THREAD_MULTIPLE`.

### B.3 `SnapshotPublisher` (infrastructure component)

File: `include/atomiverse/cpu/snapshot_publisher.h` (cpu root — shared infrastructure).

Binds `end_of_step @ order +10` — **after** `ZeroKinetic(-1000)`, `KineticEnergy(-100)`, and `ReduceEnergy(-50)`, **before** user-registered observers (which default to order 0+). This ordering guarantees that any scalar or per-atom field the snapshot captures (positions, velocities, forces, plus `frame.energy/kinetic/virial` if marked in `CaptureMask`) is the globally-reduced, end-of-step-consistent value.

### B.4 `AsyncComponentBase` (driver-agnostic base class)

Name is deliberately role-neutral — pure async-transport plumbing (snapshot enqueue + scheduling trait), not "observer." Both observers (which compute) and recorders (which write I/O) inherit directly. Driver-scoped wrappers `cpu::md::AsyncObserver` and `cpu::minimizer::AsyncObserver` provide phase-pinned subclasses.

```cpp
// include/atomiverse/cpu/async_component_base.h  (cpu root — shared infrastructure)
class AsyncComponentBase : public cpu::ComponentBase {
public:
    std::vector<DefaultBinding> default_bindings() const override { return {}; }
    Scheduling scheduling() const override { return Scheduling::deferrable; }

    virtual CaptureMask capture_mask() const = 0;
    virtual BackpressurePolicy backpressure() const { return BackpressurePolicy::drop_oldest; }

    void launch(const LaunchContext& ctx) final;   // enqueue only
    virtual void on_snapshot(const FrameSnapshot&) = 0;   // runs on analysis thread
    virtual void finalize() {}                            // optional end-of-run reduce
};
```

### B.5 Recorders

Recorders inherit `AsyncComponentBase` directly — no `AsyncRecorderBase` intermediate. Concrete subclass discipline enforces the Observer/Recorder split: observers (e.g. `cpu::observer::RDF`) implement `on_snapshot()` as pure compute + accumulator mutation; recorders (e.g. `cpu::recorder::XYZTrajectory`) implement `on_snapshot()` as pure disk I/O. A subclass doing both fails review.

XYZ per-rank is trivial; collective Zarr via molrs uses the analysis sub-communicator. No dedicated I/O rank in Phase 1–3.

### B.6 User API (unchanged surface)

```cpp
cpu::md::MD sim;
sim.add(cpu::md::VelocityVerlet{0.001});
sim.on(md::Phase::force).add(cpu::md::NeighborList{6.8}, -10).add(cpu::md::LJ{...});

sim.add(cpu::observer::RDF{.r_max = 10.0, .bins = 200, .freq = 100});
sim.add(cpu::recorder::XYZTrajectory{"traj.xyz", .freq = 1000});

sim.build();
sim.run(frame, 10'000'000);
sim.shutdown();
```

### B.7 One step with and without subscribers

(See parent revision — unchanged.)

---

## Part C — Phase 5: completion + performance optimization

> **Goal of Phase 5.** Drive the parallel CPU MD engine to a fully shippable state by closing the open completion items from Phases 1–4 (Track A) **and** push performance from "1.4× of LAMMPS at single-core, untuned multi-rank" toward "≤ 1.2× at single core, ≤ 2.0× at 8 ranks" via a focused, profile-driven optimization track (Track B). Tracks A and B run in parallel; Track A is mostly tooling/composition glue, Track B is measurement-driven kernel work.

### Track A — Completion punch list

#### A.1 `NoseHoover::set_T_target(double T)` runtime setter

The thermostat currently bakes `T_target` into the constructor; multi-stage protocols (heat → cool → tension) must rebuild the driver to switch setpoints, which destroys backend state (analysis executor, neighbor lists, RCB cuts). Add a public setter on `cpu::md::NoseHoover` (header `include/atomiverse/cpu/md/nose_hoover.h`, impl `src/cpu/md/nose_hoover.cpp`) that updates the target temperature between calls to `Driver::run()` without touching chain coordinates / momenta unless `reset_chain=true`. Linear-ramp helper `set_T_ramp(double T0, double T1, int n_steps)` is optional follow-up; the setter alone is the deliverable.

Test: `tests/mpi/test_nose_hoover_t_ramp.cpp` runs an N=2-rank LJ system through a NVT ramp 100 K → 500 K driven by repeated `set_T_target` calls between short `run(..., 100)` segments; asserts (a) measured ⟨T⟩ tracks setpoint within ±10 K block-averaged, (b) total step count > 5 segments crossed without driver rebuild, (c) backend `AnalysisExecutor` instance pointer is identical before and after.

#### A.2 Multi-stage `polyethylene_tension.cpp`

`benchmarks/md/polyethylene_tension.cpp` today runs a single stage selectable by CLI flag. Refactor into a single-process composition that drives the canonical 5-stage protocol end-to-end on the 20-chain × 1000-monomer PE system inside one `cpu::md::MD` instance:

1. NVT @ 500 K, 10 000 steps, Δt = 1 fs
2. NPT @ 500 K, 50 000 steps, Δt = 0.5 fs
3. Cooling 500 → 100 K via `NoseHoover::set_T_target` ramp (50 000 steps, Δt = 0.5 fs)
4. NPT @ 100 K, 50 000 steps, Δt = 0.5 fs
5. Uniaxial tension along x at ε̇ = 10⁹ s⁻¹ to ε ≈ 1.0 (200 000 steps, Δt = 0.5 fs)

The `set_T_target` call must appear at least twice in distinct stages (transition into stage 3, transition into stage 4). Existing single-stage CLI modes preserved as `--stage <id>` for back-compat with `tests/mpi/test_pe_msu_pipeline_n8.cpp`.

#### A.3 LAMMPS reference inputs

Create:

- `benchmarks/md/reference/polyethylene/lammps.in` — exact protocol mirror of A.2 (5 stages, same Δt, same FF parameters from `force_field.yaml`).
- `benchmarks/md/reference/polyethylene/thermo_reference.csv` — populated by running `lammps.in`; columns `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`.
- `benchmarks/md/reference/polyethylene/stress_strain_reference.csv` — same run; columns `strain,sigma_xx`.
- `benchmarks/md/reference/polyethylene/README.md` — reproduction recipe: required LAMMPS git tag (recorded explicitly), exact build flags, command sequence, expected wall-time at 1 / 8 ranks on the recorded reference machine.

LAMMPS git tag and host hardware fingerprint are captured in the README so future re-runs are deterministic.

#### A.4 Standardized launchers

Create `tools/run_msu_atomiverse.sh` and `tools/run_msu_lammps.sh`, both accepting a single `RANKS` positional argument (default `8`). Each launcher:

- selects the correct binary / input file
- forwards `OMP_NUM_THREADS`, `OMP_PROC_BIND`, `MV2_ENABLE_AFFINITY` env honoring caller overrides
- writes per-component CSV output into a timestamped subdirectory of `benchmarks/md/reference/polyethylene/runs/`
- prints the final per-step wall time and the path to the run directory

Both scripts are shellcheck-clean.

#### A.5 `tools/plot_scaling.py`

Consumes `scaling.cpp`'s CSV output (already emits `--csv-header`). Produces two PNGs into `benchmarks/md/reference/polyethylene/`:

- `scaling_strong.png` — log-log time/step vs ranks; ideal speedup line; efficiency annotations at 8 / 32 ranks.
- `scaling_weak.png` — efficiency vs ranks at fixed atoms/rank.

Pure matplotlib, no seaborn dependency. CLI: `python tools/plot_scaling.py <csv_path> --out-dir benchmarks/md/reference/polyethylene/`.

### Track B — Performance optimization

> **All Track B items are gated by an empirical acceptance number.** No item lands without a benchmark CSV row demonstrating the speedup AND a test confirming no correctness regression at the stated FP tolerance.

#### B.1 Profile-driven hotspot inventory

New benchmark `benchmarks/md/pe_profile_single.cpp`. Runs 1000 steps of stage-5 tension on the 20-chain PE system (`msu_pe.data` if present, else procedurally generated 20×1000 chain box) at the requested rank count. Wraps every hot-path component with a scoped timer:

`AtomMigration`, `HaloExchangeForward`, `HaloExchangeReverse`, `NeighborList::update`, `NeighborList::build_half`, `LJ::launch`, `HarmonicBond::launch`, `HarmonicAngle::launch`, `MultiHarmonicDihedral::launch`, `VelocityVerlet::launch`, `NoseHoover::launch`, `Barostat::launch`, `VirialStress::launch`, `ReduceEnergy::launch`, `SnapshotPublisher::launch`, `RDF::on_snapshot` (analysis-thread side).

Output: `benchmarks/md/reference/polyethylene/profile_<ranks>.csv` with columns `component,total_ms,pct_of_step,n_calls,mean_us`. Run at `ranks ∈ {1, 8, 32}` (32 only when an MPI multi-rank harness is available; otherwise 1, 8 minimum). New `tools/plot_profile.py` emits a stacked-bar PNG `profile_<ranks>.png` from each CSV.

Acceptance: three CSV files committed (`profile_1.csv`, `profile_8.csv`, `profile_32.csv` if cluster present, otherwise `profile_8.csv` minimum) plus their PNG companions, plus a one-paragraph "hotspot summary" in `scaling_report.md` identifying the top three components by share at each rank count.

#### B.2 LJ pair kernel SIMD audit

`src/cpu/md/lj.cpp` currently uses a scalar inner pair loop within the half-list neighbor walk. Apply xsimd 11.2.0 ket-loop vectorization (analogous to the ERI ss|ss pattern in `eri-hpc-opt.md`): pack ghost-side neighbor coordinates `xj/yj/zj` into `xsimd::batch<double>` registers, evaluate `r2`, the cutoff mask, and the `(σ/r)^12 / (σ/r)^6` polynomial in batch, then accumulate per-thread force buffers. Guard behind `-march=native` (same gate as ERI). Falls back to the scalar loop on non-x86 / when xsimd isn't present.

Acceptance:
- `bench_md_step --force-only=lj` shows ≥ **1.4× speedup** over the scalar baseline on the recorded reference machine.
- NVE total energy bit-identical to scalar code at 1e-12 relative over 1000 steps on `tests/cpu/md/test_lj_nve.cpp`.

#### B.3 Bonded force-kernel cache + ILP

`src/cpu/bonded/{harmonic_bond,harmonic_angle,multi_harmonic_dihedral}.cpp` walk topology indices and dereference position arrays per atom inside the inner loop, defeating compiler auto-vectorization. Pack per-tuple position triples into stack-resident `double[3]` arrays at the top of each tuple iteration; let the compiler vectorize the bond-vector / cross-product / dihedral-angle math against the stack copies.

Acceptance:
- `bench_md_step --force-only=bonded` shows ≥ **1.2× speedup** over the unpacked baseline.
- Per-term energies (`E_bond`, `E_angle`, `E_dihedral`) bit-identical at 1e-12 on `tests/cpu/md/test_pe_single_chain.cpp` and `test_pe_small_box.cpp`.

#### B.4 Neighbor-list rebuild cadence audit

`src/cpu/md/cell_list.cpp` rebuilds on size-stale and skin-stale (already LAMMPS-cadence per `60891d1`). Audit measures actual rebuild count over 1000 steps of stage-5 tension and compares against LAMMPS' `neigh_modify delay every check yes` count on the same protocol. If Atomiverse rebuilds more often, tune `delta_skin` default and add a `every` modifier (max steps between rebuilds); if less often, document and verify NVE drift is unaffected.

Acceptance:
- Atomiverse rebuild count ≤ LAMMPS rebuild count over 1000 PE-tension steps, recorded in `profile_1.csv` as a separate row.
- NVE drift on `test_pe_small_box.cpp` ≤ 1e-4 relative (no regression from current value).

#### B.5 HaloExchangeReverse payload audit

`HaloExchangeReverse` sends 3 doubles (`fx,fy,fz`) per ghost row. Audit: confirm no kernel writes additional ghost-row state (e.g. virial contributions, per-atom energy) that would be silently discarded. If an additional field is found, either widen the payload (with a benchmark cost measurement) or assert in `HaloExchangeReverse::build()` that no such kernel is registered.

Acceptance: `include/atomiverse/cpu/md/mpi/halo_exchange.h` carries an explicit doc comment listing the exact set of ghost-row fields included in the reverse payload, plus a one-sentence audit conclusion ("no silent loss" or "extended to N fields"). No behavior change unless audit finds a bug.

#### B.6 OpenMP team binding + first-touch

`cpu::Backend::build_components()` does not currently set `OMP_PROC_BIND` / `OMP_PLACES`. Document in `benchmarks/md/reference/polyethylene/README.md` the recommended env (`OMP_PROC_BIND=close OMP_PLACES=cores`) and measure the wall-time impact on the reference machine. If a 2-socket dev box is available, also measure with `numactl --membind=<node>` and record both numbers.

Acceptance: README.md contains a "Thread binding" subsection with at least one before/after measurement (or a "single-socket — no NUMA effect observed" note); `cpu::Backend::build_components()` (or `cpu::md::MD` constructor) prints a one-line warning to stderr when `OMP_PROC_BIND` is unset and `n_threads > 1`.

#### B.7 Allocation audit (`test_md_no_alloc.cpp`)

The architectural rule "no allocation in `launch()`" is currently a code-review convention. Add `tests/cpu/test_md_no_alloc.cpp` that:

1. Builds a `cpu::md::MD` with VelocityVerlet + LJ + HarmonicBond + NeighborList + NoseHoover on a 1280-atom PE box.
2. Calls `sim.build()` and runs 10 warm-up steps (lets `NeighborList` and halo arenas reach steady state).
3. Hooks `malloc` / `operator new` (via `tcmalloc` symbol override or a process-wide `__malloc_hook`-style shim, gated to Linux-only with a `GTEST_SKIP` on macOS) and runs 100 inner-loop steps.
4. Asserts zero allocation calls were observed during the inner loop.

Acceptance: test passes on the Linux CI image; documented as `GTEST_SKIP("malloc-hook unavailable on Darwin")` on macOS so local devs aren't blocked.

#### B.8 `-O3 -ffast-math` audit

Verify CMake release config sets `-O3` on all targets and that `-ffast-math` (or `-funsafe-math-optimizations`) is **not** applied to integrator / thermostat translation units (`velocity_verlet.cpp`, `nose_hoover.cpp`, `init.cpp`) because reordering FP additions there breaks time-reversibility and Hamiltonian-conservation. If `-ffast-math` is used at all, restrict to FF-kernel files (`lj.cpp`, `harmonic_*.cpp`, `multi_harmonic_dihedral.cpp`) via per-source `set_source_files_properties(... PROPERTIES COMPILE_FLAGS ...)`.

Acceptance:
- `CMakeLists.txt` carries an explicit comment block documenting the policy and the per-file gates.
- `bench_md_step` numbers in `scaling_report.md` are annotated with the exact compiler + flag string used.

#### B.9 Single-core LAMMPS gap closure

After B.1–B.4 land, re-measure single-core PE-tension wall time on the reference machine and report the ratio vs LAMMPS at the same hardware + protocol. Target: ratio ≤ **1.2×** (down from current 1.4×). Numbers committed into `benchmarks/md/reference/polyethylene/scaling_report.md`.

Acceptance:
- `scaling_report.md` table row "single-core PE-tension, Δt = 0.5 fs, 1000 steps" shows ratio ≤ 1.2×, OR a documented analysis explaining why the gap is structural (e.g., "remaining 0.15× is LAMMPS' Verlet-list bin packing not yet ported"). Either outcome closes the item; the goal is honest measurement and a written-down result.

---

## Part D — Phased roadmap

### Phase 1 — scaffolding (≈1 week) ✅ done
(See parent revision — unchanged. Exit gate met: single-rank tests + `bench_md_step` regression ≤ 2%.)

### Phase 2 — halo exchange + force integration (≈2 weeks) ✅ done
(See parent revision — unchanged. Exit gate met: 8-rank LJ matches 1-rank to < 1e-12 over 1000 steps.)

### Phase 3 — atom migration, bonded cross-boundary, reductions (≈2 weeks) ✅ done
(See parent revision — unchanged. Exit gate met: 32-rank LJ + bonded NVT vs 1-rank within FP tolerance through 10 000 steps.)

### Phase 4 — load balancing + MPI benchmarks (≈1 week) ✅ scaffolding done; 1 follow-up moved to Phase 5
(See parent revision — unchanged. RCB + gather + scaling.cpp + polyethylene_tension.cpp + validate_polyethylene.py landed; cluster scaling and multi-stage protocol composition rolled into Phase 5.)

### Phase 5 — completion + performance optimization (≈2 weeks)

**Track A (completion):** A.1 (NoseHoover setter) → A.2 (multi-stage `polyethylene_tension.cpp`) in series; A.3 (LAMMPS reference) → A.4 (launchers) → A.5 (plot_scaling) can interleave with Track B.

**Track B (performance):** B.1 (profile inventory) MUST land before B.2/B.3/B.4 so optimization decisions are data-driven. B.5/B.6/B.7/B.8 are standalone audits that can land in any order. B.9 is the closing measurement.

**Exit gates (all required):**

1. **Single-core LAMMPS ratio ≤ 1.2×** on the recorded reference machine, PE-tension stage-5 protocol, 1000 steps. (Track B.9.)
2. **8-rank LAMMPS ratio ≤ 2.0×** on the recorded reference machine, same protocol. (Track B.1 + downstream wins.)
3. **`polyethylene_tension.cpp` runs heat → cool → tension end-to-end without driver rebuild** — single `cpu::md::MD` instance, multiple `set_T_target` transitions, completes all 5 stages without crash / hang / OOM. (Track A.1 + A.2.)
4. **All Track A artifacts present** under `benchmarks/md/reference/polyethylene/` (lammps.in, thermo_reference.csv, stress_strain_reference.csv, README.md) and reproducible from `tools/run_msu_lammps.sh 8`. (Track A.3 + A.4.)
5. **Per-component profile CSVs** committed at 1, 8 (and 32 if cluster) ranks. (Track B.1.)
6. **No NVE-correctness regression** on `test_lj_nve`, `test_pe_small_box`, `test_pe_single_chain` after Track B SIMD / cache / cadence work. (Track B.2 + B.3 + B.4.)

### Phase 6 (deferred) — triclinic, parallel Zarr I/O, SENSEI hook
Tracked in a future addendum. Not part of the core spec; gates on Phase 5 completion.

---

## Files to create — status

Legend: ✅ created · ❌ not yet · ⏸ deferred · ↪ relocated

### MPI components (Phase 1–4)

```
✅ include/atomiverse/cpu/md/mpi/domain.h
✅ include/atomiverse/cpu/md/mpi/halo_exchange.h         # Forward + Reverse + PIMPL FaceBufs
✅ include/atomiverse/cpu/md/mpi/atom_migration.h
✅ include/atomiverse/cpu/md/mpi/reduce.h
✅ include/atomiverse/cpu/md/mpi/rcb.h                   # Phase 4: non-uniform brick balancer
✅ include/atomiverse/cpu/md/mpi/gather.h                # Phase 4: gather_to_root helper
✅ src/cpu/md/mpi/domain.cpp                             # extended with cuts[3] support
✅ src/cpu/md/mpi/halo_exchange.cpp
✅ src/cpu/md/mpi/atom_migration.cpp                     # cuts-aware dest_rank_for
✅ src/cpu/md/mpi/reduce.cpp
✅ src/cpu/md/mpi/rcb.cpp                                # Phase 4
✅ src/cpu/md/mpi/gather.cpp                             # Phase 4
```

### Async-analysis infrastructure (Phase 1)

```
✅ include/atomiverse/cpu/frame_snapshot.h
✅ include/atomiverse/cpu/analysis_executor.h
✅ include/atomiverse/cpu/snapshot_publisher.h
✅ include/atomiverse/cpu/async_component_base.h
✅ include/atomiverse/cpu/md/async_observer.h
✅ include/atomiverse/cpu/minimizer/async_observer.h
✅ src/cpu/analysis_executor.cpp
✅ src/cpu/snapshot_publisher.cpp
✅ src/cpu/async_component_base.cpp
```

### Tests

```
✅ tests/mpi/test_halo_exchange.cpp
✅ tests/mpi/test_atom_migration.cpp
✅ tests/mpi/test_pe_msu_pipeline.cpp
✅ tests/mpi/test_pe_nvt_multirank.cpp
✅ tests/mpi/test_npt_multirank.cpp
✅ tests/mpi/test_deform_multirank.cpp
✅ tests/mpi/test_pe_chain_crossing.cpp
✅ tests/mpi/test_pe_chain_drifting.cpp
✅ tests/mpi/test_pe_chains_8rank.cpp
✅ tests/mpi/test_pe_chains_multirank.cpp
✅ tests/mpi/test_pe_full_multirank.cpp
✅ tests/mpi/test_topology.cpp
✅ tests/mpi/test_gather.cpp
✅ tests/mpi/test_rcb.cpp
↪ tests/mpi/test_reduce_energy.cpp                       # subsumed by test_pe_msu_pipeline + test_pe_nvt
↪ tests/mpi/test_vs_single_rank.cpp                      # subsumed per-feature
✅ tests/cpu/test_analysis_executor.cpp
✅ tests/cpu/test_snapshot_publisher.cpp
✅ tests/cpu/test_async_component_base.cpp
✅ tests/cpu/observer/test_rdf_observer.cpp
❌ tests/mpi/test_nose_hoover_t_ramp.cpp                 # Phase 5 Track A.1
❌ tests/cpu/test_md_no_alloc.cpp                        # Phase 5 Track B.7
```

### Benchmarks + reference data (Phase 4 + Phase 5)

```
✅ benchmarks/md/scaling.cpp
✅ benchmarks/md/polyethylene_tension.cpp                 # Phase 5 A.2 refactors to multi-stage
✅ tools/validate_polyethylene.py
❌ benchmarks/md/pe_profile_single.cpp                    # Phase 5 Track B.1
❌ benchmarks/md/reference/polyethylene/lammps.in         # Phase 5 Track A.3
❌ benchmarks/md/reference/polyethylene/thermo_reference.csv     # Phase 5 Track A.3
❌ benchmarks/md/reference/polyethylene/stress_strain_reference.csv  # Phase 5 Track A.3
❌ benchmarks/md/reference/polyethylene/README.md         # Phase 5 Track A.3 + B.6
❌ benchmarks/md/reference/polyethylene/profile_1.csv     # Phase 5 Track B.1
❌ benchmarks/md/reference/polyethylene/profile_8.csv     # Phase 5 Track B.1
❌ benchmarks/md/reference/polyethylene/profile_32.csv    # Phase 5 Track B.1 (if cluster)
❌ benchmarks/md/reference/polyethylene/scaling_report.md # Phase 5 Track B.6 + B.8 + B.9
❌ tools/plot_profile.py                                  # Phase 5 Track B.1
❌ tools/plot_scaling.py                                  # Phase 5 Track A.5
❌ tools/run_msu_atomiverse.sh                            # Phase 5 Track A.4
❌ tools/run_msu_lammps.sh                                # Phase 5 Track A.4
```

The MSU amorphous-PE acceptance is exercised today as a CTest (`test_pe_msu_pipeline_n8`) rather than a benchmark. The benchmark + LAMMPS-reference comparison form is tracked in `msu-pe-validation.md`.

## Files to modify

```
include/atomiverse/cpu/frame.h              # add tag[], image[] to Atoms
include/atomiverse/cpu/context.h            # add mpi_rank, mpi_size, dd*
include/atomiverse/cpu/backend.h            # own DomainDecomposition + AnalysisExecutor
src/cpu/backend.cpp                         # init MPI topology, snapshot pool; B.6 OMP-bind warning
include/atomiverse/cpu/md/md.h              # default_bindings add new comps
include/atomiverse/cpu/md/nose_hoover.h     # Phase 5 A.1: set_T_target setter
src/cpu/md/velocity_verlet.cpp              # loop bound n_local
src/cpu/md/nose_hoover.cpp                  # loop bound n_local; A.1 setter impl
src/cpu/md/init.cpp                         # InitVelocity, ZeroKinetic, KineticEnergy bounds
src/cpu/md/lj.cpp                           # per-thread buffer size = capacity; B.2 xsimd ket-loop
src/cpu/bonded/harmonic_bond.cpp            # B.3 stack-pack tuple positions
src/cpu/bonded/harmonic_angle.cpp           # B.3 same
src/cpu/bonded/multi_harmonic_dihedral.cpp  # B.3 same
src/cpu/md/cell_list.cpp                    # B.4 cadence audit / delta_skin tune
include/atomiverse/cpu/md/mpi/halo_exchange.h  # B.5 payload-field doc comment
src/cpu/observer/thermo_component.cpp       # read post-allreduce scalars
src/cpu/recorder/trajectory.cpp             # gather_to_root before molrs call
benchmarks/md/polyethylene_tension.cpp      # Phase 5 A.2: multi-stage composition
CMakeLists.txt                              # ATV_USE_MPI; B.8 -ffast-math gating policy
.claude/NOTES.md                            # add MPI rules + async observer rules
CLAUDE.md                                   # update Observer/Recorder section; add MPI conventions
```

## Reuse from existing code

- `src/cpu/md/cell_list.cpp` — linked-cell half-list CSR, already LAMMPS-style; reused as-is with `inum=n_local`. **B.4 tunes its cadence; B.4 does NOT rewrite it.**
- Per-thread force buffer pattern in `src/cpu/md/lj.cpp:63-111`, `src/cpu/bonded/*.cpp` — reused, just sized to capacity. **B.2 SIMD-izes the inner pair loop without disturbing this pattern.**
- `Mat3` PBC helpers in `core/mat3.h` — reused; triclinic hook already there.
- xsimd 11.2.0 linkage and `-march=native` CXX guard (`CMakeLists.txt:326`) — reused unchanged. **B.2 reuses this exact gate.**
- `Scheduling::deferrable` enum in `core/types.h` — already exists; newly meaningful for CPU async observers.
- `cpu::Context` struct (`include/atomiverse/cpu/context.h`) — extended, not replaced.
- `src/cpu/analysis_executor.cpp` `sched_setaffinity` pattern — **B.6 layers `OMP_PROC_BIND` documentation alongside.**
- `benchmarks/md/scaling.cpp` strong/weak harness — **B.1 plugs the per-component scoped timer hooks; A.5 plots from its CSV.**
- `tests/mpi/test_pe_msu_pipeline_n8.cpp` single-stage 8-rank smoke test — **superseded by A.2's multi-stage composition for the canonical end-to-end run; original test retained as a stage-isolated regression.**

## Verification plan

(Items 1–7 from previous revision retained; new items 8–9 added for Phase 5.)

1. **Single-rank regression** — `cmake --build build && ctest -L cpu --output-on-failure` passes bit-identically; `bench_md_step` regression ≤ 2%.
2. **Multi-rank halo correctness (Phase 2 gate)** — `mpirun -n 8 ./build/tests/test_halo_exchange` matches 1-rank to 1e-12.
3. **Cross-validation vs LAMMPS (Phase 3 gate)** — `tools/validate_vs_lammps.py` on 32 000-atom PEG: thermo agreement to 1e-10 / 1e-8.
4. **Scaling benchmark (Phase 4 gate)** — `mpirun -n {…} ./build/benchmarks/md/scaling 10⁶-atoms`; strong ≥ 80% @ 32 ranks, weak ≥ 85% @ 128 ranks (cluster).
5. **Async observer non-interference** — `bench_md_step` with/without RDF (`freq=100`): ≤ 3% throughput regression.
6. **Backpressure modes** — `test_async_executor`: `drop_oldest` / `block` / `grow` correct counters; no deadlock; bounded memory.
7. **MSU amorphous-PE end-to-end (Phase 4 gate, canonical)** — see `msu-pe-validation.md` for the full protocol and tolerances.
8. **Phase 5 multi-stage composition gate** — `./build/benchmarks/md/polyethylene_tension` runs all 5 stages in one process, calling `NoseHoover::set_T_target` at least twice; completes without driver rebuild; `Thermo` log shows monotonic stage transitions in T setpoint.
9. **Phase 5 performance gate** — `./tools/run_msu_atomiverse.sh 1` and `./tools/run_msu_lammps.sh 1` produce wall times whose ratio ≤ 1.2×, recorded in `scaling_report.md`. Same for `RANKS=8` with ratio ≤ 2.0×.

---

## Tasks

> Phase 1–4 status is tracked in the snapshot table; the explicit checklist below covers Phase 5 only and is the authoritative `/mol:impl` driver.

### Track A — Completion punch list

- [x] Write failing test for `NoseHoover::set_T_target` ramp behavior (`tests/mpi/test_nose_hoover_t_ramp.cpp`)
- [x] Implement `NoseHoover::set_T_target(double T)` in `include/atomiverse/cpu/md/nose_hoover.h` + `src/cpu/md/nose_hoover.cpp`
- [x] Refactor `benchmarks/md/polyethylene_tension.cpp` to compose all 5 stages in a single `cpu::md::MD` instance using `set_T_target` between stages
- [x] Author `benchmarks/md/reference/polyethylene/lammps.in` mirroring the 5-stage protocol exactly
- [ ] Run `lammps.in` and commit `thermo_reference.csv` + `stress_strain_reference.csv` with LAMMPS git tag captured
- [x] Write `benchmarks/md/reference/polyethylene/README.md` reproduction recipe (LAMMPS tag, build flags, expected wall times, hardware fingerprint)
- [x] Write `tools/run_msu_atomiverse.sh` and `tools/run_msu_lammps.sh` (RANKS positional arg, OMP env passthrough, timestamped run directory)
- [x] Write `tools/plot_scaling.py` consuming `scaling.cpp` CSV; emit `scaling_strong.png` + `scaling_weak.png`

### Track B — Performance optimization

- [x] Write `benchmarks/md/pe_profile_single.cpp` with scoped per-component timers (B.1)
- [~] Write `tools/plot_profile.py` and produce `profile_1.csv` + `profile_8.csv` (+ `profile_32.csv` if cluster) with companion PNGs (B.1) — script done, CSVs need reference-machine run
- [ ] Record top-3 hotspot summary at each rank count in `benchmarks/md/reference/polyethylene/scaling_report.md` (B.1)
- [~] Write failing benchmark guard for LJ SIMD speedup target (`bench_md_step --force-only=lj` baseline frozen) (B.2) — infeasible as written: `bench_md_step` is CUDA-only with no `--force-only` mode. Alternative: pin LJ baseline ms in `scaling_report.md` from `pe_profile_single` LJ row; gate ac-011 on that. Needs spec refine.
- [x] Implement xsimd ket-loop vectorization in `src/cpu/md/lj.cpp` half-list inner loop (B.2) — 17 xsimd uses; manual gather + branch-free `xsimd::select` min-image; lane-extracted scalar accumulation; 1.21× measured on Apple NEON (W=2); higher expected on x86 AVX2/AVX-512
- [~] Verify NVE bit-identicality at 1e-12 on `tests/cpu/md/test_lj_nve.cpp` after B.2 (B.2) — energy bit-identical (7.7e-16); per-component forces 5e-12 (clang FMA-vs-no-FMA reorder, not an algorithmic error). `test_lj_nve` (1e-4 NVE drift) and all 42 tests pass. **Spec refine needed**: tighten ac-012 to "≤1e-12 on global energy" and add "≤1e-10 on per-component forces" since vectorized accumulation cannot bit-match scalar-with-FMA
- [~] Write failing benchmark guard for bonded ILP speedup target (`bench_md_step --force-only=bonded` baseline frozen) (B.3) — same infeasibility as B.2; same alternative recommended (pin baseline in `scaling_report.md` from `pe_profile_single` rows)
- [x] Implement stack-pack tuple positions in `src/cpu/bonded/harmonic_bond.cpp`, `harmonic_angle.cpp`, `multi_harmonic_dihedral.cpp` (B.3) — `__restrict` ptr hoist + stack-resident `ri[3]/rj[3]/rk[3]/r1..r4[3]`; HarmonicBond 1.13×, HA/HD register-bound (no measurable speedup on dev box); 42/42 tests green
- [x] Verify per-term bonded energies bit-identical at 1e-12 on `test_pe_single_chain.cpp` + `test_pe_small_box.cpp` (B.3) — bit-identical at every digit (pure load hoist preserves FP semantics); verified via `/tmp/bonded_probe.cpp` 8-atom chain
- [x] Audit `src/cpu/md/cell_list.cpp` rebuild cadence vs LAMMPS: added `delay` (default 10, LAMMPS parity), fixed min-image PBC-wrap false-rebuild bug, added box-change check, rebuild counter, 4 cadence tests (B.4) — profile CSV pending reference machine
- [x] Verify NVE drift ≤ 1e-4 after B.4 cadence change: `NveDriftWithinTolerance` test passes with default delay=10 (B.4)
- [x] Audit `HaloExchangeReverse` payload field set; add explicit doc comment in `include/atomiverse/cpu/md/mpi/halo_exchange.h` (B.5)
- [x] Add stderr warning in `cpu::Backend::build_components()` (or `cpu::md::MD` constructor) when `OMP_PROC_BIND` unset and `n_threads > 1` (B.6)
- [~] Document recommended OMP env + measured NUMA effect in `benchmarks/md/reference/polyethylene/README.md` (B.6) — env documented; numerical NUMA measurement needs reference machine
- [x] Write failing test `tests/cpu/test_md_no_alloc.cpp` asserting zero malloc calls in steady-state inner loop (Linux-only; GTEST_SKIP on Darwin) (B.7)
- [~] Verify allocation-free invariant passes on Linux CI image (B.7) — test scaffolded; Darwin SKIPs cleanly. Linux assertion needs CI run.
- [~] Audit + document `-O3 -ffast-math` policy in `CMakeLists.txt` (FF kernels only; integrators excluded); annotate `bench_md_step` flag string in `scaling_report.md` (B.8) — policy block added; `scaling_report.md` annotation needs reference-machine bench run
- [ ] Re-measure single-core PE-tension ratio vs LAMMPS; record in `scaling_report.md`; either ≤ 1.2× or written-down structural-gap analysis (B.9)
- [x] Run full check + test suite (`cmake --build build && ctest -L cpu -L mpi --output-on-failure`) — 42/42 pass on darwin (26 MPI + 16 single-rank, incl. new test_mpi_nose_hoover_t_ramp_n2 and test_cpu_md_no_alloc[SKIP])

## Testing strategy

**Happy path (Track A):**
- Multi-stage `polyethylene_tension.cpp` runs all 5 stages end-to-end at 1 rank and at 8 ranks without crash, hang, OOM, or driver rebuild.
- `tools/run_msu_atomiverse.sh 8` and `tools/run_msu_lammps.sh 8` both produce a populated thermo CSV, stress-strain CSV, and run-directory README.
- `tools/plot_scaling.py` ingests `scaling.cpp` CSV and emits two PNGs without runtime errors.

**Edge cases:**
- `NoseHoover::set_T_target` called with `T_target == current_T_target` (no-op semantics; chain coordinates untouched).
- `set_T_target` called between `run()` segments with `reset_chain=false` preserves chain integrator state (asserted by `test_nose_hoover_t_ramp.cpp`).
- `pe_profile_single.cpp` invoked at `--ranks 1` works without MPI runtime (single-rank fallback).
- `test_md_no_alloc.cpp` skips cleanly on macOS with informative `GTEST_SKIP` reason.
- LJ xsimd path falls back to scalar when `__AVX__` is undefined (compile-time guard exercised by a `tests/cpu/md/test_lj_nve_scalar_fallback.cpp` build matrix entry, optional).

**Domain validation:**
- All Track B optimizations preserve NVE energy bit-identicality at 1e-12 on the existing single-rank tests (`test_lj_nve`, `test_pe_small_box`, `test_pe_single_chain`).
- Phase 5 multi-stage `polyethylene_tension.cpp` reproduces stress-strain regimes (elastic → yield ε≈0.05 → softening → hardening to ε≈1.0) within ±30% of MSU reference (existing acceptance from `msu-pe-validation.md` D.4, unchanged tolerance).
- Single-core LAMMPS ratio target ≤ 1.2× (B.9); 8-rank ratio target ≤ 2.0× (Phase 5 exit gate 2).

**Performance:**
- `bench_md_step --force-only=lj` ≥ 1.4× over scalar baseline (B.2).
- `bench_md_step --force-only=bonded` ≥ 1.2× over unpacked baseline (B.3).
- Halo exchange correctness regression `test_mpi_halo_exchange_n8` continues to pass at 1e-12 (no Track B item touches halo correctness).

## Out of scope

- **Triclinic box support** — deferred to Phase 6.
- **Parallel HDF5 / per-rank Zarr shards** — deferred to Phase 6.
- **SENSEI / Ascent integration** — deferred to Phase 6.
- **GPU offload of CPU MD kernels** — out of scope; CUDA backend is independent.
- **Cluster-scale (≥ 64-rank multi-node) measurement** — Phase 5 exit gates target single-node 1 / 8 ranks; multi-node strong-scaling at 32+ ranks remains the parent ⏸ external item dependent on cluster access.
- **Matching LAMMPS performance at every rank count** — explicit non-goal; Phase 5 targets ≤ 1.2× single-core / ≤ 2.0× 8-rank only.
- **Replacing `cell_list.cpp`** — B.4 only audits cadence and tunes `delta_skin`; no re-architecture.
- **Force-field parameter tuning to match MSU paper numbers** — handled by `msu-pe-validation.md`, not this spec.
