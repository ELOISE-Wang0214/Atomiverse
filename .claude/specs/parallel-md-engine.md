---
title: Parallel CPU MD engine — MPI + OpenMP + xsimd with hot-path/analysis separation
status: in-progress
created: 2026-04-16
last-updated: 2026-04-29 (Phase 4 scaffolding: RCB + gather + scaling/PE benchmarks landed)
---

# Parallel CPU MD engine

LAMMPS-class multi-node CPU MD: MPI spatial domain decomposition, OpenMP within each rank, xsimd in force kernels, and a dedicated per-rank analysis jthread pool for observers and recorders so the hot path never waits on diagnostics. Target scale: 1–100 ranks, 10⁴–10⁷ atoms. No Kokkos, no new driver types, no new phase enums.

## Status snapshot (2026-04-29)

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
| 4 | Multi-stage protocol composition (heat → cool → tension via persistent backend) | ⏸ followup | needs `NoseHoover::set_T_target` or persistent driver — present benchmark runs single stage |
| 5 | Triclinic / parallel Zarr / SENSEI hook | ⏸ deferred | future addendum |

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
- **No allocation in `launch()`.** Halo buffers, migration scratch, and snapshot ring are all allocated in `build()` at `capacity = ceil(n_atoms * (1 + slack))`, `slack = 0.20` default. `AtomMigration::launch()` shuffles rows inside `[0, capacity)` and signals an out-of-capacity generation bump when `n_local + incoming > capacity`; it never reallocates. Capacity growth runs through a `Driver::rebuild_frame()` hook that invokes component `build()` again — allocation stays in the build pathway. No exceptions.
- **Observers are components, not a separate API.** Same `sim.add(obs)` surface. Async vs sync is a component-internal property, not a user-facing distinction.

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

Phase 1–3: orthorhombic only (matches current `minimum_image` in `cell_list.cpp`). Triclinic support deferred to a Phase 5 addendum.

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
    // RAII thread wrapper: std::thread + std::atomic<bool> stop + std::condition_variable.
    // No std::jthread (C++17 only). shutdown() sets stop, notifies cv, joins each worker.
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

Binds `end_of_step @ order +10` — **after** `ZeroKinetic(-1000)`, `KineticEnergy(-100)`, and `ReduceEnergy(-50)`, **before** user-registered observers (which default to order 0+). This ordering guarantees that any scalar or per-atom field the snapshot captures (positions, velocities, forces, plus `frame.energy/kinetic/virial` if marked in `CaptureMask`) is the globally-reduced, end-of-step-consistent value. Running earlier would hand async observers stale-zero KE or per-rank energy totals.

Runs only when the set of registered async components, given each one's `frequency()`, indicates at least one subscriber at the current step. Parallel-copies requested fields from `cpu::Frame[0..n_local)` into a free ring slot using the hot path's OpenMP team (~5 µs per 10⁵ atoms). Zero subscribers → 5 ns early return.

### B.4 `AsyncComponentBase` (new base class) — driver-agnostic

Name is deliberately role-neutral — the class is pure async-transport plumbing (snapshot enqueue + scheduling trait), not "observer." Both observers (which compute) and recorders (which write I/O) inherit directly. The Observer-vs-Recorder distinction is preserved in subclasses by responsibility, not by an intermediate base class.

**Driver-agnostic by design.** The base class declares **no default bindings**. Each concrete driver's phase enum (`md::Phase`, `minimizer::Phase`, `scf::Phase`) has its own integer-indexed `end_of_step`, so hardcoding one in the base would silently misbind under any non-MD driver. Two usage paths:

1. **Cross-driver user observer** — inherit `cpu::AsyncComponentBase` directly; bind explicitly in each context:
   ```cpp
   md_sim.on(md::Phase::end_of_step).add(obs, 20);
   lbfgs_sim.on(minimizer::Phase::end_of_step).add(obs, 20);
   ```

2. **Driver-scoped user observer** — inherit a driver-specific thin wrapper:
   - `cpu::md::AsyncObserver`        (binds `md::Phase::end_of_step @ +20`)
   - `cpu::minimizer::AsyncObserver` (binds `minimizer::Phase::end_of_step @ +20`)
   - (add `cpu::scf::AsyncObserver` when an `scf::SCF` driver class lands)

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

Inside `launch()`: skip if `step == 0 || step % frequency() != 0`, read `exec_->last_published_slot()` (set by `SnapshotPublisher::publish()` earlier in the same `end_of_step`), verify the slot's step tag matches, call `exec_->submit(this, slot)`. ~30 ns.

`on_snapshot()` runs on the analysis worker thread and may call `MPI_Allreduce(exec_->analysis_comm(), ...)` freely (MPI builds).

**`SnapshotPublisher` is also driver-agnostic.** Each driver's constructor wires it explicitly:
```cpp
// cpu::md::MD
this->on(md::Phase::end_of_step).add(cpu::SnapshotPublisher{}, 10);
// cpu::minimizer::LBFGS
this->on(minimizer::Phase::end_of_step).add(cpu::SnapshotPublisher{}, 10);
```
No per-driver `SnapshotPublisher` subclass is needed — the shared class has empty default bindings, and drivers select phase through `on()`.

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

sim.build();        // starts jthread pool, allocates snapshot ring
sim.run(frame, 10'000'000);
sim.shutdown();     // drain analysis queue, finalize observers, join threads
```

### B.7 One step with and without subscribers

**Empty subscribers (production throughput run):**
```
main:  init_integrate → pre_force → force → post_force → final_integrate
     → end_of_step:
         ZeroKinetic(-1000) → KineticEnergy(-100) → ReduceEnergy(-50)
       → SnapshotPublisher(+10) (≤ 5 ns early return, zero subscribers)
```

**With async RDF (freq=100) at step 100:**
```
main: … → end_of_step:
    ZeroKinetic → KineticEnergy → ReduceEnergy(allreduce KE/PE/virial)
    SnapshotPublisher(+10): acquire slot, OMP-copy x,y,z,box + scalars (~5 µs at 10⁵ atoms)
    RDF::launch(0): submit to executor                                  (~30 ns)
analysis worker thread (concurrent with step 101's init_integrate):
    dequeue → RDF::on_snapshot(snap) → MPI_Allreduce(analysis_comm)
```

Hot path blocking time at `end_of_step` is bounded by the snapshot copy, *not* by observer compute. Observer can lag by multiple MD steps without stalling.

---

## Part C — Phased roadmap

### Phase 1 — scaffolding (≈1 week)
- Add `tag[]`, `image[]` to `cpu::Frame::Atoms`. Default `tag[i]=i`, `image=0`.
- Create `DomainDecomposition` owned by `cpu::Backend`; single-rank construction only (`dims={1,1,1}`, `n_local=frame.atoms.size()`, `n_ghost=0`, `capacity=n_local`, no oversize slack). Rank ≥ 2 exercised starting in Phase 2.
- Phase 1 single-rank invariant: `n_local == n_ghost == 0's complement == atoms.size() == capacity`. Integrator / force loops unchanged in Phase 1; the `n_local` loop-bound audit is deferred to Phase 2 where it becomes observable.
- `cpu::Context` gains `md::mpi::DomainDecomposition* dd = nullptr;` (opaque forward-decl in `context.h`, full type pulled in only by MD components).
- Add no-op stub components: `HaloExchangeForward`, `HaloExchangeReverse`, `AtomMigration`, `ReduceEnergy` (all `cpu::ComponentBase` subclasses). Wire into `cpu::md::MD` constructor alongside existing `ZeroForces`/`ZeroEnergy`/`ZeroKinetic`/`KineticEnergy`.
- Add `AnalysisExecutor`, `SnapshotPublisher`, `AsyncComponentBase`, `FrameSnapshot` at the **cpu root**. Wire `SnapshotPublisher` into `cpu::md::MD` default bindings at `end_of_step @ +10` (after `KineticEnergy @ -100` and `ReduceEnergy @ -50`).
- Port `cpu::observer::RDF` from free-function APIs to a `cpu::md::AsyncObserver` subclass (landed as `cpu::observer::RDF`; free functions `rdf_accumulate`/`rdf_write` preserved for CUDA callers).
- `cpu::observer::Density` deferred: SCF-scoped, requires `FrameSnapshot::custom_payload` to carry P from `SCFEvaluator`, plus a concrete SCF driver class for auto-wiring. Tracked under AIMD density-observation workstream (Phase 1b / Phase 5).
- CMake: `ATV_USE_MPI` option. `find_package(MPI)` already present (line 133). Keep `CMAKE_CXX_STANDARD 17`.
- **Exit gate:** all existing single-rank tests and benchmarks pass bit-identically; `bench_md_step` regression ≤ 2%.

### Phase 2 — halo exchange + force integration (≈2 weeks)
- Implement `HaloExchange` forward (positions) and reverse (forces) passes with `MPI_Isend/Irecv` over 6 face neighbors.
- Loop bound `n_local` in `VelocityVerlet`, `NoseHoover`, `InitVelocity`, `ZeroKinetic`, `KineticEnergy`.
- Per-thread force buffers in `LJ`, `HarmonicBond`, `HarmonicAngle`, `MultiHarmonicDihedral` grow to `n_local + n_ghost_capacity`.
- New test: `tests/mpi/test_halo_exchange.cpp` on 2×2×2 brick with 1000 LJ atoms; verify `sum(fx)` bit-equal to single-rank after allreduce.
- **Exit gate:** LJ liquid at 8 ranks matches 1-rank energy to < 1e-12 relative over 1000 steps.

### Phase 3 — atom migration, bonded cross-boundary, reductions (≈2 weeks)
- Implement `AtomMigration` with `MPI_Alltoallv`. Bump `neighbor_generation`. Rebuild `tag_to_local`.
- Implement `ReduceEnergy` (batched allreduce). Add `ReduceVirial` as a separate component (order −190).
- Bonded ownership: lowest-tag for bonds, center for angles, 2nd-atom for dihedrals. Topology migration alongside atoms.
- Cross-validate against LAMMPS on a 32-rank LJ + bonded PEG chain benchmark. Tolerance 1e-10 relative on energy.
- New CTest label `mpi`; `add_test(… ${MPIEXEC_EXECUTABLE} -n N …)`.
- **Exit gate:** full NVT run of 10⁵-atom LJ + harmonic bonds on 32 ranks matches 1-rank trajectory energy to within FP tolerance through step 10 000.

### Phase 4 — load balancing + MPI benchmarks (≈1 week)
- `RCBBalancer` component (in-house RCB, ~300 LOC; Zoltan optional).
- `benchmarks/md/scaling.cpp` — strong and weak scaling at 1/2/4/8/16/32/64/128 ranks.
- Per-rank analysis thread affinity hinting (via `sched_setaffinity`) — optional.
- **Exit gate:** 80% strong-scaling efficiency at 32 ranks on 10⁶-atom LJ; weak scaling efficiency > 85% at 128 ranks.
- **Canonical end-to-end gate (see Verification §7):** the full MSU amorphous-PE uniaxial tension protocol at 8 ranks completes with stress-strain regimes and hot-path / async-observer criteria listed there.

### Phase 5 (deferred) — triclinic, parallel Zarr I/O, SENSEI hook
Tracked in a future addendum. Not part of the core spec; gates on Phase 4 completion.

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
✅ tests/mpi/test_pe_msu_pipeline.cpp                    # MSU end-to-end smoke (n=8)
✅ tests/mpi/test_pe_nvt_multirank.cpp                   # NVT correctness vs single-rank
✅ tests/mpi/test_npt_multirank.cpp                      # NPT correctness vs single-rank
✅ tests/mpi/test_deform_multirank.cpp                   # uniaxial-deform correctness
✅ tests/mpi/test_pe_chain_crossing.cpp                  # chain crossing rank boundary
✅ tests/mpi/test_pe_chain_drifting.cpp                  # multi-step drift / migration
✅ tests/mpi/test_pe_chains_8rank.cpp                    # 8-rank PE bonded NVE
✅ tests/mpi/test_pe_chains_multirank.cpp                # 2-rank PE bonded NVE
✅ tests/mpi/test_pe_full_multirank.cpp                  # full pipeline 2-rank
✅ tests/mpi/test_topology.cpp                           # MPI Cartesian + partition
✅ tests/mpi/test_gather.cpp                             # Phase 4: gather_to_root single + multi rank
✅ tests/mpi/test_rcb.cpp                                # Phase 4: RCB rebalance + atom redistribution
↪ tests/mpi/test_reduce_energy.cpp                       # subsumed by test_pe_msu_pipeline + test_pe_nvt
↪ tests/mpi/test_vs_single_rank.cpp                      # subsumed per-feature
✅ tests/cpu/test_analysis_executor.cpp                  # at tests/cpu/, not tests/observer/
✅ tests/cpu/test_snapshot_publisher.cpp
✅ tests/cpu/test_async_component_base.cpp
✅ tests/cpu/observer/test_rdf_observer.cpp              # observer-side RDF async test
```

### Benchmarks + reference data (Phase 4)

```
✅ benchmarks/md/scaling.cpp                              # strong/weak scaling, --csv-header
✅ benchmarks/md/polyethylene_tension.cpp                 # NVT/NPT/tension modes, CSV recorder
✅ tools/validate_polyethylene.py                         # modulus/yield/hardening pass-fail
❌ benchmarks/md/reference/polyethylene/lammps.in         # LAMMPS reference run
❌ benchmarks/md/reference/polyethylene/thermo_reference.csv     # populated by ref run
❌ benchmarks/md/reference/polyethylene/stress_strain_reference.csv
```

The MSU amorphous-PE acceptance is exercised today as a CTest (`test_pe_msu_pipeline_n8`) rather than a benchmark. The benchmark + LAMMPS-reference comparison form is tracked in `msu-pe-validation.md`.

## Files to modify

```
include/atomiverse/cpu/frame.h              # add tag[], image[] to Atoms
include/atomiverse/cpu/context.h            # add mpi_rank, mpi_size, dd*
include/atomiverse/cpu/backend.h            # own DomainDecomposition + AnalysisExecutor
src/cpu/backend.cpp                         # init MPI topology, snapshot pool
include/atomiverse/cpu/md/md.h              # default_bindings add new comps
src/cpu/md/velocity_verlet.cpp              # loop bound n_local
src/cpu/md/nose_hoover.cpp                  # loop bound n_local
src/cpu/md/init.cpp                         # InitVelocity, ZeroKinetic, KineticEnergy bounds
src/cpu/md/lj.cpp                           # per-thread buffer size = capacity
src/cpu/bonded/harmonic_bond.cpp            # topology iteration uses local rows
src/cpu/bonded/harmonic_angle.cpp           # same
src/cpu/bonded/multi_harmonic_dihedral.cpp  # same
src/cpu/md/cell_list.cpp                    # inum = n_local; ghost indices in neighbors[]
src/cpu/observer/thermo_component.cpp       # read post-allreduce scalars
src/cpu/recorder/trajectory.cpp             # gather_to_root before molrs call
CMakeLists.txt                              # ATV_USE_MPI option, link MPI::MPI_CXX, ctest_mpi wiring
.claude/NOTES.md                            # add MPI rules + async observer rules
CLAUDE.md                                   # update Observer/Recorder section; add MPI conventions
```

## Reuse from existing code

- `src/cpu/md/cell_list.cpp` — linked-cell half-list CSR, already LAMMPS-style; reused as-is with `inum=n_local`.
- Per-thread force buffer pattern in `src/cpu/md/lj.cpp:63-111`, `src/cpu/bonded/*.cpp` — reused, just sized to capacity.
- `Mat3` PBC helpers in `core/mat3.h` — reused; triclinic hook already there.
- xsimd 11.2.0 linkage and `-march=native` CXX guard (`CMakeLists.txt:326`) — reused unchanged.
- `Scheduling::deferrable` enum in `core/types.h` — already exists; newly meaningful for CPU async observers.
- `cpu::Context` struct (`include/atomiverse/cpu/context.h`) — extended, not replaced.

## Verification plan

1. **Single-rank regression**
   `cmake --build build && ctest -L cpu --output-on-failure`
   Must pass bit-identically to pre-spec baseline. Benchmark `bench_md_step` must regress ≤ 2%.

2. **Multi-rank halo correctness (Phase 2 gate)**
   `mpirun -n 8 ./build/tests/test_halo_exchange`
   Compares `sum(fx)` after allreduce against single-rank force sum; tolerance 1e-12 relative.

3. **Cross-validation vs LAMMPS (Phase 3 gate)**
   Script at `tools/validate_vs_lammps.py` runs both engines on a 32 000-atom LJ + harmonic-bond PEG system, compares thermo output at step 100 / 1 000 / 10 000. Tolerance 1e-10 relative on energy, 1e-8 on temperature.

4. **Scaling benchmark (Phase 4 gate)**
   `mpirun -n {1,2,4,8,16,32,64,128} ./build/benchmarks/md/scaling 10^6-atoms`
   Reports time per MD step. Strong-scaling efficiency ≥ 80% at 32 ranks; weak-scaling ≥ 85% at 128 ranks.

5. **Async observer non-interference**
   `bench_md_step` with and without registered RDF (`freq=100`): hot-path throughput difference ≤ 3%.

6. **Backpressure modes**
   `test_async_executor` injects a slow observer and verifies `drop_oldest` / `block` / `grow` policies behave per spec (counters, no deadlock, bounded memory).

7. **Final end-to-end acceptance — MSU amorphous-PE uniaxial tension** (Phase 4 gate, canonical)
   Reference protocol: <https://icme.hpc.msstate.edu/mediawiki/index.php/Deformation_of_Amorphous_Polyethylene.html>
   - System: 20 chains × 1 000 united-atom CH₂ monomers = **20 000 atoms** on an FCC seed at 1.53 Å spacing, Dreiding-style force field (`HarmonicBond` + `HarmonicAngle` + `MultiHarmonicDihedral` + `LJ` for non-bonded).
   - Protocol (see `benchmarks/md/polyethylene_tension.cpp`):
     1. 10 000 steps NVT @ 500 K, `Δt = 1 fs` (`InitVelocity` + `NoseHoover`).
     2. 50 000 steps NPT @ 500 K, 0 bar, `Δt = 0.5 fs` (`NoseHoover` + `Barostat`).
     3. 50 000 steps cooling 500 K → 100 K (linear ramp via `NoseHoover` T-schedule).
     4. 50 000 steps NPT relax @ 100 K.
     5. 200 000 steps uniaxial tensile load along x at `ε̇ = 10⁹ s⁻¹` to ε ≈ 1.0, transverse faces zero-pressure (decoupled `Barostat`); `Deform` drives x-strain, `VirialStress` logs σ_xx(ε).
   - Observables: `Thermo` (T, P, E_bond, E_angle, E_dihedral, E_vdW, E_kinetic), `VirialStress` σ(t), `XYZTrajectory` every 10 000 steps.
   - Acceptance:
     * **Runs to completion at 8 ranks** on a single node without hang, deadlock, memory growth > 5%/10⁶ steps, or any MPI error.
     * **Stress-strain shape matches published qualitative regimes** (elastic → yield near ε ≈ 0.05, softening, hardening to ε ≈ 1.0). Peak yield stress within ±30% of the MSU reference at 100 K / `ε̇ = 10⁹ s⁻¹` (stringent quantitative match against LAMMPS is tracked separately; Atomiverse uses Dreiding-style, not identical parameter set).
     * **Energy-conservation sanity**: in the NVE portion (optional probe run without thermostat) total energy drift ≤ 1 × 10⁻⁴ eV/atom over 10 000 steps.
     * **Async observers non-blocking**: running the tension stage with `VirialStress` (sync, every step) + `RDF` (async, every 1 000 steps) + `XYZTrajectory` (async, every 10 000 steps) shows ≤ 3% hot-path throughput regression vs the same run with observers disabled.
     * **Analysis-thread MPI collective isolation**: RDF's `MPI_Allreduce` runs on the analysis sub-communicator without blocking step `k+1`'s `HaloExchange`.
   - Storage: reference LAMMPS input/output under `benchmarks/md/reference/polyethylene/`. Pass/fail script: `tools/validate_polyethylene.py` compares σ_xx(ε) curve and log quantities against reference CSV with tolerances above.

---

## Rationale vs alternatives (for future maintainers)

- **Not Kokkos.** Preserves C++17, existing xsimd ERI work, and every working force kernel. Trade-off: no CPU→GPU portability from one source; acceptable because the CUDA backend is already separate.
- **Not SENSEI/Ascent.** Target is 1–100 ranks — the jthread + triple-buffer pattern is sufficient and ~10× cheaper to own. SENSEI is a Phase 5 hook point if ever needed.
- **No new driver.** `cpu::md::MD` gains components via `default_bindings()`; the Driver template signature is untouched. Upholds the CLAUDE.md "Driver is the single canonical facade" invariant.
- **No new phase enum.** MPI exchange and async dispatch slot into existing `pre_force/post_force/end_of_step` via ordered bindings. Keeps the CUDA backend, AIMD (`SCFEvaluator`), and minimizer drivers unchanged.
