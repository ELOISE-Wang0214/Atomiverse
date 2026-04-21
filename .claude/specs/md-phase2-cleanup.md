---
title: CPU MD Phase 2 — Architecture Cleanup
status: draft
created: 2026-04-21
---

## Summary

`av-architect` review of the Phase 2 hardening work (commit `a9f9230`) flagged
**3 BLOCKERs**, **7 HIGHs**, and several MEDIUM/LOW drift items in the CPU MD
module. This spec inventories each issue against its line anchor, proposes a
concrete fix, and groups the fixes into three execution phases so that
correctness-critical items (lifecycle violations, cross-instance state
corruption) land before any structural refactor or cosmetic polish.

**Scope:** `include/atomiverse/cpu/md/`, `src/cpu/md/` (including `mpi/`
subdirectory), plus one touch in `include/atomiverse/cpu/frame.h`.

**Non-scope:** any new MD feature, any CUDA-side change, Phase 4 work.

## Issues inventory

### BLOCKERs (hard rule violations)

| # | File | Line | Issue |
|---|------|------|-------|
| B1 | `src/cpu/md/mpi/halo_exchange.cpp` | 215–220 | `bufs_for(HaloExchangeForward*)` returns `static thread_local FaceBufs` — shared across every instance of the class in the process. Two drivers in the same process, or two `HaloExchangeForward` components in one driver, corrupt each other's per-face buffers. Same bug in the `HaloExchangeReverse` overload. |
| B2 | `include/atomiverse/cpu/md/neighbor_list.h` | 39–45 | `launch()` does `pairs_ = strategy_.build_half(...)` — the strategy returns a freshly-constructed `PairList` by value, dropping the previous xtensor storage and allocating new arrays of length `n_active` and `total_pairs` every step. Violates "no allocation in `launch()`". |
| B3 | `src/cpu/md/cell_list.cpp` | 336, 448 | `build_half()` / `build_full()` upstream of B2 — `result.neighbors = xt::zeros<int>({result.total_pairs})` reallocates the neighbors buffer unconditionally each call. |

### HIGHs (structural drift from CLAUDE.md rules)

| # | File | Line | Issue |
|---|------|------|-------|
| H1 | `include/atomiverse/cpu/frame.h` | 160 | `const void* pair_list = nullptr` — Frame carries component-communication state that is neither a fundamental physical quantity nor required by all drivers (minimizer and SCF never touch it). Violates "hot-path data structures carry only fundamental physical quantities". |
| H2 | `src/cpu/md/{nve_limit,velocity_verlet,nose_hoover_npt}.cpp` | see § Fix plan H2 | Three identical `wrap_pbc()` bodies in anonymous namespaces. Violates "no duplication across methods". |
| H3 | `src/cpu/md/{nve_limit,fix_langevin,nose_hoover_npt,velocity_verlet}.cpp` | see § Fix plan H3 | Four identical `owned_count(ctx)` helpers. |
| H4 | `src/cpu/md/nose_hoover_npt.cpp` | 32–41 | `n_global()` has both `if`/`else` branches returning `ctx.frame.atoms.size()`. Function has no caller; dead code. |
| H5 | `src/cpu/md/fix_langevin.cpp` | 66–82 | Box-Muller `GaussStream` constructed **inside** the per-atom loop every step. The `log/cos` pair dominates the cost; should be hoisted or replaced with the polar method. |
| H6 | `include/atomiverse/cpu/md/nose_hoover_npt.h` | 73–75 | Ctor takes 5 positional doubles `(T, tauT, P, tauP, dt)` — easy to flip `tauT`/`tauP` silently. Compare to `NoseHoover(T, tau, dt)` (3 args, unambiguous). |
| H7 | `src/cpu/md/mpi/halo_exchange.cpp` | 118 | `FaceBufs` is declared inside the anonymous namespace of the .cpp — paired with B1, the fix needs a PIMPL (or explicit forward-declared type in the header) to keep the buffer definition out of the public API. |

### MEDIUM / LOW (polish)

| # | File | Line | Issue |
|---|------|------|-------|
| M1 | `include/atomiverse/cpu/md/nose_hoover_npt.h` | 70 | `@param P_target` omits units. Both `Barostat` and `VirialStress` use `eV/Å³`; state it. |
| M2 | `include/atomiverse/cpu/md/fix_langevin.h` | 52–53 | `@param damp` doc narrative explains `γ = 1/damp` but the param line doesn't; copy the relationship into the `@param`. |
| M3 | `src/cpu/md/nose_hoover_npt.cpp` | 203 | `eps_dt_max = 0.1` hardcoded; promote to file-scope named constant with a doc comment or ctor arg. |
| M4 | `include/atomiverse/cpu/md/nose_hoover_npt.h` | 65 | Class-level "replaces VV + NoseHoover + Barostat" is a runtime contract with no enforcement. Add a `build()`-time assert that none of the three are also registered. |
| L1 | `include/atomiverse/cpu/md/mpi/halo_exchange.h` | 122–138 | `HaloExchangeReverse` missing doxygen `@brief` on `launch()` and `@see HaloExchangeForward`. |
| L2 | `src/cpu/md/fix_langevin.cpp` | 38 | Bare `3.14159265358979323846`; use `std::numbers::pi` (C++17 via `<numbers>` in C++20 — tree is C++17, so add `kPi` to `core/constants.h`). |

## Fix plan — phased

### Phase A: BLOCKER fixes (correctness critical)

Must land first. Each sub-step ships with a targeted test.

**A.1 — Halo buffer per-instance ownership (B1, H7)**

- Pick PIMPL: forward-declare `struct FaceBufs` in `include/atomiverse/cpu/md/mpi/halo_exchange.h`, hold `std::unique_ptr<FaceBufs> bufs_` as a member of both `HaloExchangeForward` and `HaloExchangeReverse`. The full `FaceBufs` definition stays in the .cpp.
- Delete `bufs_for(...)` helpers. `build()` allocates `bufs_ = std::make_unique<FaceBufs>()` and calls `prealloc_bufs(*bufs_, dd->capacity)`. `launch()` uses `*bufs_`.
- Test: `tests/mpi/test_halo_exchange_two_instances.cpp` — construct two `HaloExchangeForward` instances in the same driver (or two drivers), run each for a step, assert neither corrupts the other's send buffers. Requires MPI 2-rank to exercise the real comm path.

**A.2 — Persistent `PairList` arena (B2, B3)**

- Change `PairList::neighbors` from `xt::xtensor<int, 1>` to `std::vector<int>`: preserves capacity across `resize()`, gives us steady-state allocation-free rebuilds. `numneigh` and `firstneigh` stay xtensor (they resize to `n_active`, bounded by `dd->capacity`, so one reserve in `build()`).
- Add `CellList::refresh_half(const Frame&, int n_active, PairList& out)` that reuses `out.neighbors`, `out.numneigh`, `out.firstneigh` buffers. `out.neighbors.resize(total_pairs)` uses std::vector's geometric growth (amortized zero alloc once steady state is reached).
- `NeighborList::build(ctx)` preallocates `pairs_.numneigh.resize({capacity})`, `pairs_.firstneigh.resize({capacity})`, `pairs_.neighbors.reserve(capacity * typical_neighbors_per_atom)`. `launch(ctx)` calls `strategy_.refresh_half(frame, n_active, pairs_)`.
- Keep the existing `CellList::build_half(Frame, n_active) -> PairList` overload as a thin wrapper that delegates to `refresh_half` on a local `PairList` — tests and one-shot callers continue to compile unchanged.
- Test: `tests/cpu/md/test_neighbor_list_no_alloc.cpp` — record `malloc` count across 100 `launch()` calls with a stable `n_active`; assert no growth after step 5.
- Force kernels that iterate `neighbors[...]` (`lj.cpp`, etc.) already use `neighbors[off + jj]` — `std::vector<int>::operator[]` returns the same type as the old xtensor version. No kernel change.

### Phase B: HIGH structural drift

**B.1 — Frame depollution (H1)**

- Remove `const void* pair_list` from `cpu::Frame`.
- Pass the active `PairList*` via `cpu::LaunchContext` (new field `pair_list_view`, set by `NeighborList::launch` before returning, consumed by any component that binds after neighbor list in the force phase). Set semantics: `NeighborList` writes it; force components read it; outside the force phase it is a dangling view (documented).
- Alternative: a `cpu::PairListPublisher` singleton component (mirrors `SnapshotPublisher` if that exists). Lower churn but more plumbing — defer unless `LaunchContext` option causes phase-ordering headaches.
- Touch: `NeighborList::launch`, every force component that reads `frame.pair_list` (grep list below), tests.
- Test: existing `test_neighbor_list.cpp` covers the contract; add assert that `ctx.pair_list_view != nullptr` inside force components.

**B.2 — Dedup shared helpers (H2, H3)**

- New file: `include/atomiverse/cpu/md/pbc.h` (header-only):
  ```cpp
  inline void wrap_orthorhombic(cpu::Frame& frame, int n_owned) { ... }
  ```
- New file: `include/atomiverse/cpu/md/owned_range.h` (header-only):
  ```cpp
  inline int owned_count(const cpu::LaunchContext& ctx) { ... }
  ```
  Or promote to an inline method on `cpu::LaunchContext` itself — cleaner, touches one header.
- Sweep: replace 3 × `wrap_pbc` and 4 × `owned_count` with the shared helper. Delete the anonymous-namespace duplicates.

**B.3 — Dead code purge (H4)**

- Delete `n_global()` in `src/cpu/md/nose_hoover_npt.cpp:32-41`. No caller; both branches identical anyway.

**B.4 — Langevin Box-Muller hoist (H5)**

- Move `GaussStream` construction out of the per-atom loop. One stream per step, seeded with `(seed, step)`; per-atom uniqueness comes from advancing the stream by `3 * i + component`.
- Optional follow-up: switch to Marsaglia polar method (avoids trig). Defer to HPC pass — correctness first.
- Test: existing Langevin T-target test should still pass within statistical tolerance.

**B.5 — `NoseHooverNPT` ctor ergonomics (H6)**

- Replace the 5-positional-double ctor with an aggregate-initialized config struct:
  ```cpp
  struct NoseHooverNPT::Config {
      double T_target;   // [K]
      double tau_T;      // [fs]
      double P_target;   // [eV/Å³]
      double tau_P;      // [fs]
      double dt;         // [fs]
  };
  explicit NoseHooverNPT(Config cfg);
  ```
  Call sites become `NoseHooverNPT{{.T_target = 300.0, .tau_T = 100.0, ...}}` — designated initializers (C++20) or named-argument idiom.
- Project is C++17; use a builder pattern or require-all-fields-named via struct with comments. Prefer the struct + comment approach to avoid a C++ version bump. Either way it costs one commit in the test suite.

### Phase C: MEDIUM / LOW polish

**C.1** — doxygen unit fix: `P_target [eV/Å³]` (M1).
**C.2** — doxygen clarify `damp = 1/γ` in `@param` line (M2).
**C.3** — `eps_dt_max` → file-scope `kEpsDtMaxFraction` with doc (M3).
**C.4** — `NoseHooverNPT::build()` asserts absence of VV + NoseHoover + Barostat (M4). Requires `cpu::md::MD` exposing a "component registered" query — add if missing (small facade touch).
**C.5** — `HaloExchangeReverse` doxygen polish (L1).
**C.6** — bare π literal → `ATV::kPi` in `core/constants.h` if not already (L2).

## Verification

Per phase:

```bash
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
mpirun -np 2 ./build/tests/test_pe_chain_crossing
mpirun -np 2 ./build/tests/test_halo_exchange_two_instances   # new (A.1)
./build/tests/test_neighbor_list_no_alloc                     # new (A.2)
```

Acceptance criteria:

1. All existing tests still green.
2. New A.1 test: two `HaloExchangeForward` instances do not clobber each other (can be verified by checksum of `send_pos[f]` buffers after a known workload on each instance).
3. New A.2 test: `malloc` / new-op counter stable (±0) across steps 10..100 after warmup.
4. Post-Phase B: `grep -r "wrap_pbc\|owned_count" src/cpu/md/` returns only the canonical definitions (one each).
5. Post-Phase B: `grep "const void\* pair_list" include/atomiverse/cpu/frame.h` returns nothing.
6. Post-Phase B: `cpu::Frame` size (sizeof) shrinks by the size of a pointer (not enforced, but a nice sanity check in a debug print).

## Design decisions

1. **PIMPL over header-exposed `FaceBufs` (A.1):** keeps the public halo header lean and MPI-type-free. `std::unique_ptr<FaceBufs>` in the member list pays a 24B per-instance cost, negligible vs the per-face buffers themselves.
2. **`std::vector<int>` for `PairList::neighbors` (A.2):** xtensor does not guarantee capacity preservation across `resize()`. `std::vector` does. Force kernels use `operator[]` on both — no API break.
3. **`LaunchContext` over publisher component (B.1):** existing `LaunchContext` already flows into every component's `launch()`; adding a field is lower churn than introducing a new component and wiring its phase ordering. If a future driver (minimizer-with-neighbor-list) has different phase semantics, promote to publisher.
4. **Config struct over builder (B.5):** builder is more verbose and less discoverable than `NoseHooverNPT{{T, tauT, P, tauP, dt}}`. C++17 allows `-Wmissing-field-initializers` to catch incomplete aggregates, which is the safety we want.

## Execution order and estimated scope

| Phase | Blocks | Est. diff | Risk |
|-------|--------|-----------|------|
| A.1 — halo buffer | B1, H7 | ~80 LOC + 1 test (~100 LOC) | low — purely mechanical, covered by new test |
| A.2 — PairList arena | B2, B3 | ~150 LOC (CellList + NeighborList) + 1 test | medium — touches hot path; must verify no perf regression |
| B.1 — Frame depollution | H1 | ~30 LOC spread across ~8 force components | low — grep-and-replace |
| B.2 — dedup | H2, H3 | ~50 LOC deleted, 2 headers added (~20 LOC) | low — pure refactor |
| B.3 — dead code | H4 | -10 LOC | trivial |
| B.4 — Langevin hoist | H5 | ~20 LOC | low — tested by existing T-target check |
| B.5 — NPT ctor | H6 | ~15 LOC + test updates | low |
| C.1–C.6 — polish | M1–M4, L1, L2 | ~40 LOC | trivial |

**Total:** ~500 LOC churn, 2 new tests, zero behavior change on green paths.

## Out of scope

- Any Phase 4 work (RCBBalancer, scaling harness).
- Halo exchange for non-orthorhombic boxes.
- Bonded-tuple exchange (handled separately in halo + ghost mechanism already).
- CUDA-side equivalents of these fixes.
- Switching `PairList::numneigh` / `firstneigh` away from xtensor — their resize pattern is fixed-size per call, and they are consumed via xtensor expressions elsewhere.
