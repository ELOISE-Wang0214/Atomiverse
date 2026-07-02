# Cell-list OMP rebuild: row-owned finalize + candidate-loop hoisting

**Status:** draft (review-derived implementation plan, 2026-07-02).
**Origin:** multi-axis CPU review — perf axis 🔴 `cell_list.cpp:736-779`: the
OMP neighbor rebuild pays a 3-pass **serial** thread-pairs finalize and its
candidate loop reads coordinates via xtensor `operator[]` +
`frame.box.min_image()` per candidate, unlike the hoisted serial path.
Measured: Neigh 2.2× slower than LAMMPS, 35-44% of the step at scale
(`.claude/notes/cpu_md_perf_assessment_2026-06.md`) — the #1 open workstream.
**Gate:** cluster benchmarking required (see Validation); serial path must show
zero regression at T=1.

## Context (verified 2026-07-02)

- MT path `refresh_half` `src/cpu/neighbor/cell_list.cpp:621-780`, duplicated
  in `refresh_full` :788-945. Finalize = 3 serial passes: count :753-755,
  prefix :757-760, scatter through `csr_write_cursor_` :770-779. Candidate
  loop: xtensor `operator[]` coordinate reads + per-candidate
  `frame.box.min_image()` (:736-739); full 27-cell stencil + `j <= i` reject
  (:678, :728-729) — ~2× the serial path's candidate volume.
- Serial path :357-619: `__restrict` raw pointers (:401-413), half stencil
  ~14 cells (:437-439), per-cell PBC shift (:489-525, :575-580) with multi-wrap
  `while` fallback (:590-601), one-pass CSR cursor (:444-445, :610-615).
  Serial gate: `n_threads == 1 || n < n_threads` (:378).
- `Box::min_image` (frame/box.h:117-150): 3 divides + 3 `nearbyint` +
  `this`-indirection per candidate — the exact pattern CLAUDE.md documents as
  defeating hoisting (~7%/step in LJ).
- Structural fact: with `omp for` over `i`, every CSR row is produced by
  exactly one thread, contiguous in that thread's buffer in generation order.

## Plan — (b) row-owned direct finalize (recommended over minimal parallel 3-pass)

1. Candidate loop counts each row's survivors and writes `result.numneigh[i]`
   **directly** at end-of-row (what serial does at :615). Deletes the
   zero-init loop (:657-658) and finalize pass 1 (:753-755). Ghost/skipped
   rows write 0 explicitly.
2. Replace `thread_pairs_` (`vector<vector<pair<int,int>>>`, cell_list.h:237)
   with per-thread flat `vector<int>` of j-values + `vector<RowSpan>{int i;
   int count;}` row descriptors (rows contiguous by construction). Same
   monotone-capacity/`clear()` arena discipline + the 32-per-atom reserve
   heuristic (:644-651). Halves buffer bytes/pair.
3. Serial exclusive scan `numneigh → firstneigh` (O(n) ints — noise; blocked
   parallel scan is an optional follow-up per 1M-atom profiles).
4. `result.neighbors.resize(total_pairs)` — steady-state no-op (pair_list.h:28-31).
5. `omp parallel`: each thread `memcpy`s its rows' j-runs to
   `neighbors.data() + firstneigh[i]`. Race-free by row ownership; **delete
   `csr_write_cursor_`** (cell_list.h:240-245) and update `bytes()` (:953-972).
6. Apply identically to `refresh_full` — extract the finalize (ideally the
   candidate loop too, templated on half/full + tag_dedup) into one shared
   private helper so the twins can't drift.

## Candidate-loop hoisting (mirror the serial idiom)

- Hoist before `#pragma omp parallel`: `const double* __restrict x_ptr/...`,
  `const int* __restrict atom_cell_ptr/cell_count_ptr/cell_offset_ptr/binned_ptr`
  (serial block :401-413).
- Hoist `lx/ly/lz/half_l*`; replace `frame.box.min_image(...)` (:739) with the
  serial full-stencil multi-wrap `while` form (:590-601 — must be `while`, not
  single-wrap: ±2L ghost case :585-589). Makes serial/MT arithmetic identical.
- **Half stencil in MT** (biggest single win, ~2× candidate reduction; LAMMPS
  newton-on bin path). Import as a package, gated on
  `use_half = ncx_>=3 && ncy_>=3 && ncz_>=3` (:437): per-cell shift (:489-525),
  ghost-row skip (:464-467), home-cell-only tag tie-break `reject_le_i`
  (:548, :564-572). The MT everywhere-tie-break (:732-734) is only correct for
  symmetric full-stencil scans. If schedule risk: land hoisting+finalize first,
  half-stencil as a follow-up commit.
- Replace `schedule(dynamic)` chunk-1 (:669) with `schedule(dynamic, 32)` or
  `static`; A/B on the cluster.

## Invariants + coverage

1. Half-list exactness (each unordered pair within cutoff+skin exactly once,
   incl. multi-wrap cross-PBC) — pinned vs O(N²) by `CpuMdTest.CellListPairCount`
   (tests/cpu/md/test_md.cpp:197-235).
2. Row placement: consumed pairs land in rows `i < inum` after
   `trim_inum_to_local` (neighbor_list.h:219-223). Consumers (PairLJCut
   :164-227/:269-277; `PairList::for_each_pair`) need only self-consistent
   rows — serial/MT row *assignment* legitimately differs, so parity tests
   compare pair **sets**, not CSR layout.
3. Cross-rank dedup: owner-ghost pair emitted on exactly one rank (:554-572;
   lattice-alignment precondition :77-89). `test_mpi_npt_multirank_n8` is a
   documented coin-flip — not a sole gate.
4. Index bounds: `CellListShrinkingNActive` (test_md.cpp:243-273).
5. Zero steady-state allocation: `test_neighbor_list_arena.cpp:52-100`
   (pointer+capacity stability over 100 refreshes) + `test_md_no_alloc.cpp`
   (glibc-only). New arenas must be `mutable` members, same reserve heuristics.
6. `has_shifts = false` stays on the MT path (:623).

Coverage gap: ctest runs with `OMP_NUM_THREADS=8` (tests/CMakeLists.txt:17-18)
so existing tests do exercise the MT path, but there is **no explicit
serial-vs-MT parity/determinism test**. Add one: same frame,
`omp_set_num_threads(1)` vs `(8)`, compare sorted per-atom pair sets + assert
per-row j-order is thread-count-invariant (single-thread-owned rows +
contiguous scatter must preserve this).

## Size / risks / validation

~250-350 lines (cell_list.cpp ~120-160 or ~90 + ~60-line shared helper,
cell_list.h ~25, parity test ~70).

1. Half-stencil in MT changes multi-rank dedup semantics — highest risk;
   separate commit, gate with the parity test in a fake-ghost configuration
   (like `CellListMultiRankMigrationCycle`, test_md.cpp:332).
2. `refresh_full` drift — extract the shared helper.
3. Direct `numneigh[i]` writes: every row in `[0, n)` written exactly once
   (skipped ghost rows explicitly 0) or stale counts leak from the previous
   rebuild.
4. Arena warm-up transient fine; pointer-stability asserts must hold from
   refresh #2.

**Validation (cluster):** `benchmarks/md/pe/performance/polyethylene_tension.cpp`
stage-5 Neigh Matom-step/s (`stage_component_timing.csv`; baselines AV 3.77 vs
LAMMPS 8.33 single-rank; Neigh 35-44% of step at 8-64 ranks). Caveats:
(i) the 2026-06 matrix ran 1 thread/rank — production Neigh cost measured there
is the *serial* path; first re-profile stage 5 with an
`OMP_NUM_THREADS ∈ {1,2,4,8,16}` single-rank sweep to quantify the MT deficit;
serial path must show zero regression at T=1. (ii) same-allocation multi-rep
A/B — node noise ±7%.
