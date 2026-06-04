# CPU MD ↔ LAMMPS performance assessment — decision report (2026-06)

Deliverable of spec `cpu-md-lammps-perf-assessment` (ac-004). Question posed:
*is Atomiverse's CPU MD architecture wrong / far behind LAMMPS, or just
under-optimized — refactor or optimize?* System: MSU amorphous PE, 10 000 UA
atoms, the `lammps.in` 5-stage protocol. All runs 1 thread/rank, MPICH 4.3.2,
`mpiexec`. Data under `benchmarks/md/reference/polyethylene/perf_assessment_2026-06-04/`.

## TL;DR — **OPTIMIZE, do not refactor.**

Single-rank Atomiverse is **~1.8× slower** than LAMMPS — an *implementation*
gap (kernels 1.2–2.2× off, not 10×), and AV already uses the same core
algorithms LAMMPS does (half neighbor lists + Newton's 3rd law, spatial
decomposition, per-atom topology, one batched `MPI_Allreduce`). The architecture
is sound. The 32→64 non-scaling is communication + load imbalance, which is the
*regime LAMMPS is also bound by here* (over-decomposition: 156 atoms/rank).

## Evidence

### 1. Single-rank kernel throughput (Matom-step/s, higher = faster)

| Section | LAMMPS 1r | Atomiverse 1r | AV slowdown |
|---------|-----------|---------------|-------------|
| Pair    | 2.74      | 1.38          | **2.0×**    |
| Bonded  | 12.24     | 10.03         | 1.2×        |
| Neigh   | 8.33      | 3.77          | **2.2×**    |
| Overall | 1.638     | 0.906         | **1.8×**    |

LAMMPS: `lammps_perf.in` (nve/limit+langevin, 5000 steps). AV: `av_1r` stage 5
(17 180 steps), per-section seconds from `stage_component_timing.csv`. Both
single rank, 10 000 atoms. The gap is kernel efficiency — **not** algorithmic:
AV's neighbor list is half + newton (`NeighborList::refresh_half`,
`cell_list.cpp half_stencil_`), exactly like LAMMPS (`attributes: half, newton on`).

### 2. Per-component attribution vs rank count (AV, stage 5 %)

| ranks | Pair | Bonded | **Neigh** | Comm | nlocal min/max (ave 156) |
|-------|------|--------|-----------|------|--------------------------|
| 1     | 65.7 | 9.0    | **24.0**  | 0.0  | 10000/10000              |
| 8     | 36.6 | 9.6    | **34.9**  | 17.2 | 1110/1374                |
| 32    | 28.9 | 12.2   | **44.3**  | 10.4 | 90/379                   |
| 64    | 19.4 | 14.2   | **39.6**  | 20.9 | 16/195                   |

LAMMPS for comparison (its own breakdown, stage 5): 32r Comm 60% / Neigh 17% /
Pair 16%; 64r Comm 65% / Neigh 17% / Pair 11%; nlocal min/max 4/345 (64r).

### 3. Whole-run wall (5 stages)

LAMMPS 32r=164s / 64r=86s (1.9× strong scaling). AV 32r=324s / 64r=336s (no
scaling). av/lmp ≈ 2.0× (32r) / 3.9× (64r).

## Diagnosis — algorithmic vs implementation

| Gap | Classification | Evidence |
|-----|----------------|----------|
| **Neighbor build 2.2× slow + 35–44% of step at scale** (vs LAMMPS 17%) | **Implementation** | half-list algo already correct; cost is the rebuild itself — matches `cpu_md_perf_review_2026-05-24` HIGH#3 (`cell_list.cpp` 3-pass `emplace_back(i,j)` into thread-local `vector<pair>` then count/prefix/scatter). |
| **Pair 2.0× slow** | **Implementation** | ket-loop is well-engineered but the exclusion-filtered CSR + scalar reorder accumulation + 5-pass per-evaluator force reduction add overhead (perf-review HIGH#2, LOW LJ items). |
| **No 32→64 strong scaling** | **Architecture-adjacent (shared with LAMMPS)** | Comm rises to ~21%; load imbalance nlocal 16↔195 at 64r. LAMMPS is *itself* 65% Comm here — 156 atoms/rank is past the useful-decomposition point. |
| Bonded 1.2× | none material | already near parity. |

## Recommendation — ranked optimization workstream (no refactor)

1. **Neighbor-list rebuild (biggest win).** Replace the 3-pass
   `thread_pairs_.emplace_back` build in `cell_list.cpp` with a parallel
   direct-CSR bin scatter (perf-review HIGH#3). Expect to close most of the
   2.2× Neigh gap → it is the largest single cost at scale.
2. **Fused force arena + single reduction** (perf-review HIGH#2/#5): one shared
   thread-local force scratch + one `ForceReduce` instead of 5 per-evaluator
   zero+reduce passes; also fixes the per-rank global-capacity over-alloc.
   Targets the Pair/Bonded reduction overhead.
3. **Resolved-local-index bonded loops** (perf-review HIGH#1): drop the
   per-tuple `tag_to_local.find()` from bonded hot loops.
4. **Load balance for scaling**: `rcb_balancer` exists — verify it is enabled
   for this workload and improve it; nlocal 16↔195 shows the static slab
   partition is poorly balanced for the amorphous blob. This is the lever for
   32→64. (Diminishing returns past ~32 ranks for a 10k-atom system regardless;
   the honest scaling target is "match LAMMPS's 1.9× to 64r", not linear.)

Realistic target after 1–3: ~1.2–1.4× of LAMMPS single-node throughput (from
1.8×). Item 4 restores strong scaling toward LAMMPS's curve.

## Caveats / method notes

- `--scale 0.1` runs (≈33k steps) used for the matrix; per-step fractions are
  scale-independent (verified consistent across stages). Absolute kernel
  throughput taken from the longest sample (stage 5 / 5000-step LAMMPS run).
- Two bad compute nodes (n215/n217) hung MPI init mid-assessment; reran on
  healthy nodes. Login-node test runs need `FI_PROVIDER=tcp` (see build memory).
- Instrumentation is additive (timing only); physics unchanged by construction
  (ac-005) — thermo on every run reproduced the expected T/stage behaviour.
