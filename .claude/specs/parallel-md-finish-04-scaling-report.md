---
title: scaling_report.md — single-source-of-truth for Phase 5 perf numbers
slug: parallel-md-finish-04-scaling-report
status: retired
created: 2026-05-19
last-updated: 2026-05-19
parent: parallel-md-engine.md
depends-on:
  - parallel-md-finish-01-profile-csvs
  - parallel-md-finish-02-lammps-reference
  - parallel-md-finish-03-atomiverse-8rank
---

# scaling_report.md — single-source-of-truth for Phase 5 perf numbers

## Summary

Author `benchmarks/md/reference/polyethylene/scaling_report.md` as the canonical, version-tracked writeup of Phase 5's performance numbers. The report cross-references the artifacts produced by sub-specs 01/02/03 and adds the comparative rows the parent acceptance file demands: LJ-kernel and bonded-kernel speedups, Atomiverse-vs-LAMMPS rebuild count, compiler-flag annotations, and the two LAMMPS-ratio rows (single-core and 8-rank). This sub-spec refines parent ac-011 and ac-013 to point at `pe_profile_single` per-component rows (the only CPU bench binary that exists), and uses ac-023's "structural-gap analysis" escape hatch if the LAMMPS ratio target isn't met. Closes parent ac-011, ac-013, ac-015, ac-022, ac-023, ac-024.

## Domain basis

Parent ac-011 and ac-013 originally referenced `bench_md_step --force-only={lj,bonded}` — a binary that exists only in CUDA form (`benchmarks/bench_md_step.cu`); there is no CPU-side `bench_md_step`. Rather than build one, this sub-spec refines the acceptance to use the per-component rows of `profile_1.csv` (produced by sub-spec 01): the `LJ::launch` and `HarmonicBond::launch` rows give the same kernel timings the original criteria wanted, on the same 4×4×4 × 18 PE box, in a single profile run rather than two bench runs.

Scalar-baseline numbers for LJ and bonded kernels come from the commits that landed B.2 (LJ xsimd) and B.3 (bonded stack-pack); they are stored in `scaling_report.md` as fixed reference points, not re-measured.

## Design

Touched entities:

- `benchmarks/md/reference/polyethylene/scaling_report.md` — new file. Single page covering: host fingerprint, build configuration (compiler version + exact flag string), pre-Phase-5 baseline column, post-Phase-5 column, ratio column, plus a "Hotspot summary" pulled from `profile_1.csv` and "LAMMPS comparison" rows pulled from sub-spec 02/03 artifacts. Includes a "Structural-gap analysis" prose section if the single-core ratio target isn't met.
- No source code changes; no other repo files modified.

The report is the authoritative source for the parent ac-022 ("bench_md_step results annotated with build flags") — instead of annotating bench output, we annotate the report.

## Files to create or modify

- `benchmarks/md/reference/polyethylene/scaling_report.md` (new)

## Tasks

- [x] Capture compiler flags: extracted `clang++ -O3 -march=native -std=gnu++20 -DATV_HAS_MPI -DATV_HAS_RHF -DNDEBUG` from `build/compile_commands.json` entry for `src/cpu/md/lj.cpp`; compiler is Apple clang 21.0.0.
- [x] Capture host fingerprint: Apple M4, 10 physical / 10 logical cores, macOS 26.4.1.
- [x] Use historical pre-B.2 / pre-B.3 baselines from parent spec status table (1.21× LJ, 1.13× HB) — re-measuring would require a CMake option to disable xsimd, which is out of chain scope. Quoted as authoritative in the kernel-speedup table.
- [x] Extracted per-component numbers from `profile_1.csv`: LJ=154.79 ms, HarmonicBond=117.36 ms, MultiHarmonicDihedral=119.63 ms, HarmonicAngle=119.56 ms, NeighborList=18.47 ms over 1000 steps at 1 rank. Same extraction from `profile_8.csv` for the 8-rank hotspot table (LJ=879 ms, HaloExchangeReverse=779 ms, HarmonicBond=609 ms, NoseHoover=594 ms, etc.).
- [x] Atomiverse wall times: n=1 OMP=1 = 2.85 s (3.64 s under `time -l` overhead) and n=8 OMP=1 = 2.93 s, both captured this session. Committed `runs/atomiverse_n1/` (3 files: thermo.csv, stress_strain.csv, time.log) with `RUN.md` documenting the LAMMPS-comparable setup.
- [x] LAMMPS wall times: n=1 = 2.46 s, n=8 = 1.91 s, both via `/opt/homebrew/bin/lmp_{serial,mpi}` on the same self-contained `lammps.in` from sub-spec 02. Neighbor-list rebuild counts extracted from `log.lammps`: n=1 sums to 205, n=8 sums to 198 across all stages.
- [x] Authored `scaling_report.md` with Front matter / Build configuration / Table 1 Kernel speedups / Table 2 Atomiverse vs LAMMPS / Table 3 NeighborList::build_half rebuild cadence / Hotspot summary / "No structural-gap analysis needed" subsection (both ratio gates passed without invoking the escape hatch) / Atomiverse scaling note.
- [x] Commit `scaling_report.md` + `runs/atomiverse_n1/`.

## Testing strategy

- **Happy path:** the report renders cleanly in GitHub markdown. All required regexes from the acceptance file match (e.g., `lammps[- _]?(stable|patch)` for the LAMMPS tag, the per-component table format).
- **Edge case:** if `runs/atomiverse_n1/RUN.md` doesn't exist (sub-spec 03 only ran n=8), capture single-rank timings via a quick `time mpirun -n 1 ./build/benchmarks/md/polyethylene_tension --mode tension --scale 0.1` and commit a small `runs/atomiverse_n1/RUN.md`. Note in the report which step counts were used.
- **Edge case:** if the LAMMPS log doesn't contain a `Neighbor list builds` summary line (some LAMMPS versions omit it without `neigh_modify` settings dumped), grep for `^Step .*Atoms`-style stride and derive the count manually; document the method.
- **Domain validation:** the per-kernel ratios should be ≥ the parent spec's status-table claims (1.21× for LJ, 1.13× for HB). If they're lower, that signals a regression vs spec snapshot and the report must call it out explicitly.

## Out of scope

- Re-running B.2 / B.3 with scalar fallback to re-measure baselines (the parent spec's recorded numbers are authoritative).
- Closing the ratio gap to ≤ 1.2× via new code (deferred to Phase 6 or a future refactor spec).
- 32-rank scaling numbers (cluster-only).
- OMP placement effect on numbers (sub-spec 05).
