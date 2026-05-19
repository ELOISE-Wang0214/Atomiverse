---
title: scaling_report.md — single-source-of-truth for Phase 5 perf numbers
slug: parallel-md-finish-04-scaling-report
status: approved
created: 2026-05-19
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

- [ ] Capture compiler version + active flags from the existing build: `cmake --build build --target atomiverse -v 2>&1 | head -1` (the first compile line shows the full `clang++ -O3 -march=native … -ffast-math …` invocation). Save as `BUILD_FLAGS` for use in the report.
- [ ] Capture host fingerprint: `sysctl -n machdep.cpu.brand_string`, `sysctl -n hw.ncpu`, `sw_vers -productVersion`. Save as `HOST_FP`.
- [ ] Find pre-B.2 and pre-B.3 baseline numbers from git history: `git log --oneline -- src/cpu/md/lj.cpp src/cpu/potentials/harmonic_bond.cpp | head -20`. The commits where B.2 (1.21× NEON) and B.3 (HB 1.13×) landed are recorded in the parent spec status table — quote those ratios as the documented baseline-to-post comparison.
- [ ] Extract per-component numbers from `profile_1.csv` (committed by sub-spec 01): `awk -F, '$1=="LJ::launch" || $1=="HarmonicBond::launch" || $1=="NeighborList::build_half"' benchmarks/md/reference/polyethylene/profile_1.csv`. Record `total_ms` and `n_calls` for each.
- [ ] Extract Atomiverse wall time from `runs/atomiverse_n8/RUN.md` and `runs/atomiverse_n1/RUN.md` (if a single-rank run was also captured — if not, run it here as well via the launcher and commit alongside).
- [ ] Extract LAMMPS wall time + rebuild count from `runs/lammps_n1/log.lammps` and `runs/lammps_n8/log.lammps` (the LAMMPS log records the neighbor-rebuild count in lines matching `^Neighbor list builds = `).
- [ ] Author `scaling_report.md` with these sections:
  - Front matter: host fingerprint, compiler/flag string, git SHA at measurement.
  - Table 1 "Kernel speedups": rows for `LJ::launch (scalar baseline / xsimd current / ratio)`, `HarmonicBond::launch (unpacked baseline / stack-packed current / ratio)`, with `total_ms` from profile_1.csv and historical baselines quoted from the parent spec status table.
  - Table 2 "Atomiverse vs LAMMPS": rows for `single-core PE-tension, Δt=0.5fs, 1000 steps` and `8-rank PE-tension` with `Atomiverse wall time / LAMMPS wall time / ratio`.
  - Table 3 "Neighbor-list rebuild cadence": `Atomiverse rebuild count over 1000 PE-tension steps (from profile_1.csv `NeighborList::build_half` n_calls) / LAMMPS rebuild count (from log.lammps "Neighbor list builds")`.
  - Subsection "Hotspot summary": top-3 components by `pct_of_step` at 1 rank and 8 ranks.
  - Subsection "Structural-gap analysis": included if single-core ratio > 1.2× — names the remaining bottleneck (citing the matching `profile_1.csv` row) and the architectural cost of closing it (e.g., "LAMMPS Verlet-list bin packing not yet ported, would require a `cell_list` rewrite; deferred to Phase 6").
  - Subsection "Build configuration": pasted `BUILD_FLAGS` line plus the `-ffast-math` per-source allowlist from `CMakeLists.txt`.
- [ ] Commit `scaling_report.md` and a tiny README cross-reference: `docs(parallel-md-finish-04): scaling report with kernel speedups + LAMMPS ratios`.

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
