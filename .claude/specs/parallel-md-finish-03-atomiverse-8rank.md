---
title: Atomiverse 8-rank polyethylene_tension end-to-end run
slug: parallel-md-finish-03-atomiverse-8rank
status: approved
created: 2026-05-19
parent: parallel-md-engine.md
---

# Atomiverse 8-rank polyethylene_tension end-to-end run

## Summary

Run `mpirun -n 8 ./build/benchmarks/md/polyethylene_tension --mode full` through all 5 stages on this host, capture exit status, capture peak RSS, and commit the run artifacts (`thermo.csv`, `stress_strain.csv`, plus a small `RUN.md` describing the run) into a fixed subdirectory `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/` for downstream consumption by sub-spec 04 (scaling-report). Closes parent ac-004.

## Domain basis

Same protocol as sub-spec 02, but executed by Atomiverse via the canonical `polyethylene_tension` benchmark binary. Geometry generated procedurally by `make_pe_box()` (4×4×4 chains × 18 carbons = 1152 atoms). The "5% memory growth" gate from ac-004 is the architectural rule from spec §B.2 ("no allocation in `launch()`"): after warmup, the inner loop must hold steady-state memory.

## Design

Touched entities:

- `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/thermo.csv` — written by the binary.
- `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/stress_strain.csv` — written by the binary.
- `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/RUN.md` — new, manually authored: records the exact `mpirun` command, `git rev-parse HEAD`, OMP env, peak RSS (from `/usr/bin/time -l` on macOS), wall time, exit code.

No source code changes. The `runs/atomiverse_n8/` directory is fixed (not timestamped) so sub-spec 04 can reference stable paths. The launcher's default behavior creates timestamped paths; this sub-spec uses the launcher's `--out-dir` override path or copies the launcher's output into the canonical location.

## Files to create or modify

- `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/thermo.csv` (new)
- `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/stress_strain.csv` (new)
- `benchmarks/md/reference/polyethylene/runs/atomiverse_n8/RUN.md` (new)

## Tasks

- [ ] Ensure release build with MPI is current: `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_BENCHMARKS=ON -DATV_USE_MPI=ON` and `cmake --build build -j --target polyethylene_tension`.
- [ ] Run with explicit fixed-path output, capturing exit code, wall time, peak RSS: `mkdir -p benchmarks/md/reference/polyethylene/runs/atomiverse_n8 && /usr/bin/time -l mpirun -n 8 ./build/benchmarks/md/polyethylene_tension --mode full --out-dir benchmarks/md/reference/polyethylene/runs/atomiverse_n8 2> benchmarks/md/reference/polyethylene/runs/atomiverse_n8/time.log` and confirm exit 0.
- [ ] Verify both CSVs exist and are non-empty under `runs/atomiverse_n8/`; each ≥ 11 lines (header + ≥ 10 data rows).
- [ ] Compute peak RSS from `time.log` (look for `maximum resident set size`) and the wall time (look for `real`).
- [ ] Author `RUN.md` with: the exact command, `git rev-parse HEAD`, `uname -mrs`, `sysctl -n machdep.cpu.brand_string`, `OMP_NUM_THREADS`/`OMP_PROC_BIND` settings used, wall time, peak RSS, exit code (0), and the binary's printed final per-step wall time.
- [ ] Re-run the same command a second time and confirm peak RSS does not exceed 5% of the first run's number (steady-state inner-loop allocation gate from parent ac-004). Record both numbers in `RUN.md`.
- [ ] Commit the three files in one commit: `bench(parallel-md-finish-03): Atomiverse 8-rank end-to-end PE-tension run`.

## Testing strategy

- **Happy path:** exit 0, both CSVs land, RSS diff between two consecutive runs ≤ 5%.
- **Edge case:** if `mpirun -n 8` oversubscribes on a host with fewer than 8 performance cores, retry with `OMP_NUM_THREADS=1 mpirun -n 8 --oversubscribe …` and document in `RUN.md`.
- **Edge case:** if Open-MPI complains about insufficient slots, set `OMPI_MCA_rmaps_base_oversubscribe=1`.
- **Domain validation:** end-of-stage-1 T should be within 50 K of 500 K (verify by reading the `T_K` column at step ≈ stage1_steps); end-of-stage-4 T should be within 30 K of 100 K. Final `strain` in stress_strain.csv ≥ 0.3.

## Out of scope

- Comparing Atomiverse output to LAMMPS output (sub-spec 04).
- Scaling beyond single-node 8 ranks.
- Performance-tuning OMP placement (sub-spec 05).
- Editing `polyethylene_tension.cpp` source.
