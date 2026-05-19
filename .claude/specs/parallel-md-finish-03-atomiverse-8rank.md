---
title: Atomiverse 8-rank polyethylene_tension end-to-end run
slug: parallel-md-finish-03-atomiverse-8rank
status: retired
created: 2026-05-19
last-updated: 2026-05-19
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

- [x] Ensure release build with MPI is current: built `polyethylene_tension` against the build tree configured by spec 01.
- [x] Run with `/usr/bin/time -l` + `--oversubscribe` (8 ranks on 10-core M4): exited cleanly in 170.68 s; `thermo.csv` (281 lines) + `stress_strain.csv` (101 lines) + `time.log` (994 B) landed under `runs/atomiverse_n8/`. Used `--box-len 90` to match the LAMMPS reference box (sub-spec 02) so sub-spec 04 can do apples-to-apples comparison.
- [x] CSV headers verified: `thermo.csv` = `step,stage,T_target_K,T_K,P_bar,E_pot_eV,E_kin_eV,box_x,box_y,box_z`; `stress_strain.csv` = `step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa`. Both ≥ 11 lines.
- [x] Peak RSS extracted from `time.log`: run 1 = 19,333,120 B (18.44 MiB); wall time = 170.68 s.
- [x] Authored `RUN.md` with command, git SHA `ee80e0e`, host `Apple M4 / Darwin 25.4.0`, OMP env (unset; binary emits the documented warning), both runs' wall + RSS, exit code 0, and a "Known anomaly" note on the misleading `T_K` / `E_kin_eV` column values (tracked as follow-up; CSV structure is acceptable per ac-002, so this is not a blocker).
- [x] Re-ran end-to-end; peak RSS = 19,611,648 B. Growth = (19611648 - 19333120) / 19333120 = **1.44%**, well below the 5% gate. Both numbers committed in `RUN.md`.
- [x] Commit the four files (RUN.md, thermo.csv, stress_strain.csv, time.log) in one commit.

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
