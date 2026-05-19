---
title: Per-component profile CSVs at 1 and 8 ranks
slug: parallel-md-finish-01-profile-csvs
status: retired
created: 2026-05-19
last-updated: 2026-05-19
parent: parallel-md-engine.md
---

# Per-component profile CSVs at 1 and 8 ranks

## Summary

Run the existing `benchmarks/md/pe_profile_single` binary on this host (Apple Silicon, Homebrew Open-MPI) at 1 rank and 8 ranks of the 20-chain × 1000-monomer PE system for 1000 timed steps after a 100-step warm-up. Commit the two resulting CSVs (`profile_1.csv`, `profile_8.csv`) plus their stacked-bar PNGs (`profile_1.png`, `profile_8.png`) into `benchmarks/md/reference/polyethylene/`. No source code changes — this sub-spec is an acceptance-verification run that closes parent ac-008 (binary + timers exist), ac-009 (CSVs committed), ac-010 (plot script works), and produces the input data that downstream sub-specs (04) consume for scaling-report ratios.

## Domain basis

No new physics. The protocol matches Phase 5 spec §B.1 — stage-5 tension at 300 K, Δt = 0.5 fs, NVT-only (NoseHoover τ = 100 fs), to isolate force-kernel timing from barostat/deform overhead. The 20-chain × 1000-monomer (`4×4×4` chains × 18 UA carbons per chain = 1152 atoms) size matches MSU's small validation box; large enough that LJ neighbor walks dominate, small enough that runs finish in seconds.

## Design

Touched entities:

- `benchmarks/md/reference/polyethylene/profile_1.csv` — new artifact, written by binary
- `benchmarks/md/reference/polyethylene/profile_8.csv` — new artifact
- `benchmarks/md/reference/polyethylene/profile_1.png` — new artifact, written by `tools/plot_profile.py`
- `benchmarks/md/reference/polyethylene/profile_8.png` — new artifact
- `.gitattributes` (optional) — register `*.csv` as `text eol=lf` so CRLF doesn't corrupt header line

No new C++ code; no CMake changes. The binary already accepts `--out PATH`, `--chains-per-axis N`, `--n-atoms-per-chain N`, `--steps N`, `--warmup N`. MPI launch via Homebrew `mpirun -n <ranks> ...`.

## Files to create or modify

- `benchmarks/md/reference/polyethylene/profile_1.csv` (new)
- `benchmarks/md/reference/polyethylene/profile_8.csv` (new)
- `benchmarks/md/reference/polyethylene/profile_1.png` (new)
- `benchmarks/md/reference/polyethylene/profile_8.png` (new)

## Tasks

- [x] Configure release build with benchmarks + MPI: `cmake -B build -DCMAKE_BUILD_TYPE=Release -DATV_BUILD_BENCHMARKS=ON -DATV_BUILD_TESTS=ON -DATV_USE_MPI=ON` and build `cmake --build build -j --target pe_profile_single`
- [x] Run 1-rank profile: `./build/benchmarks/pe_profile_single --chains-per-axis 4 --n-atoms-per-chain 18 --steps 1000 --warmup 100 --out benchmarks/md/reference/polyethylene/profile_1.csv` and verify exit-0 and CSV present with ≥ 10 data rows (resolved: binary actually lives at `./build/benchmarks/pe_profile_single`, not `…/benchmarks/md/`)
- [x] Run 8-rank profile: `/opt/homebrew/bin/mpirun -n 8 --oversubscribe ./build/benchmarks/pe_profile_single --chains-per-axis 4 --n-atoms-per-chain 18 --steps 1000 --warmup 100 --out benchmarks/md/reference/polyethylene/profile_8.csv` and verify same (used `--oversubscribe` because Open-MPI sees 10 logical / 8 performance cores on this Apple-Silicon host)
- [x] Generate PNGs: `python3 tools/plot_profile.py benchmarks/md/reference/polyethylene/profile_1.csv` and same for `profile_8.csv`
- [x] Spot-check both CSVs: profile_1.csv `LJ` row total_ms=154.79 (non-zero); profile_8.csv `HaloExchangeForward` row total_ms=118.98 with n_calls=1000 (non-zero, MPI collectives ran); pct_of_step sums to 99.77% (1r) and 99.80% (8r) — both within [90, 105].
- [x] Commit all four artifacts in one commit message `bench(parallel-md-finish-01): commit profile_{1,8}.csv + PNGs`

## Testing strategy

- **Happy path:** all four artifacts land, file sizes > 0, headers exactly `component,total_ms,pct_of_step,n_calls,mean_us`.
- **Edge case:** if `mpirun` complains about oversubscription on Apple Silicon (8 ranks on 8/10/12 cores), retry with `--oversubscribe` or pin to performance cores via `taskset` equivalent (`OMP_NUM_THREADS=1` is the simplest mitigation — record which path was taken in the commit message).
- **Domain validation:** sanity check pct_of_step column sums to ≈ 100% (allow ±5% from un-timed gaps); LJ::launch + NeighborList::* + HaloExchange* should account for ≥ 60% of step at 1 rank, ≥ 70% at 8 ranks.

## Out of scope

- 32-rank profile (cluster-only; tracked under `parallel-md-finish-06` as documented-pending).
- Editing `pe_profile_single.cpp` source — if a needed timer is missing, file a separate spec.
- LAMMPS reference data (sub-spec 02).
- Translating CSV numbers into ratios / scaling_report.md (sub-spec 04).
