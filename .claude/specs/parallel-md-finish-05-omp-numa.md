---
title: OMP_PROC_BIND measurement + README Thread binding section
slug: parallel-md-finish-05-omp-numa
status: approved
created: 2026-05-19
parent: parallel-md-engine.md
---

# OMP_PROC_BIND measurement + README Thread binding section

## Summary

On this host (Apple Silicon, single-socket — no NUMA), run the small profile binary twice — once with `OMP_PROC_BIND` unset (default), once with `OMP_PROC_BIND=close OMP_PLACES=cores` — and record both wall-times in a new "Thread binding" subsection of `benchmarks/md/reference/polyethylene/README.md`. If the delta is below noise (≤ 2% — likely on single-socket Apple Silicon), record that explicitly with a "single-socket — no NUMA effect observed" note as the parent acceptance permits. Closes parent ac-019.

## Domain basis

Pure measurement-and-documentation work. The C++ stderr warning ("OMP_PROC_BIND unset and n_threads > 1") is already implemented and verified by parent ac-018; this sub-spec only adds the numerical measurement that parent ac-019 demands.

## Design

Touched entities:

- `benchmarks/md/reference/polyethylene/README.md` — append a "## Thread binding" subsection between "Reproducing the Atomiverse run" and any later section. Records the measurement command, both wall times, the host CPU brand, and the conclusion.

No source code changes.

## Files to create or modify

- `benchmarks/md/reference/polyethylene/README.md` (modify)

## Tasks

- [ ] Confirm the release build is current: `cmake --build build -j --target pe_profile_single`.
- [ ] Run measurement A (unset): `unset OMP_PROC_BIND OMP_PLACES; OMP_NUM_THREADS=4 /usr/bin/time -p ./build/benchmarks/md/pe_profile_single --steps 200 --warmup 50 --out /tmp/profile_unset.csv 2>&1 | tee /tmp/time_unset.log`.
- [ ] Run measurement B (bound): `OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=4 /usr/bin/time -p ./build/benchmarks/md/pe_profile_single --steps 200 --warmup 50 --out /tmp/profile_bound.csv 2>&1 | tee /tmp/time_bound.log`.
- [ ] Extract the `real` line from each `time.log` file; record both wall times.
- [ ] Compute the relative delta `(unset - bound) / unset`. If `|delta| < 0.02` (2%), the conclusion is "single-socket — no NUMA effect observed"; otherwise record the actual numbers and the bound version's relative speedup.
- [ ] Append a "## Thread binding" subsection to `README.md` with: the exact two commands, host CPU brand (`sysctl -n machdep.cpu.brand_string`), number of physical cores (`sysctl -n hw.physicalcpu`), the two wall times, the delta, and the conclusion.
- [ ] Commit: `docs(parallel-md-finish-05): record OMP_PROC_BIND wall-time measurement on dev box`.

## Testing strategy

- **Happy path:** both runs complete cleanly; README has a new "Thread binding" subsection with both numerical values.
- **Edge case:** if `OMP_PROC_BIND=close` fails to apply (e.g., Apple Silicon Clang LLVM-OpenMP without `OMP_PLACES` support), document the failure and the fallback ("setting honoured by GNU libgomp on Linux only — not measurable on this host's libomp"). The parent acceptance permits this outcome.
- **Domain validation:** the difference, in either direction, should be small (< 5%) on single-socket Apple Silicon. A large delta would be suspicious and warrants re-running with `--steps 1000` for noise reduction.

## Out of scope

- Editing `cpu::Backend` to set `OMP_PROC_BIND` programmatically (parent decision: warn but don't override).
- NUMA measurement on multi-socket Linux (out of host scope).
- Comparing thread-binding effect on LAMMPS (out of scope; LAMMPS owns its own threading defaults).
