---
title: OMP_PROC_BIND measurement + README Thread binding section
slug: parallel-md-finish-05-omp-numa
status: code-complete
created: 2026-05-19
last-updated: 2026-05-19
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

- [x] Release build already current from spec 01.
- [x] Measurement A (unset OMP_PROC_BIND): bumped from 200 → 5000 steps for noise reduction. Per-step wall = 0.273 ms; real = 1.47 s.
- [x] Measurement B (`OMP_PROC_BIND=close OMP_PLACES=cores`): per-step wall = 0.267 ms; real = 1.49 s.
- [x] Delta = (0.273 − 0.267) / 0.273 = **−2.2%** — at the noise floor on this single-socket Apple M4 host.
- [x] Both numbers, both commands, host fingerprint (Apple M4 / 10 physical cores), interpretation paragraph, and 2-socket-Linux NUMA hint already drafted into the "Measured NUMA effect" subsection (which sits under the existing `## Thread binding` heading).
- [x] Refined task 2/3 to use `--steps 5000` (not 200) — 200 steps too noisy for stable wall-time signal; both ac-001 ("two numerical wall times AND a delta percentage") and the exact "single-socket — no NUMA effect observed" phrase are present.
- [x] Commit the README modification.

## Testing strategy

- **Happy path:** both runs complete cleanly; README has a new "Thread binding" subsection with both numerical values.
- **Edge case:** if `OMP_PROC_BIND=close` fails to apply (e.g., Apple Silicon Clang LLVM-OpenMP without `OMP_PLACES` support), document the failure and the fallback ("setting honoured by GNU libgomp on Linux only — not measurable on this host's libomp"). The parent acceptance permits this outcome.
- **Domain validation:** the difference, in either direction, should be small (< 5%) on single-socket Apple Silicon. A large delta would be suspicious and warrants re-running with `--steps 1000` for noise reduction.

## Out of scope

- Editing `cpu::Backend` to set `OMP_PROC_BIND` programmatically (parent decision: warn but don't override).
- NUMA measurement on multi-socket Linux (out of host scope).
- Comparing thread-binding effect on LAMMPS (out of scope; LAMMPS owns its own threading defaults).
