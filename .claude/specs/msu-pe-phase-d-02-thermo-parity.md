---
title: MSU PE validation Phase D.2 — per-step thermo parity vs LAMMPS
status: in-progress
created: 2026-05-20
parent: msu-pe-validation.md
---

# MSU PE validation Phase D.2 — per-step thermo parity vs LAMMPS

## Summary

This sub-spec delivers the per-step thermo-comparison portion of `tools/validate_polyethylene.py`: a pure-Python analysis tool that reads the LAMMPS reference thermo CSV and the Atomiverse thermo CSV produced by chain link `msu-pe-phase-d-01-reference-runs`, and over the first 1000 steps of the MSU PE-tension protocol asserts cross-engine agreement step-by-step — 1e-4 relative on `T`, `E_pot`, `E_kin` and 1e-3 relative on `Pxx`, `Pyy`, `Pzz`. The check runs for both the 1-rank and the 8-rank CSV pairs. Invoked as `validate_polyethylene.py --check thermo`, it exits non-zero on any tolerance breach and prints the offending step number and the max deviation for each breached column, so it is a machine-verifiable parity gate once link 01's reference artifacts are committed. No C++ or production-source change.

## Domain basis

Over the first ~1000 steps of the deterministic regime, two correct MD engines integrating the same initial state with the same potential and timestep track each other to integration-order floating-point error, before thermostat/barostat stochastics and chaotic Lyapunov divergence dominate. This justifies the tight 1e-4 relative bound on the smoothly-varying state quantities `T`, `E_pot`, `E_kin`. Instantaneous pressure components `Pxx/Pyy/Pzz` are derived from the virial and carry larger step-to-step fluctuations, so they get the looser 1e-3 relative bound. Both tolerances and the 1000-step window are taken verbatim from `.claude/specs/msu-pe-validation.md` Phase D.2. Units follow the project convention (T in K, energies in eV, pressure in bar); the comparison is unit-agnostic because it is purely relative, but both CSVs must report each column in the same unit — a constraint owned by chain link 01 and re-stated under Testing strategy.

## Design

`tools/validate_polyethylene.py` already exists from earlier work and contains a generic block-average thermo check plus stress-strain analysis. This sub-spec **replaces the thermo-comparison logic** with a per-step comparison and adds a `--check` mode selector so D.3 and D.4 can extend the same file. Entities:

- `ThermoTolerances` — frozen dataclass holding the six per-column relative thresholds; defaults are the D.2 values (`T/E_pot/E_kin = 1e-4`, `Pxx/Pyy/Pzz = 1e-3`). Replaces the `--strict`-gated tolerance mutation; D.2 tolerances are now the contract default for the thermo check.
- `StepDeviation` — frozen dataclass `{column, step, rel_err}` describing the single worst per-step breach for one column.
- `compare_thermo_per_step(av_rows, ref_rows, tol, n_steps)` — pure function: aligns the two CSVs by the `step` column over the first `n_steps`, computes per-step relative error for each of the six columns, returns the list of `StepDeviation` (one worst entry per breached column). Empty list ⇒ pass.
- `run_thermo_check(pairs, tol, n_steps)` — orchestrates the comparison over a list of `(label, av_csv, ref_csv)` pairs (label is `"1rank"` / `"8rank"`); aggregates breaches into a `Report`.
- `main` — gains a `--check {thermo,all}` argument; `--check thermo` resolves the four reference CSVs under `benchmarks/md/reference/polyethylene/` (two ranks × {atomiverse, lammps}) and runs `run_thermo_check`. Existing stress-strain code is left intact under `--check all` so D.4 can wire to it later.

Step alignment is by explicit `step` value, not row index, so that a mismatched logging cadence between engines is detected (a missing/extra step is itself a failure) rather than silently skewing the comparison. Missing input CSVs exit code 2; tolerance breach exits code 1; clean exits 0.

## Files to create or modify

- `tools/validate_polyethylene.py` — modify: replace `_check_thermo` and `Tolerances` thermo fields with per-step logic; add `ThermoTolerances`, `StepDeviation`, `compare_thermo_per_step`, `run_thermo_check`; add `--check` mode and the link-01 CSV path resolution to `main`.
- `tools/test_validate_polyethylene.py` (new) — pytest module: unit tests for `compare_thermo_per_step` (within-tolerance pass, single-step breach detection with correct step/deviation, P vs T threshold separation, step-misalignment failure) and an integration test invoking `main(["--check","thermo", ...])` against tiny synthetic CSV fixtures.
- `tools/testdata/pe_thermo/` (new) — small synthetic CSV fixtures used by the test module: a matching pair (passes) and a breaching pair (fails at a known step), for both `1rank` and `8rank` labels.

## Tasks

- [x] Write failing tests for `compare_thermo_per_step` (`tools/test_validate_polyethylene.py`)
- [x] Add synthetic thermo CSV fixtures — matching pair and breaching pair (`tools/testdata/pe_thermo/`)
- [x] Implement `ThermoTolerances` and `StepDeviation` dataclasses in `tools/validate_polyethylene.py`
- [x] Implement `compare_thermo_per_step` (step-aligned per-step relative-error scan) in `tools/validate_polyethylene.py`
- [x] Implement `run_thermo_check` aggregating 1-rank and 8-rank pairs in `tools/validate_polyethylene.py`
- [x] Add `--check {thermo,all}` mode and link-01 reference-CSV path resolution to `main` in `tools/validate_polyethylene.py`
- [x] Add docstrings with units and exit-code semantics to all new functions per project doc style
- [x] Write failing integration test invoking `main(["--check","thermo", ...])` on fixtures (`tools/test_validate_polyethylene.py`)
- [x] Run full check + test suite

## Testing strategy

- **Happy path** — `compare_thermo_per_step` on a matching CSV pair (relative error below threshold on every step) returns an empty breach list; `main(["--check","thermo"])` on matching fixtures exits 0.
- **Breach detection** — a fixture pair that deviates above tolerance at exactly one known step yields a single `StepDeviation` naming that step and the correct `rel_err`; `main` exits 1 and prints the step number and max deviation.
- **Threshold separation** — a deviation of ~5e-4 on `Pxx` passes (under 1e-3) while the same magnitude on `T` fails (over 1e-4), confirming the two tolerance bands are applied per column.
- **Step misalignment** — a reference CSV missing a step present in the Atomiverse CSV (or vice versa) over the first 1000 steps is reported as a failure, not silently skipped.
- **Missing input** — `main(["--check","thermo"])` with absent reference CSVs exits 2.
- **Domain validation** — once link 01 commits the real reference CSVs, `validate_polyethylene.py --check thermo` passes for both the 1-rank and 8-rank pairs, confirming Atomiverse and LAMMPS agree to 1e-4 (T/E_pot/E_kin) and 1e-3 (Pxx/Pyy/Pzz) per step over the first 1000 steps. This run is **blocked until `msu-pe-phase-d-01-reference-runs` commits its four reference CSVs under `benchmarks/md/reference/polyethylene/`** — a hard chain dependency, not a deferred nicety. The unit and integration tests above use synthetic fixtures and are not blocked.

## Out of scope

- The reference runs themselves (LAMMPS and Atomiverse thermo CSV generation) — owned by chain link `msu-pe-phase-d-01-reference-runs`.
- Equilibrium-state observables: density, RDF peak, end-to-end distance — Phase D.3, a later extension of the same file.
- Stress-strain curve analysis: elastic modulus, yield stress, hardening slope — Phase D.4; the existing `_check_stress_strain`/`_compute_stress_strain_stats` code is left untouched and runs only under `--check all`.
- Plotting / visualization of the thermo traces — not required for the parity gate.
- Any non-stdlib dependency in the script itself; the script stays stdlib-only (`csv`, `argparse`). pytest is used only by the test module, consistent with project test tooling.

## Tasks-to-files cross-check

Every file in Files to create or modify appears in a Task: `tools/validate_polyethylene.py` (tasks 3–7), `tools/test_validate_polyethylene.py` (tasks 1, 8), `tools/testdata/pe_thermo/` (task 2).
