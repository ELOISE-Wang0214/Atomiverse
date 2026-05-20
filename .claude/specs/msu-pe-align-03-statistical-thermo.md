---
title: MSU PE align 03 — statistical block-average thermo parity (supersedes D.2)
status: code-complete
created: 2026-05-20
parent: msu-pe-validation.md
---

# MSU PE align 03 — statistical block-average thermo parity (supersedes D.2)

## Summary

This sub-spec supersedes `msu-pe-phase-d-02-thermo-parity`. The superseded spec compared the LAMMPS and Atomiverse thermo CSVs **per step** with tight 1e-4 / 1e-3 relative tolerances over the first 1000 steps. That contract is physically wrong for this protocol: stage 1 of the MSU reference is a stochastic Langevin thermostat (`fix langevin ... 904297`, an RNG-seeded run), so the two engines' trajectories are uncorrelated from step one and a per-step comparison can never pass. This sub-spec replaces the per-step logic in `tools/validate_polyethylene.py` with a **statistical block-average comparison**: each thermo column is reduced to block averages, and the two engines are required to agree on the *ensemble means* of `⟨T⟩` and `⟨E_pot⟩` to within ~2%, which is the physically meaningful parity statement for a stochastic equilibration. Invoked as `validate_polyethylene.py --check thermo`, it exits non-zero on any breach and prints the breaching column, the two means and the relative gap. No C++ or production-source change; the superseded spec body and acceptance file are rewritten in place.

## Domain basis

Two MD engines running the *same* stochastic-thermostat protocol with different RNG streams do not share a trajectory — they sample the *same statistical ensemble*. The correct cross-engine parity statement is therefore agreement of ensemble-averaged observables, not of instantaneous values. Block averaging (Flyvbjerg & Petersen 1989, *J. Chem. Phys.* 91, 461, doi:10.1063/1.457480) partitions a time series into contiguous blocks, takes the per-block mean, and treats the block means as the sampling units — this both yields the ensemble mean and gives a variance estimate that accounts for serial correlation in the raw series. For an equilibrated NVT/NPT segment the mean temperature is fixed by the thermostat setpoint and the mean potential energy by the force field, so two correct engines must agree on `⟨T⟩` and `⟨E_pot⟩` to within the combined statistical error of the two block-averaged estimates. The ~2% relative tolerance on `⟨T⟩` and `⟨E_pot⟩` is the contract bar — loose enough to absorb finite-sample block-mean scatter, tight enough that a genuine force-field or thermostat-setpoint discrepancy fails. Units: T in K, energies in eV, consistent with project convention; the comparison is relative so unit-agnostic provided both CSVs use the same unit per column (a constraint owned by the reference-run link).

## Design

`tools/validate_polyethylene.py` currently carries the per-step thermo logic from the superseded D.2 spec (`ThermoTolerances`, `StepDeviation`, `compare_thermo_per_step`, `run_thermo_check`) plus stress-strain code. This sub-spec **replaces the thermo-comparison logic** with statistical block-average comparison; the `--check` mode selector and the stress-strain code are retained. Entities:

- `BlockStat` — frozen dataclass `{column, mean, sem, n_blocks}`: the block-averaged mean of one thermo column, its standard error across blocks, and the block count.
- `StatTolerances` — frozen dataclass holding the per-column relative tolerances on the means; defaults are ~2% on `⟨T⟩` and `⟨E_pot⟩` (the columns that carry a physically fixed ensemble mean). Replaces `ThermoTolerances`.
- `block_average(values, block_size)` — pure function: partitions a numeric series into contiguous blocks of `block_size`, returns a `BlockStat` (mean of block means, standard error across blocks, block count). A trailing partial block shorter than `block_size` is dropped so every block has equal weight.
- `compare_thermo_statistical(av_rows, ref_rows, tol, block_size)` — pure function: block-averages each compared column for both engines, computes the relative gap between the two means, returns the list of breaches (one entry per column whose relative mean gap exceeds its `StatTolerances` threshold). Empty list ⇒ pass.
- `run_thermo_check(pairs, tol, block_size)` — orchestrates the comparison over `(label, av_csv, ref_csv)` pairs (`"1rank"` / `"8rank"`); aggregates breaches into a `Report`.
- `main` — keeps the `--check {thermo,all}` argument; `--check thermo` resolves the reference CSVs under `benchmarks/md/reference/polyethylene/` and runs `run_thermo_check`. Stress-strain code stays under `--check all`.

The comparison no longer aligns engines step-by-step — it compares reduced statistics, so a differing logging cadence between engines is tolerated as long as each series has enough samples to form blocks. Missing input CSVs exit code 2; tolerance breach exits code 1; clean exit 0.

## Files to create or modify

- `tools/validate_polyethylene.py` — modify: remove `ThermoTolerances`, `StepDeviation`, `compare_thermo_per_step`; add `BlockStat`, `StatTolerances`, `block_average`, `compare_thermo_statistical`; rewrite `run_thermo_check` for statistical comparison; update `main`'s thermo branch.
- `tools/test_validate_polyethylene.py` — modify: replace per-step unit tests with tests for `block_average` (known-series mean/SEM, partial-block drop) and `compare_thermo_statistical` (within-tolerance pass, mean-gap breach, per-column threshold separation); keep the `main(["--check","thermo", ...])` integration test, retargeted at statistical fixtures.
- `tools/testdata/pe_thermo/` — modify: replace the per-step matching/breaching fixtures with statistical fixtures — a pair whose block-averaged means agree within tolerance and a pair whose `⟨T⟩` mean gap exceeds it, for both `1rank` and `8rank` labels.
- `.claude/specs/msu-pe-phase-d-02-thermo-parity.md` — rewrite: replace the per-step spec body with the statistical-comparison spec (this superseding content).
- `.claude/specs/msu-pe-phase-d-02-thermo-parity.acceptance.md` — rewrite: regenerate the acceptance contract for the statistical comparison, criteria restarting at `ac-001`, all `status: pending`.

## Tasks

- [x] Add `--check {thermo,all}` mode and reference-CSV path resolution to `main` in `tools/validate_polyethylene.py`
- [x] Write failing tests for `block_average` (known mean/SEM, partial-block drop) (`tools/test_validate_polyethylene.py`) (rework: per-step test suite replaced by block-average tests)
- [x] Write failing tests for `compare_thermo_statistical` (mean-gap pass/breach, per-column thresholds) (`tools/test_validate_polyethylene.py`) (rework: replaces `compare_thermo_per_step` tests)
- [x] Replace per-step thermo fixtures with statistical agree/breach fixtures in `tools/testdata/pe_thermo/` (rework: per-step CSV fixtures invalidated by the statistical contract)
- [x] Implement `BlockStat` and `StatTolerances` dataclasses in `tools/validate_polyethylene.py` (rework: replaces `ThermoTolerances`/`StepDeviation`)
- [x] Implement `block_average` and `compare_thermo_statistical` in `tools/validate_polyethylene.py` (rework: replaces `compare_thermo_per_step`)
- [x] Rewrite `run_thermo_check` for statistical comparison and update `main`'s thermo branch in `tools/validate_polyethylene.py` (rework: per-step orchestration replaced)
- [x] Rewrite `.claude/specs/msu-pe-phase-d-02-thermo-parity.{md,acceptance.md}` as the statistical-comparison spec
- [x] Run full check + test suite

## Testing strategy

- **Happy path** — `block_average` on a known constant series returns that constant as the mean and zero SEM; `compare_thermo_statistical` on a CSV pair whose block-averaged `⟨T⟩` and `⟨E_pot⟩` agree within tolerance returns an empty breach list; `main(["--check","thermo"])` on agreeing fixtures exits 0.
- **Breach detection** — a fixture pair whose `⟨T⟩` block-mean gap exceeds the ~2% threshold yields a single breach naming `T`, the two means and the relative gap; `main` exits 1.
- **Per-column threshold** — a mean gap that breaches `⟨T⟩` but stays under the `⟨E_pot⟩` tolerance fails on `T` only, confirming per-column tolerances are applied independently.
- **Partial-block drop** — `block_average` on a series whose length is not a multiple of `block_size` drops the trailing partial block, so the block count is `floor(len/block_size)`.
- **Missing input** — `main(["--check","thermo"])` with absent reference CSVs exits 2.
- **Domain validation** — once `msu-pe-phase-d-01-reference-runs` commits the real reference CSVs, `validate_polyethylene.py --check thermo` confirms Atomiverse and LAMMPS agree on the block-averaged `⟨T⟩` and `⟨E_pot⟩` within ~2% for both the 1-rank and 8-rank pairs. This run is **blocked until link 01 commits its reference CSVs** — a hard chain dependency. The unit and integration tests above use synthetic fixtures and are not blocked.

## Out of scope

- The reference runs themselves (LAMMPS and Atomiverse thermo CSV generation) — owned by `msu-pe-phase-d-01-reference-runs`.
- Equilibrium-state observables: density, RDF peak, end-to-end distance — Phase D.3.
- Stress-strain curve analysis — Phase D.4; the existing stress-strain code is left untouched under `--check all`.
- Plotting / visualization of the thermo traces.
- Any non-stdlib dependency in the script; it stays stdlib-only (`csv`, `argparse`, `statistics`). pytest is used only by the test module.
- Autocorrelation-time estimation / optimal block-size selection — a fixed `block_size` argument is sufficient for this parity gate; adaptive block sizing is a possible later refinement.
