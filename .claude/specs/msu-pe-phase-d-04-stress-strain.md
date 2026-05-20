---
title: MSU-PE Phase D.4 — stress-strain curve comparison and Phase D exit gate
status: approved
created: 2026-05-20
parent: msu-pe-validation.md
---

# MSU-PE Phase D.4 — stress-strain curve comparison and Phase D exit gate

## Summary

This is chain link 4 of 4 (final) of the `msu-pe-phase-d-*` chain under `.claude/specs/msu-pe-validation.md`. It gives `tools/validate_polyethylene.py` a first-class stress-strain check (`--check stress-strain`) that overlays the Atomiverse σ_xx(ε) uniaxial-tension curve on the LAMMPS reference curve from chain link `msu-pe-phase-d-01-reference-runs` (stage 5), extracts three mechanical descriptors — elastic modulus, yield stress, strain-hardening slope — and exits non-zero if any of the three breaches ±5% relative tolerance. It also adds a `--check all` entry point that runs the D.2, D.3, and D.4 checks together at a given rank count, so an operator can verify the complete Phase D exit gate with one invocation at 1 rank and at 8 ranks. When this link lands and both rank counts pass, Phase D of the parent `msu-pe-validation.md` is complete; the operator then updates the Phase D status row of `msu-pe-validation.md` to reflect that. This is a pure-analysis change in `tools/` — no C++ or production-code change.

## Domain basis

Uniaxial tension of a glassy amorphous polymer produces a σ_xx(ε) stress-strain curve whose shape encodes three standard engine-validation descriptors:

- **Elastic modulus** `E = dσ/dε` — small-strain linear slope, fitted by least squares over ε < 0.02 (Hookean regime before any plastic event).
- **Yield stress** `σ_y = max σ_xx` over ε ∈ [0, 0.1] — the stress peak before plastic flow sets in.
- **Strain-hardening slope** `dσ/dε` — large-strain linear slope, fitted by least squares over ε ∈ [0.5, 1.0], where entangled chains align and the curve re-stiffens.

Units: strain ε is dimensionless engineering strain (`box_x / box_x0 − 1`); stress σ_xx is reported in GPa (`sigma_xx_GPa` column; conversion `1 eV/Å³ = 160.21766 GPa` applied in the benchmark). `E` and the hardening slope therefore carry units of GPa (slope of GPa-vs-dimensionless); `σ_y` is in GPa.

Tolerance: ±5% relative on each of the three descriptors. This is the standard engine-validation tolerance and matches the primary success criterion of `msu-pe-validation.md` ("elastic modulus, yield stress, and hardening slope agree with LAMMPS within ±5% relative error"). The σ_xx(ε) curve is the engine's headline mechanical prediction; the tolerances are taken verbatim from `msu-pe-validation.md` Phase D.4. References: MSU CAVS ICME amorphous-PE tension protocol (<https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>); Tschopp 2010 Dreiding-UA parameter attribution recorded in `benchmarks/md/reference/polyethylene/SOURCES.md`.

## Design

`tools/validate_polyethylene.py` already contains the descriptor-extraction core: `_compute_stress_strain_stats` (the three least-squares / max fits) and `_check_stress_strain` (the ±5% comparison and failure reporting), plus `StressStrainStats`, `_linear_fit`, and `Tolerances` with the `modulus_rel / yield_rel / hardening_rel = 5e-2` fields. This sub-spec does **not** rewrite that physics; it wires it into an explicit, machine-driven CLI surface and adds the Phase D gate.

New / changed entities in `tools/validate_polyethylene.py`:

- A `--check {thermo,equilibrium,stress-strain,all}` argument (default `all`) selecting which Phase D check(s) run. `thermo` → D.2, `equilibrium` → D.3, `stress-strain` → D.4, `all` → D.2+D.3+D.4 together.
- A `run_checks(check, av_dir, ref_dir, tol, thermo_rows, rep)` dispatch function that calls the existing `_check_thermo` / `_check_stress_strain` (and the D.3 equilibrium check from link 03) according to the selected `--check` value, so `main` becomes a thin argument-parse-then-dispatch shell.
- An optional `--ranks {1,8}` argument recorded into the `Report` notes and echoed in output, so an operator can label which rank count a run covers. It does not change the analysis (the CSVs already encode a single run); it makes the "pass at 1 rank and 8 ranks" gate auditable from the printed report.
- Exit-code contract preserved and made explicit per check: `0` all selected checks pass, `1` any tolerance breach, `2` required input missing. `--check stress-strain` with no `stress_strain_reference.csv` present is exit `2` (a missing reference is an operator error for this gate, not a silent skip) — distinct from the current "no reference → note and skip" behavior used by the broader `all` run when only some references exist.

Lifecycle: the script is stateless and side-effect-free apart from stdout/stderr and its exit code. It consumes `thermo.csv` / `stress_strain.csv` produced by `benchmarks/md/polyethylene_tension` (already emitting `step,stage,strain,sigma_xx_eV_per_A3,sigma_xx_GPa`) and the `*_reference.csv` files committed by chain link 01. No new files are written.

Chain dependency (stated plainly): this link depends on chain links `msu-pe-phase-d-01-reference-runs` (commits `stress_strain_reference.csv` and `thermo_reference.csv` for 1 rank and 8 ranks), `msu-pe-phase-d-02` (creates `validate_polyethylene.py` with the D.2 thermo check), and `msu-pe-phase-d-03` (adds the D.3 equilibrium check). The `--check all` dispatch requires the D.2 and D.3 functions to exist. If links 02/03 have not yet landed when this is implemented, `run_checks` for `all` must fail loudly (clear error naming the missing check function) rather than silently degrade — but authoring the D.2/D.3 checks themselves is **out of scope** for this link.

## Files to create or modify

- `tools/validate_polyethylene.py` — add `--check` and `--ranks` arguments, `run_checks` dispatch, stress-strain-specific missing-reference exit code; reuse existing `_compute_stress_strain_stats` / `_check_stress_strain`.
- `tools/test_validate_polyethylene.py` (new) — pytest coverage for the stress-strain descriptor extraction, the `--check` dispatch, and the exit-code contract, using synthetic CSV fixtures.
- `.claude/specs/msu-pe-validation.md` — update the Phase D status row to reflect Phase D complete once D.2/D.3/D.4 pass at 1 and 8 ranks (operator action on landing).

## Tasks

- [ ] Write failing tests for stress-strain descriptor extraction and `--check` dispatch (tools/test_validate_polyethylene.py)
- [ ] Implement `--check {thermo,equilibrium,stress-strain,all}` and `--ranks` arguments in tools/validate_polyethylene.py
- [ ] Implement `run_checks` dispatch routing each `--check` value to its D.2/D.3/D.4 check function in tools/validate_polyethylene.py
- [ ] Implement stress-strain-specific exit code 2 when `stress_strain_reference.csv` is absent in tools/validate_polyethylene.py
- [ ] Add docstrings with units for the new `--check`/`--ranks` surface and `run_checks` per project doc style
- [ ] Verify `--check stress-strain` flags a >5% modulus/yield/hardening deviation on a synthetic breach fixture and exits non-zero
- [ ] Update .claude/specs/msu-pe-validation.md Phase D status row to "Phase D complete" reflecting D.2/D.3/D.4 pass at 1 and 8 ranks
- [ ] Run full check + test suite

## Testing strategy

- **Happy path**: synthetic `stress_strain.csv` + `stress_strain_reference.csv` whose three descriptors agree within ±5% → `--check stress-strain` exits `0`; `--check all` with all three reference files present exits `0`.
- **Edge cases**:
  - `stress_strain_reference.csv` missing → `--check stress-strain` exits `2` (required-input error), distinct from the `all`-run skip-with-note path.
  - Atomiverse `stress_strain.csv` missing or with mismatched `strain` / `sigma_xx_GPa` column lengths → exit `1` with a malformed-input failure.
  - Degenerate strain window (fewer than 2 points with ε < 0.02) → `_linear_fit` returns slope 0 without raising; covered so a sparse curve does not crash the gate.
  - Unknown `--check` value → argparse rejects with exit `2`.
- **Domain validation**: synthetic curves with a deliberately injected +6% elastic-modulus offset, a +6% yield-stress offset, and a +6% hardening-slope offset each independently produce a `stress_strain:<descriptor>` failure with `severity > 1.0`; a +4% offset on each passes. Confirms the ±5% boundary is enforced per descriptor, matching `msu-pe-validation.md` Phase D.4.

## Out of scope

- Authoring the D.2 per-step thermo check or the D.3 equilibrium-observables check — those are owned by chain links `msu-pe-phase-d-02` and `msu-pe-phase-d-03` respectively; this link only dispatches to them.
- Producing the LAMMPS reference CSVs (`stress_strain_reference.csv`, `thermo_reference.csv`) — owned by chain link `msu-pe-phase-d-01-reference-runs`.
- Running the actual 1-rank and 8-rank LAMMPS and Atomiverse simulations — an operator action; this link delivers the comparison tooling only.
- Any change to `benchmarks/md/polyethylene_tension.cpp` or other C++/production code — the benchmark already emits the required `stress_strain.csv` schema.
- Plotting / curve overlay images — the "overlay" requirement is satisfied numerically by the three-descriptor comparison; visual plotting is `tools/plot_*.py` territory and not gated.
- Phase E performance benchmarks and the Phase E exit gate.
