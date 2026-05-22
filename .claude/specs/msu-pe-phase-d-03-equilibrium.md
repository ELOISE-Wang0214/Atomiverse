---
title: MSU-PE Phase D.3 — equilibrium-state observable comparison vs LAMMPS
status: code-complete
created: 2026-05-20
parent: msu-pe-validation.md
---

# MSU-PE Phase D.3 — equilibrium-state observable comparison vs LAMMPS

## Summary

Phase D.3 of the MSU amorphous-polyethylene validation closes the equilibrium-state arm of the LAMMPS cross-check: over stage 4 of the protocol (100 K NPT, last 25 000 steps — the equilibrated glass), Atomiverse and LAMMPS must agree on five observables — density ρ, mean potential energy ⟨E_pot⟩, RDF first-peak position r\*, RDF first-peak height g(r\*), and the end-to-end chain distance ⟨R²⟩. This spec extends `tools/validate_polyethylene.py` (the comparison script started by chain link msu-pe-phase-d-02) with an equilibrium-observables section reachable as `validate_polyethylene.py --check equilibrium`, which exits non-zero on any tolerance breach. It is pure analysis tooling: no C++ or production change. The user-visible outcome is a single machine-verifiable gate that confirms the two engines produce the same equilibrated polymer glass.

## Chain dependencies

This is chain link 3 of the `msu-pe-phase-d-*` chain (parent: `msu-pe-validation.md`, Phase D). It has two hard ordering dependencies that are *chain ordering*, not vague defers:

- **msu-pe-phase-d-01-reference-runs** must land first. It produces the reference artifacts this script consumes: the stage-4 thermo log (volume + E_pot per step), the finalized RDF curve, and the trajectory XYZ — for both engines, at both 1 and 8 ranks. This spec does **not** generate any reference data; it only reads what link 01 commits.
- **msu-pe-phase-d-02** creates the `tools/validate_polyethylene.py` file and its per-step thermo section (Phase D.2). This spec **edits the same file** — it adds an equilibrium section and the `--check` dispatch. If link 02 has not landed, this spec's edits have nothing to attach to.

Authoring/landing order is therefore strict: 01 → 02 → 03 (this) → 04. This spec authors only the D.3 equilibrium section; the per-step thermo (D.2) and stress-strain (D.4) checks are other links' jobs and are not (re)specified here.

## Domain basis

Equilibrium observables of an amorphous glassy polymer carry different statistical noise, which sets the per-observable tolerances (from `msu-pe-validation.md` Phase D.3):

- **Density ρ and ⟨E_pot⟩** are well-averaged extensive/intensive thermodynamic quantities; over 25 000 steps they are tight — 1% and 0.5% relative respectively. ρ is derived from the `Thermo` cell volume V and the fixed total mass M as ρ = M / V (report in g/cm³; convert from amu/Å³ via the `core/units.h` factor used elsewhere in the polyethylene benchmarks).
- **RDF first-peak position r\*** is a structural signature — the nearest-neighbour shell of the united-atom CH₂ packing — reproducible to one histogram bin. With the Phase-D RDF binning the bin width is < 0.05 Å, so the tolerance is "within 1 bin".
- **RDF first-peak height g(r\*)** is noisier than its position (it depends on local density fluctuations in the first shell): 3% relative.
- **⟨R²⟩ end-to-end distance** — the squared distance between the first and last UA bead of each chain, averaged over all chains and over the sampled frames — carries large chain-conformational variance (few chains, slow conformational relaxation in a glass), hence the loose 10% relative tolerance.

Block-averaging over the last 25 000 steps of the 100 K NPT stage 4 samples the equilibrated glass, after stage 3's 500→100 K cooling has formed it. No new equations are introduced beyond ρ = M/V and ⟨R²⟩ = ⟨|r_last − r_first|²⟩; both are standard.

## Design

`tools/validate_polyethylene.py` currently runs the thermo and stress-strain checks unconditionally inside `main()`. This spec introduces a `--check` selector so the equilibrium gate is independently invokable and independently exit-coded:

- A new CLI argument `--check {all,thermo,stress-strain,equilibrium}` (default `all`) routes `main()` to the requested section(s). `equilibrium` runs only the new section; `all` preserves current behaviour plus the new section.
- A new dataclass field group on `Tolerances` holds the five D.3 thresholds: `eq_density_rel = 0.01`, `eq_epot_rel = 0.005`, `eq_rdf_peak_pos_abs = 0.05` (Å), `eq_rdf_peak_height_rel = 0.03`, `eq_r2_rel = 0.10`. Frozen-style defaults, overridable via CLI flags.
- A new dataclass `EquilibriumStats` holds the five computed scalars (`density_g_cm3`, `epot_mean_eV`, `rdf_peak_pos_A`, `rdf_peak_height`, `r2_mean_A2`) for one engine.
- A pure function `compute_equilibrium_stats(...)` takes a stage-4 thermo row list, an RDF `(r, g)` table, a trajectory frame list, and per-chain atom-index ranges; returns an `EquilibriumStats`. It reuses the existing `_block_average` / `_to_float` helpers for the thermo-derived quantities.
- Helper `_rdf_first_peak(r, g)` finds the first interior local maximum of g(r) (skipping the r≈0 noise floor) and returns `(r_star, g_star)`.
- Helper `_end_to_end_r2(frames, chain_ranges)` computes ⟨R²⟩ from XYZ frames over stage 4; `chain_ranges` is the list of `(first_atom_index, last_atom_index)` pairs identifying each chain's terminal beads, derived from the MSU topology (10 chains × 1000 UA beads).
- A new XYZ reader `_read_xyz(path)` parses the multi-frame trajectory into a list of frames (each a list of `(x, y, z)`); stdlib-only, consistent with the script's no-numpy constraint.
- A new `_check_equilibrium(av_stats, ref_stats, tol, rep)` compares the two `EquilibriumStats` and appends a `Failure` per breached observable, with `severity = rel_err / threshold`, reusing the existing `Report` / `Failure` plumbing so output formatting and exit-code semantics are unchanged.

Ownership: all new code lives in `validate_polyethylene.py`; no module split is warranted at this size. Input artifact paths follow the link-01 reference layout under `benchmarks/md/reference/polyethylene/` and the Atomiverse run directory passed via `--atomiverse-dir`.

## Files to create or modify

- `tools/validate_polyethylene.py` — add `--check` dispatch, equilibrium tolerances, `EquilibriumStats`, `compute_equilibrium_stats`, RDF/XYZ/⟨R²⟩ helpers, and `_check_equilibrium`.
- `tools/test_validate_polyethylene.py` (new) — pytest unit + integration tests for the equilibrium section, with small synthetic fixtures.

## Tasks

- [x] Write failing tests for RDF first-peak and end-to-end-R² helpers (tools/test_validate_polyethylene.py)
- [x] Implement _read_xyz, _rdf_first_peak, and _end_to_end_r2 in tools/validate_polyethylene.py
- [x] Implement EquilibriumStats and compute_equilibrium_stats in tools/validate_polyethylene.py
- [x] Write failing tests for _check_equilibrium tolerance gating (tools/test_validate_polyethylene.py)
- [x] Implement equilibrium Tolerances fields and _check_equilibrium in tools/validate_polyethylene.py
- [x] Implement --check {all,thermo,stress-strain,equilibrium} dispatch in main()
- [x] Write integration test asserting non-zero exit on a seeded equilibrium breach (tools/test_validate_polyethylene.py)
- [x] Run full check + test suite

## Testing strategy

- **Happy path**: with synthetic stage-4 thermo, RDF, and XYZ fixtures whose Atomiverse and reference values sit inside every tolerance, `--check equilibrium` reports `ALL CHECKS PASSED` and exits 0.
- **Edge cases**: missing reference RDF or trajectory file → graceful note + skip (no crash, consistent with the script's existing optional-reference behaviour); malformed XYZ (atom-count mismatch) → a `Failure`, not an exception; RDF with no interior peak → a `Failure`, not a silent zero; `--check thermo` and `--check stress-strain` still behave exactly as before (no regression of the link-02 sections).
- **Domain validation**: a fixture seeded so density is off by 1.5% (> 1% tol) and ⟨R²⟩ off by 8% (< 10% tol) yields exactly one `Failure` (density), exit code 1, with the breached observable named and its `severity` ≈ 1.5. RDF first-peak-position breach is asserted at the 1-bin (0.05 Å) boundary.

## Out of scope

- Generating the reference trajectory, RDF, or stage-4 thermo — that is `msu-pe-phase-d-01-reference-runs`; this script only consumes those artifacts. (Considered folding a generation fallback in; rejected — it would duplicate link 01 and confound "comparison" with "data production".)
- The per-step thermo comparison (D.2) and the stress-strain curve fits (D.4) — those are other links' sections of the same file; this spec touches only the equilibrium section and the shared `--check` dispatch.
- Any change to C++ production code, the `cpu::observer::RDF` component, or the polyethylene benchmark — the RDF and trajectory are taken as given outputs.
- Plotting/visualization of the equilibrium observables (no `plot_*` companion is in scope here).
