---
title: Migrate Atomiverse PE benchmark + tests to Dreiding-UA force field (MSU paper)
slug: pe-ff-migrate-dreiding-ua
status: approved
created: 2026-05-19
parent: msu-pe-validation.md
---

# Migrate Atomiverse PE force field to Dreiding-UA (MSU paper)

## Summary

Replace the OPLS-AA-like angle (θ₀ = 112.7°, k = 58.35 kcal/mol/rad²) and OPLS-style dihedral (K₁..K₄ = 2.0058, −2.5128, 0.2710, 6.2900) constants currently used in `benchmarks/md/polyethylene_tension.cpp`, `benchmarks/md/pe_profile_single.cpp`, and 11 `tests/{cpu,mpi}/*_pe_*.cpp` files with the **Dreiding united-atom (Tschopp 2010)** parameters published in the MSU CAVS ICME `in.deform_dreiding.txt` reference deck: θ₀ = 109.5°, k = 60 kcal/mol/rad², K₁..K₄ = 1.73, −4.49, 0.776, 6.99. Regenerate `benchmarks/md/reference/polyethylene/pe_chain.mol` initial positions from the new half-apex angle (35.25° vs the old 33.65°). Re-run the small-scale `lammps.in` self-contained reference and re-commit `thermo_reference.csv` + `stress_strain_reference.csv`. Update `SOURCES.md` discrepancy table to "resolved". Outcome: Atomiverse and the MSU paper LAMMPS deck `lammps_msu_paper.in` now exercise the same physics, **unblocking `msu-pe-validation.md` Phase D** (the ±5% LAMMPS stress-strain comparison).

## Domain basis

The Dreiding force field for united-atom polyethylene (Tschopp 2010, Hossain–Tschopp–Horstemeyer 2010 *Polymer*) treats each CH₂ as a single bead with sp³-tetrahedral covalent geometry:

- **Bond** (harmonic): `V = k_b · (r − r₀)²` with `k_b = 350 kcal/mol/Å²`, `r₀ = 1.53 Å`. (Atomiverse uses 1.529 Å — sub-Å rounding, no change needed; both round to the same eV/Å² when converted via `kcal/mol → eV` at FP precision.)
- **Angle** (harmonic): `V = k_θ · (θ − θ₀)²` with `k_θ = 60 kcal/mol/rad²`, `θ₀ = 109.5°`. Sp³ tetrahedral. Atomiverse currently uses 58.35 / 112.7° — a 3.2° offset that compounds across the 17 dihedrals per chain.
- **Dihedral** (LAMMPS `multi/harmonic`, Σ Kₙ · cosⁿ⁻¹(φ)): `K₁..K₄ = 1.73, −4.49, 0.776, 6.99`. Atomiverse currently uses an OPLS-style set `2.0058, −2.5128, 0.2710, 6.2900` from a different paper family; the value sets are mutually inconsistent and produce different gauche-trans populations.
- **Pair** (LJ 12-6 with 10.5 Å cutoff, 1–2 / 1–3 / 1–4 exclusion `special_bonds lj 0 0 0`): `ε = 0.112 kcal/mol`, `σ = 4.01 Å`. (Already match — no change.)

All values verbatim from `benchmarks/md/reference/polyethylene/SOURCES.md` § Force field, sourced from the upstream `in.deform_dreiding.txt`.

Initial chain geometry change: the half-apex angle `(π − θ₀) / 2` shifts from 33.65° to 35.25°. With `r₀ = 1.529` Å, the zig-zag spacing changes from `(dy, dz) = (1.272798, 0.847247)` Å to `(1.248645, 0.882455)` Å — a 1.9% contraction along the chain axis and a 4.2% widening of the transverse offset. This affects the initial coordinate generation in `make_pe_box()` and the position table in `pe_chain.mol`, but does **not** affect the eventual equilibrated structure (the NVT-500 K stage in the protocol fully randomises the chain).

## Design

The 15 affected files duplicate the same five constants across each file (a recognised smell — see "Out of scope" below). This spec does literal value replacement, not refactoring. Touched files:

| File | What changes |
|---|---|
| `benchmarks/md/polyethylene_tension.cpp` | `kAngleK_kcal` 58.35 → 60, `kAngleTheta0` rad-encoding `112.7` → `109.5`, dihedral coefficient array. |
| `benchmarks/md/pe_profile_single.cpp` | Same three constants — file is a near-duplicate of the polyethylene_tension constants block. |
| `tests/cpu/md/test_pe_single_chain.cpp` `test_pe_multi_chain.cpp` `test_pe_nve.cpp` `test_pe_nvt.cpp` `test_lj_exclusions.cpp` | Same three constants. The tests use formula-based oracles (`V = k · Δθ²`, NVE drift bounds, T equipartition) — no numerical reference re-calibration is needed beyond the constant swap. |
| `tests/mpi/test_pe_*.cpp` (8 files) `test_nose_hoover_t_ramp.cpp` | Same three constants. Multi-rank correctness tests compare 1-rank vs N-rank, so the absolute values irrelevant; the constants just need to be consistent. |
| `benchmarks/md/reference/polyethylene/pe_chain.mol` | 18 atom-position lines regenerated for new `(dy, dz)`. The bond/angle/dihedral connectivity tables are unchanged (topology, not geometry). |
| `benchmarks/md/reference/polyethylene/lammps.in` | `angle_coeff 1 58.35 112.7` → `60 109.5`; `dihedral_coeff 1 2.0058 -2.5128 0.2710 6.2900 0.0` → `1.73 -4.49 0.776 6.99 0.0`. |
| `benchmarks/md/reference/polyethylene/thermo_reference.csv` | Re-run LAMMPS on the updated `lammps.in`, re-commit. |
| `benchmarks/md/reference/polyethylene/stress_strain_reference.csv` | Same. |
| `benchmarks/md/reference/polyethylene/SOURCES.md` | "Known discrepancies" section: change "DOES NOT match" text + table to "matches as of `<this-commit-SHA>`"; preserve historical context as a one-line footnote. |
| `.claude/specs/msu-pe-validation.md` | Status snapshot row for Phase D: `🟡 blocked on FF migration` → `✅ unblocked` (or `🟡 ready to run` if paper-scale runs haven't been driven yet). |

No new header files; no API surface changes; no new components; no architectural impact. Component lifecycle, MPI layout, Driver taxonomy, and benchmark CLI surfaces are all untouched. The Box / Frame / Driver / Plan invariants from `CLAUDE.md` are unaffected.

### Implementation notes

1. **Single-pass replacement.** All 15 source files use exactly the same three constant names and the same literal numeric values. A `sed`-style targeted replacement of the three literal numbers (`58.35` → `60`, `112.7` → `109.5`, and the four dihedral coefficients) is mechanically safe — these literals do not appear elsewhere in the codebase by accident (grep verified during research). Hand-edit `pe_chain.mol` and `lammps.in` separately because they live outside the C++ source pattern.
2. **`pe_chain.mol` regeneration.** Compute the new positions in Python (3 lines: `dy = 1.529 * cos(35.25°)`, `dz = 1.529 * sin(35.25°)`, then `atom_i = (0, i*dy, (i%2)*dz)`). Replace the Coords block. Update the docstring constants at the top of the file.
3. **`lammps.in` re-run.** After editing the coeff lines, re-run `cd benchmarks/md/reference/polyethylene && /opt/homebrew/bin/lmp_mpi -in lammps.in -log /tmp/lammps_smoke.log` and capture the new `thermo_reference.csv` and `stress_strain_reference.csv` (these are written by `fix print` blocks already in the deck). Both should still satisfy parent ac-002 / ac-003 from `parallel-md-finish-02` (≥10 rows, final T ∈ [80, 150] K, final strain ≥ 0.3). If LAMMPS now fails (e.g. the geometry minimiser doesn't converge with the new params on the small 4×4×4×18 system), document the failure and either tighten the minimisation or rebalance the warmup — do not paper over a regression.
4. **Test recalibration cost.** The Atomiverse tests use formula-based oracles (`V_angle = k_θ · (θ − θ₀)²` for a small displacement, NVE drift `< 1e-4`, NVT temperature within `±10 K` of setpoint). None of them embed a numerical reference value that would need updating to a new LAMMPS run. The expected outcome is **every existing PE test still passes** with the new constants; if a test fails, the migration has introduced a bug, not exposed a stale reference.

## Files to create or modify

Modified only — no new files:

```
benchmarks/md/polyethylene_tension.cpp
benchmarks/md/pe_profile_single.cpp
benchmarks/md/reference/polyethylene/lammps.in
benchmarks/md/reference/polyethylene/pe_chain.mol
benchmarks/md/reference/polyethylene/thermo_reference.csv      (regenerated)
benchmarks/md/reference/polyethylene/stress_strain_reference.csv (regenerated)
benchmarks/md/reference/polyethylene/SOURCES.md
tests/cpu/md/test_lj_exclusions.cpp
tests/cpu/md/test_pe_multi_chain.cpp
tests/cpu/md/test_pe_nve.cpp
tests/cpu/md/test_pe_nvt.cpp
tests/cpu/md/test_pe_single_chain.cpp
tests/mpi/test_nose_hoover_t_ramp.cpp
tests/mpi/test_pe_chain_crossing.cpp
tests/mpi/test_pe_chain_drifting.cpp
tests/mpi/test_pe_chains_8rank.cpp
tests/mpi/test_pe_chains_multirank.cpp
tests/mpi/test_pe_full_multirank.cpp
tests/mpi/test_pe_msu_pipeline.cpp
tests/mpi/test_pe_nvt_multirank.cpp
.claude/specs/msu-pe-validation.md                              (status table only)
```

## Tasks

- [ ] Write failing test `tests/cpu/md/test_pe_ff_constants.cpp` that pin-points the Dreiding-UA values (asserts `kAngleK_kcal == 60.0`, `kAngleTheta0` corresponds to 109.5°, `kDihedralK[0..3] == {1.73, -4.49, 0.776, 6.99}` in the included header). Confirm RED.
- [ ] Replace the three constant values in `benchmarks/md/polyethylene_tension.cpp`. Confirm RED test now passes; rebuild + run unit-level PE tests.
- [ ] Replace the three constant values in `benchmarks/md/pe_profile_single.cpp` (near-duplicate of above).
- [ ] Replace the three constant values in each of the 13 `tests/{cpu,mpi}/test_pe_*.cpp` + `test_nose_hoover_t_ramp.cpp` + `test_lj_exclusions.cpp` files. Confirm full `ctest -L "cpu"` green.
- [ ] Confirm full `ctest -L "mpi"` green at 1/2/8 ranks (multi-rank PE pipeline tests are sensitive to all three constants).
- [ ] Regenerate `pe_chain.mol` atom-position table using `dy = 1.529 · cos(35.25°)`, `dz = 1.529 · sin(35.25°)`; update the docstring constants header.
- [ ] Update `benchmarks/md/reference/polyethylene/lammps.in` `angle_coeff` and `dihedral_coeff` lines to MSU values.
- [ ] Re-run `/opt/homebrew/bin/lmp_mpi -in lammps.in` and re-commit the regenerated `thermo_reference.csv` + `stress_strain_reference.csv`; verify both still ≥ 11 lines with right headers and final T ∈ [80, 150] K + final strain ≥ 0.3.
- [ ] Update `SOURCES.md` "Known discrepancies" subsection: mark resolved, preserve the original mismatch table as a historical footnote, link to this spec's commit SHA.
- [ ] Update `msu-pe-validation.md` status snapshot: flip Phase D from `🟡 blocked on FF migration` to `🟡 ready to run` (Phase D paper-scale comparison is now unblocked but not yet executed — that's a downstream task, not part of this spec).
- [ ] Run full `ctest --test-dir build --output-on-failure` (no label filter) and verify no regression in the 87-test baseline from `parallel-md-finish-06-close-out`.

### Architecture audit task (uses this migration as a probe)

The mechanical part of this spec (constant value replacement across 15 files) is itself a signal: if the codebase had a LAMMPS-style abstraction (`bond_style harmonic` + a single parameter registration site), this migration would touch **one** file, not fifteen. Use this opportunity to audit the gap.

- [ ] **Audit potential abstraction across CPU + CUDA + torch backends.** Delegate to the `av-architect` agent (read-only). Compare:
  1. **Class hierarchy**: does either `cpu::potentials::HarmonicBond` or `cuda::HarmonicBond` (or `torch::*`) inherit from a common `Potential` / `PairStyle` / `BondStyle` base class? If not, what's the closest common interface (currently `ComponentBase`)?
  2. **Cross-backend naming**: is `HarmonicBond` namespaced consistently (`ATV::cpu::HarmonicBond` vs `ATV::cuda::HarmonicBond`) or is there a name collision in the bare `ATV` namespace?
  3. **Capability/role coverage**: does the CUDA-side class declare `role()` and `capabilities()` like the CPU one (Phase 5 added these to CPU only)?
  4. **Parameter registration**: does a user `set_pair`-style helper exist that dispatches to the right backend, or must the user instantiate the backend-specific class manually?
  5. **Plugin / runtime selection**: is there any string-keyed registry that maps `"harmonic"` → factory? Or is the FF choice baked into the user's call site at compile time?
  6. **Adding a new FF form (e.g. Morse bond)**: enumerate the touch-points. How many files? How many tests need recalibration? Same scaling story as this Dreiding-UA migration?
  
  Deliverable: `.claude/notes/potential-abstraction-audit.md` (passive note, project-scope) with a per-question table of findings, gap severity (low/medium/high), and one concrete proposal per gap. Read-only audit — no source changes in this task.

- [ ] **Promote the audit findings into a follow-up spec proposal** at the bottom of the audit note: name the next spec (e.g. `pe-ff-constants-extract` or `potential-base-class-refactor`), list the smallest concrete refactor that would close the most severe gap, and reference this spec's commit SHA.

## Testing strategy

**Happy path:**

- Every existing PE test passes unchanged (formula-based oracles adapt to new constants).
- The new `test_pe_ff_constants.cpp` pins the Dreiding-UA values so future drift back to OPLS triggers a RED test.
- LAMMPS smoke run of the updated small-scale `lammps.in` completes in seconds.
- Full ctest count = pre-migration count (87) + 1 (the new pin test) = 88.

**Edge cases:**

- `test_pe_nve.cpp`'s NVE drift bound (`< 1e-4 relative`) was tuned on the OPLS-like params; a 3.2° shift in θ₀ stiffens the angle oscillation modestly. If the drift bound is now violated, the tolerance was over-tight, not the test; document a one-line tolerance bump with rationale.
- `test_pe_msu_pipeline_n8.cpp` runs the full 5-stage protocol at 8 ranks and is the most sensitive integration test. Failure here = something deeper than a constant swap and aborts the spec.
- `lammps.in`'s `minimize` + `fix nve/limit 0.1` warmup (added by `parallel-md-finish-02` for the 90 Å box) may need re-tuning if the new angle/dihedral mix overshoots harder at start.

**Domain validation** (`$META.science.required: true`):

- After migration, a fresh single-chain (C₂₀H₄₂, all-trans) energy decomposition under Atomiverse should agree with LAMMPS under `lammps_msu_paper.in` constants on the same geometry to **1e-6 eV absolute** per term. This is the test that fundamentally demonstrates "we now run the MSU paper's FF". (Spec body says this is owed to Phase D in `msu-pe-validation`; this spec only verifies the consistency check.)
- Single-chain NVT @ 500 K block-averaged ⟨T⟩, ⟨E_pot⟩, ⟨P⟩ should agree with LAMMPS to 0.5% relative over 1000 steps. (Same — owed to Phase D, this spec just enables it.)

**Performance:**

- No expected change in any kernel timing — same number of bonds/angles/dihedrals per step, only the coefficient values differ. `profile_1.csv` and `profile_8.csv` `total_ms` per component should match the pre-migration numbers to ≤ 1% (any larger drift is suspect and worth investigating).

## Out of scope

- **Refactor the duplicated FF constants into a shared header.** The 15 files all declare `kAngleK_kcal`, `kAngleTheta0`, etc. in their own scope. A single `include/atomiverse/benchmarks/pe_ff_constants.h` would prevent the next FF migration from being this painful. Tracked as a follow-up `pe-ff-constants-extract` spec; not done here to keep this spec mechanically simple and reviewable.
- **Acting on the architecture-audit findings.** This spec only *captures* the audit (Tasks "Audit potential abstraction" + "Promote the audit findings"). Introducing a `Potential` / `PairStyle` / `BondStyle` base class, a string-keyed registry, or a cross-backend `set_pair("harmonic", ...)` dispatch is the **follow-up spec the audit proposes**, not work this spec performs. The audit is the deliverable; the refactor is downstream.
- **Paper-scale Atomiverse run on `msu_pe.data`** (10 chains × 1000 monomers). That's `msu-pe-validation.md` Phase D's actual ±5% comparison; this spec only unblocks it.
- **Tools to read `msu_pe.data` directly into a `cpu::Frame`** (via `molrs` LAMMPS data reader). The MSU data file uses LAMMPS's "molecular" atom style with explicit Bonds/Angles/Dihedrals sections; the `molrs` reader needs to be confirmed compatible or extended. Tracked separately.
- **Re-tuning `parallel-md-finish-02`'s small-scale `lammps.in` step counts.** Even after the FF migration, the {1000, 5000, 5000, 5000, 5000} step schedule remains the right *acceptance-scale* size; paper-scale {10000, 50000, 50000, 50000, 200000} stays a separate Phase D activity.
- **Updating `parallel-md-finish-04`'s `scaling_report.md` LAMMPS-ratio numbers.** The kernel-speedup rows (LJ::launch, HarmonicBond::launch) and the single-core / 8-rank LAMMPS ratios should not move meaningfully; if they do, that's a separate finding to investigate, not a deliverable here.
