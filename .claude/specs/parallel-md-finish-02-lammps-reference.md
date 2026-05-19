---
title: LAMMPS reference run with self-contained geometry
slug: parallel-md-finish-02-lammps-reference
status: retired
created: 2026-05-19
last-updated: 2026-05-19
parent: parallel-md-engine.md
---

# LAMMPS reference run with self-contained geometry

## Summary

Make `benchmarks/md/reference/polyethylene/lammps.in` executable on this host without an external `msu_pe.data` file by switching it to procedural geometry (`lattice` + `create_atoms` + `create_bonds`) that mirrors `make_pe_box()` in `polyethylene_tension.cpp` — a 4×4×4 chain grid × 18 united-atom carbons. Add a tail block that writes the two reference CSVs (`thermo_reference.csv`, `stress_strain_reference.csv`) via `fix print` directly into this directory. Run `lmp_mpi -in lammps.in` on this host (Homebrew LAMMPS 22 Jul 2025 Update 4), commit both CSVs, and fill the version + hardware blocks in `README.md`. Closes parent ac-005.

## Domain basis

Same Dreiding-style UA polyethylene force field already declared in `lammps.in` (LJ 10.5 Å, harmonic bond k=350 kcal/mol/Å² at r₀=1.529 Å, harmonic angle k=58.35 kcal/mol/rad² at θ₀=112.7°, multi-harmonic dihedral OPLS coefficients). The 5-stage protocol (NVT-500 K → NPT-500 K → cool 500→100 K → NPT-100 K → uniaxial tension at ε̇=10⁻⁶/fs) matches the parent spec. Geometry size is reduced from MSU's 20-chain × 1000-monomer (paper-scale validation) to the 4×4×4 × 18 system the Atomiverse profile binary uses — paper-scale science validation is owned by `msu-pe-validation.md`, not this sub-spec.

## Design

Touched entities:

- `benchmarks/md/reference/polyethylene/lammps.in` — replace `read_data msu_pe.data` with a procedural geometry block; shrink step counts to {1 000, 5 000, 5 000, 5 000, 5 000} so the full pipeline finishes in minutes on this host; add `fix print` blocks writing the two CSVs.
- `benchmarks/md/reference/polyethylene/thermo_reference.csv` — new artifact (written by LAMMPS).
- `benchmarks/md/reference/polyethylene/stress_strain_reference.csv` — new artifact (written by LAMMPS).
- `benchmarks/md/reference/polyethylene/README.md` — fill the `<TBD>` placeholders in the version + hardware blocks with the real LAMMPS tag (`stable_22Jul2025_update4`), compiler version (`clang` from Apple toolchain), hardware fingerprint of this machine, and reproduction command.

No source-tree C++ code changes.

## Files to create or modify

- `benchmarks/md/reference/polyethylene/lammps.in` (modify)
- `benchmarks/md/reference/polyethylene/thermo_reference.csv` (new)
- `benchmarks/md/reference/polyethylene/stress_strain_reference.csv` (new)
- `benchmarks/md/reference/polyethylene/README.md` (modify)

## Tasks

- [x] Replace `read_data msu_pe.data` with procedural geometry (decided to use a `molecule` template `pe_chain.mol` + `lattice` + `create_atoms` rather than `create_bonds many` because LAMMPS does not auto-derive angles/dihedrals from bonds, and the molecule template is the cleanest way to carry the topology). Inline comments in `lammps.in` document the constants matching `make_pe_box()`.
- [x] Reduce `run` step counts to {1000, 5000, 5000, 5000, 5000}; paper-scale schedule deferred to `msu-pe-validation`.
- [x] Add the `fix thermolog all print 100 …` block writing comma-separated `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz` to `thermo_reference.csv` (placed after `thermo` so it captures all 5 stages).
- [x] In stage 5 (tension), add `variable lx0 equal $(lx)` (immediate-evaluation `$()` syntax — needed because `equal lx` re-evaluates each access and gives strain ≡ 0) + `variable strain equal "(lx - v_lx0) / v_lx0"` + `fix stressout all print 100 "${strain},$(pxx)" file stress_strain_reference.csv …`. Also added pre-stage warmup: `minimize` + 200 fs `fix nve/limit 0.1` to dissipate the 0.9 Å chain-end overlap energy before NVT @ 500 K takes over (necessary because LAMMPS' default Verlet step is stricter than Atomiverse's NoseHoover; without the warmup, bonds break at step 5).
- [x] Bump box from 60 Å to 90 Å (lattice spacing 22.5 Å) so the 21.6 Å chain fits within one sub-cell — necessary to avoid chain-on-chain end overlap that crashed the LAMMPS run on the original 60 Å box. Documented in `lammps.in` comments; sub-spec 04 will run Atomiverse at matching `--box-len 90` for the head-to-head comparison.
- [x] Bump stage-5 strain rate from 5e-7/step (= 1e-6/fs) to `erate = 1.5e-4/fs` (150× paper rate) so the 5000-step run reaches ε ≈ 0.37 — needed for ac-003 (final strain ≥ 0.3). Documented as acceptance-scale departure in `lammps.in` comments.
- [x] Run `cd benchmarks/md/reference/polyethylene && /opt/homebrew/bin/lmp_mpi -in lammps.in -log /tmp/lammps_smoke.log` — exits 0 in ~2 s. `thermo_reference.csv` lands with 231 lines, `stress_strain_reference.csv` with 51 lines.
- [x] Spot-checks: final thermo T = 100.75 K (within [80, 150] from stage-4 cooling target 100 K); final stress strain = 0.3723 (≥ 0.3 ✓).
- [x] Update `README.md`: filled `<TBD>` placeholders in the "Required versions" subsection with `lammps-stable_22Jul2025_update4`, Homebrew Cellar path, Open-MPI, AppleClang -O3 -march=native; added `#### Hardware fingerprint` subsection with Apple M4 / 10 cores / macOS 26.4.1; rewrote the "Steps" subsection to drop the moltemplate workflow (no longer needed); refreshed the "What's here" table to mark CSVs / PNGs / pe_chain.mol as committed; added a "Step-count deviation from MSU paper" subsection documenting the {1000,5000,5000,5000,5000} step counts and 150× strain rate.
- [x] Commit all five files (lammps.in, pe_chain.mol, thermo_reference.csv, stress_strain_reference.csv, README.md) in one commit.

## Testing strategy

- **Happy path:** `lmp_mpi -in lammps.in` exits 0 in ≤ 5 minutes on this host. Both CSVs land with the documented header. README has no `<TBD>` markers left in the version block.
- **Edge case:** if procedural `create_bonds` over-bonds at periodic edges, add `bond_atoms` exclusion or shrink the `region` slightly — record in a comment which path was taken.
- **Edge case:** if `lmp_mpi` segfaults at stage 5 due to extreme deformation on the small box, shorten stage 5 step count further (annotate the comment).
- **Domain validation:** `thermo_reference.csv` row at end of stage 1 shows T within 50 K of 500 K (NVT equilibration); row at end of stage 4 shows T within 30 K of 100 K (cooling worked). `stress_strain_reference.csv` shows σ_xx monotonically rising in the first 20% of rows (elastic regime) — this proves the engine actually deformed the box.

## Out of scope

- 20-chain × 1000-monomer paper-scale system (owned by `msu-pe-validation`).
- Tight force-field parameter agreement with `force_field.yaml` (the parameters in `lammps.in` are already aligned at spec creation; no re-tuning).
- Atomiverse-vs-LAMMPS numerical comparison of thermo / stress (covered by sub-spec 04).
- Cross-platform `lmp_mpi` (this run is recorded on Apple Silicon / Homebrew only).
