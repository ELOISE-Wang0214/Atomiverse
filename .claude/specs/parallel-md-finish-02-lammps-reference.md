---
title: LAMMPS reference run with self-contained geometry
slug: parallel-md-finish-02-lammps-reference
status: approved
created: 2026-05-19
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

- [ ] Replace the `read_data msu_pe.data` line in `lammps.in` with a procedural geometry block: `lattice sc 4.7`, `region` + `create_box`, `create_atoms` for 4×4×4 chains of 18 atoms each, then `create_bonds` for the 17 bonds per chain (use `create_bonds many` with a tight distance cut to bond only nearest neighbors along each chain). Document the parameters in inline comments matching `make_pe_box()` in `polyethylene_tension.cpp`.
- [ ] Reduce the `run` step counts in the 5 stages to {1000, 5000, 5000, 5000, 5000} to match acceptance-verification scale rather than paper-scale validation; leave a comment noting paper-scale runs are tracked in `msu-pe-validation`.
- [ ] Add at the top (after `thermo_style`) a `fix thermolog all print 100 "${step} ${temp} ${press} ${pe} ${ke} ${pxx} ${pyy} ${pzz}" file thermo_reference.csv title "step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz" screen no` block.
- [ ] In stage 5 (tension), add `variable strain equal "(lx - v_lx0) / v_lx0"` and `fix stressout all print 1000 "${strain} ${pxx}" file stress_strain_reference.csv title "strain,sigma_xx" screen no`; set `variable lx0 equal lx` once at the top of stage 5.
- [ ] Run `cd benchmarks/md/reference/polyethylene && /opt/homebrew/bin/lmp_mpi -in lammps.in -log log.lammps` and verify (a) exit 0, (b) both CSVs are written in this directory, (c) each CSV has ≥ 11 lines (header + ≥ 10 data rows).
- [ ] Spot-check `thermo_reference.csv` line 2: `T` should be within 50 K of 500 (stage 1 setpoint); line `(stage1_end / 100) + 1` should show `T` near 500 with `P` finite.
- [ ] Spot-check `stress_strain_reference.csv` final line: `strain` ≥ 0.5 (we run enough tension steps for non-trivial strain).
- [ ] Update `README.md` version block: replace `<TBD: e.g. stable_29Aug2024>` with `stable_22Jul2025_update4` (Homebrew Cellar path `/opt/homebrew/Cellar/lammps/20250722-update4`), `<TBD: compiler>` with the actual `lmp_mpi -h` output's compiler line, and `<TBD: CPU model...>` with `sysctl -n machdep.cpu.brand_string` + core count from `sysctl -n hw.ncpu`.
- [ ] Add a short note under "Reproducing the LAMMPS reference run" recording the step-count reduction (paper-scale runs deferred to `msu-pe-validation`).
- [ ] Commit all four files in one commit: `bench(parallel-md-finish-02): self-contained lammps.in + reference CSVs + README annotation`.

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
