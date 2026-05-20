# MSU amorphous-polyethylene reference data

This directory holds the **MSU CAVS ICME amorphous-PE LAMMPS reference**
— the canonical ground-truth input deck + initial geometry for the
Tschopp 2010 / Hossain–Tschopp et al. 2010 *Polymer* paper. Used by
`.claude/specs/msu-pe-validation.md` Phase C–D.

## What's here

| File | Description |
|---|---|
| `lammps_msu_paper.in` | Original `in.deform_dreiding.txt` from MSU CAVS ICME — the canonical 5-stage Dreiding-UA uniaxial-tension protocol. Author: Mark Tschopp 2010. Do not modify. |
| `lammps.in` | Phase D comparison deck — physics byte-for-byte identical to `lammps_msu_paper.in`, plus non-physics canonical-CSV / dump / RDF output hooks. See `SOURCES.md` § Phase D.1. |
| `msu_pe.data` | Original `PE_nc10_cl1000.dat` — 10 chains × 1000 united atoms = 10 000 atoms, 80.0586 Å cubic box. Initial geometry. Do not modify. |
| `SOURCES.md` | Provenance + attribution + the FF discrepancy table that triggered `pe-ff-migrate-dreiding-ua`. |

Both source files are bit-identical to the upstream content
(HTML wrappers stripped only). See `SOURCES.md`.

## Ground truth

**The MSU files in this directory are the reference**. When the
Atomiverse PE force-field code disagrees with `lammps_msu_paper.in`,
the Atomiverse code is wrong — that's the migration spec
`pe-ff-migrate-dreiding-ua` exists to fix.

Force-field parameters live in `lammps_msu_paper.in`:

```
bond_coeff      1  350     1.53          # k [kcal/mol/Å²], r0 [Å]
angle_coeff     1   60   109.5            # k [kcal/mol/rad²], θ0 [°]
dihedral_coeff  1  1.73 -4.49 0.776 6.99 0.0   # Dreiding-UA (Tschopp 2010)
pair_coeff      1 1  0.112   4.01          # ε [kcal/mol], σ [Å]
special_bonds   lj 0.0 0.0 0.0
neighbor        2.0 bin
neigh_modify    delay 5 every 1 check yes
```

## Reproducing the LAMMPS reference run

```bash
cd benchmarks/md/reference/polyethylene
/opt/homebrew/bin/lmp_mpi -in lammps_msu_paper.in -log log.lammps
```

The deck runs the full 5-stage protocol (Langevin NVT-500 → NPT-500 →
NVT-cool → NPT-100 → uniaxial tension along x) at paper scale on the
10 000-atom system; takes minutes to hours depending on rank count.

Provenance, citation, and the FF-parameter table are in `SOURCES.md`.

## Phase D — Atomiverse vs LAMMPS reproduction

`lammps.in` (this directory) and the Atomiverse benchmark
`benchmarks/md/polyethylene_tension.cpp` run the **same** 5-stage
protocol on `msu_pe.data`. The full ~360000-step run is ~24h-class wall
time — an **operator action, not a CI step**. Reproduce the reference
artifacts at 1 and 8 ranks with:

```bash
./tools/run_msu_lammps.sh 1
./tools/run_msu_lammps.sh 8
./tools/run_msu_atomiverse.sh 1
./tools/run_msu_atomiverse.sh 8
```

Each launcher publishes engine- and rank-suffixed canonical CSVs
(`thermo_<engine>_<N>r.csv`, `stress_strain_<engine>_<N>r.csv`) into
this directory; the Phase D.2/D.3/D.4 comparison
(`tools/validate_polyethylene.py`) consumes them by path. See
`SOURCES.md` § Phase D.1 for the protocol table, the canonical filename
scheme, and the full operator procedure.
