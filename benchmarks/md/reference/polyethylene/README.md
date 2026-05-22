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
protocol on `msu_pe.data`, stage-for-stage. The benchmark mirrors the
LAMMPS deck exactly — each stage uses the integrator, ensemble and
timestep of the matching `lammps.in` stage:

| Stage | lammps.in                        | Atomiverse components                          | dt     |
|-------|----------------------------------|------------------------------------------------|--------|
| 1     | `nve/limit 0.05` + `langevin 500`| `NVELimit` + `FixLangevin` + `InitVelocity`    | 1.0 fs |
| 2     | `npt temp 500 500 iso 0`         | isotropic `NoseHooverNPT` (`T_stop==T_start`)  | 0.5 fs |
| 3     | `npt temp 500 100 iso 0`         | `NoseHooverNPT`, continuous 500→100 K ramp     | 0.5 fs |
| 4     | `npt temp 100 100 iso 0`         | isotropic `NoseHooverNPT` (`T_stop==T_start`)  | 0.5 fs |
| 5     | `npt temp 100 100 y/z 0` + `deform x` | anisotropic `NoseHooverNPT` (`couple_axis={false,true,true}`) + `Deform` x | 1.0 fs |

Because `dt` is fundamental to each integrator kernel (no `dt` setter)
and the per-stage timesteps differ (1.0/0.5/0.5/0.5/1.0 fs), the
benchmark **builds a fresh `cpu::md::MD` driver per stage**, runs that
stage, then discards the driver. The `cpu::Frame` (positions,
velocities, box) carries state across stages; the shared
`CsvRecorder` / `TrajectoryRecorder` / optional `RDF` are re-registered
into every per-stage driver so the canonical CSV step column stays
globally monotonic across the stage boundaries.

Stage 1 uses `NVELimit` (per-step displacement cap) rather than plain
Velocity-Verlet so the unrelaxed-configuration atom overlap
(E_pair ~3.4e8 at step 0) is tamed instead of blowing the integrator
up. Stage 3 is a **single** `NoseHooverNPT` with `T_start=500,
T_stop=100` and a run span covering the stage — the integrator
interpolates the setpoint internally; there is no cooling micro-segment
loop. Stage 5 holds the x box edge rigid under the NPT barostat
(`couple_axis={false,true,true}`) while `Deform` drives the x edge at
the engineering strain rate and y/z barostat to P=0.

The full **331800-step** run (10000 + 50000 + 50000 + 50000 + 171800)
is ~24h-class wall time — an **operator action, not a CI step**.

The LAMMPS launcher `tools/run_msu_lammps.sh` defaults to the user's
single `lmp_cuda` binary run in **CPU mode** — the CPU arm passes no GPU
package args, so `lmp_cuda` runs as a plain MPI/OpenMP CPU build;
`LMP_BIN` still overrides the binary. Reproduce the reference artifacts
at 1 and 8 ranks on a local node with:

```bash
./tools/run_msu_lammps.sh 1
./tools/run_msu_lammps.sh 8
./tools/run_msu_atomiverse.sh 1
./tools/run_msu_atomiverse.sh 8
```

On **PDC Dardel**, submit each run as a SLURM job on the `main`
partition with the submission wrapper `tools/run_msu_pdc.sh` (account
`naiss2025-1-25` for the CPU arm; the `-A` flag is required, no default):

```bash
./tools/run_msu_pdc.sh lammps     1 -A naiss2025-1-25
./tools/run_msu_pdc.sh lammps     8 -A naiss2025-1-25
./tools/run_msu_pdc.sh atomiverse 1 -A naiss2025-1-25
./tools/run_msu_pdc.sh atomiverse 8 -A naiss2025-1-25
```

`run_msu_pdc.sh` generates a `main`-partition, 24h-cap SLURM script that
runs the matching launcher under `srun -n <ranks>` and `sbatch`es it.

Each launcher routes raw engine output through the shared, pytest-covered
converter `tools/canonicalize_msu_output.py` and publishes engine- and
rank-suffixed canonical CSVs (`thermo_<engine>_<N>r.csv`,
`stress_strain_<engine>_<N>r.csv`) into this directory; the Phase
D.2/D.3/D.4 comparison (`tools/validate_polyethylene.py`) consumes them
by path. See `SOURCES.md` § Phase D.1 for the protocol table, the
331800-step total, the canonical filename scheme, the
`canonicalize_msu_output.py` step, and the full operator procedure.
