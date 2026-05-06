# MSU amorphous-polyethylene reference data

This directory holds the LAMMPS reference data and reproduction recipe for
the parallel-md-engine Phase 5 acceptance gates (see
`.claude/specs/parallel-md-engine.md` Phase 5 + `.claude/specs/msu-pe-validation.md`).

## What's here

| File | Status | Description |
|---|---|---|
| `lammps.in` | committed | LAMMPS input mirroring `benchmarks/md/polyethylene_tension.cpp` 5-stage protocol. |
| `msu_pe.data` | **TBD** | Initial geometry. 20 chains × 1000 united-atom CH₂ on FCC seed. Generate via moltemplate / PSP — instructions below. |
| `thermo_reference.csv` | **TBD** | Populated by running `lammps.in`. Columns: `step T P E_pot E_kin Pxx Pyy Pzz`. |
| `stress_strain_reference.csv` | **TBD** | Populated by `lammps.in` stage 5. Columns: `strain sigma_xx`. |
| `profile_1.csv` / `profile_8.csv` / `profile_32.csv` | **TBD** | Atomiverse per-component profile. `benchmarks/md/pe_profile_single.cpp` writes these. |
| `scaling_strong.png` / `scaling_weak.png` | **TBD** | `tools/plot_scaling.py` writes these from `scaling.cpp` CSV. |
| `profile_8.png` | **TBD** | `tools/plot_profile.py` writes this from `profile_8.csv`. |
| `scaling_report.md` | **TBD** | Human-readable scaling/profile writeup; populated by Phase 5 Track B.1, B.6, B.8, B.9. |

The `**TBD**` rows require external tooling (LAMMPS install, geometry generation)
or actual benchmark runs on the recorded reference machine. They are not generated
by `make` — see the reproduction recipe below.

## Reproducing the LAMMPS reference run

### Required versions

The reference numbers in `thermo_reference.csv` and `stress_strain_reference.csv`
were produced with:

- **LAMMPS git tag**: `<TBD: e.g. stable_29Aug2024>` (record exact tag here when
  the reference run is generated).
- **Compiler**: `<TBD: e.g. gcc 13.2.0>`.
- **Build flags**: `cmake -DBUILD_MPI=on -DPKG_MOLECULE=on -DPKG_RIGID=on
  -DCMAKE_CXX_FLAGS="-O3 -march=native"`.
- **Hardware fingerprint**: `<TBD: CPU model, # sockets, # cores, RAM>`.

### Steps

```bash
# 1. Generate msu_pe.data (one-off; commit the result).
#    Two options:
#    (a) Direct download from MSU ICME wiki supplementary material:
#        https://icme.hpc.msstate.edu/mediawiki/index.php/Deformation_of_Amorphous_Polyethylene
#    (b) Regenerate via moltemplate (preferred for reproducibility):
#        molc=20  monomers=1000  rho=0.94 g/cm^3
#        FCC seed at 1.53 A, random-coil relaxation, NVT @ 500 K.
#    Place result here as msu_pe.data.

# 2. Run LAMMPS reference.
mpirun -n 8 lmp -in lammps.in -log log.lammps

# 3. The lammps.in script writes thermo_reference.csv and
#    stress_strain_reference.csv directly into this directory.

# 4. Annotate the version block at the top of this README with the exact
#    LAMMPS tag, compiler version, build flags, and hardware fingerprint
#    used. Commit the resulting CSVs alongside this README.
```

## Reproducing the Atomiverse run

```bash
cmake -B build -DATV_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Single-node multi-rank smoke (default scale=1; production runs use --scale 1
# for the spec-canonical step counts):
mpirun -n 8 ./build/benchmarks/polyethylene_tension --mode full --out-dir runs/atomiverse_n8

# Standardised launchers (see tools/run_msu_*.sh):
./tools/run_msu_atomiverse.sh 8
./tools/run_msu_lammps.sh     8
```

## Thread binding

Atomiverse Phase 5 Track B.6 audits OpenMP thread placement. Recommended env
for both engines on a single-socket box:

```bash
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

For a 2-socket NUMA box, also pin the MPI ranks to NUMA nodes:

```bash
mpirun --bind-to numa --map-by numa --rank-by node -n 8 ...
# or under SLURM:
srun --cpu-bind=cores --hint=memory_bound -n 8 ...
```

### Measured NUMA effect

`<TBD: record before/after wall-time numbers from the reference machine here.
On a single-socket dev box this entry should read "single-socket — no NUMA
effect observed", with the measurement command preserved for later replay.>`

## Post-processing

```bash
# Stress-strain comparison + modulus / yield / hardening fits:
python tools/validate_polyethylene.py \
    --reference benchmarks/md/reference/polyethylene/stress_strain_reference.csv \
    --candidate runs/atomiverse_n8/stress_strain.csv

# Strong/weak scaling plots:
python tools/plot_scaling.py runs/scaling_n*.csv \
    --out-dir benchmarks/md/reference/polyethylene/

# Per-component profile breakdown:
python tools/plot_profile.py \
    benchmarks/md/reference/polyethylene/profile_8.csv
```

## Tolerances (Phase 5 acceptance)

| Observable | Tolerance vs LAMMPS |
|---|---|
| Per-step thermo (first 1000 steps) | 1e-4 relative on T, E_pot, E_kin |
| Stage-4 ⟨ρ⟩, ⟨E_pot⟩ | 1% / 0.5% |
| Stress-strain modulus E | ±5% |
| Stress-strain yield σ_y | ±5% |
| Stress-strain hardening dσ/dε | ±5% |
| Single-core wall time / step | ≤ 1.2× LAMMPS |
| 8-rank wall time / step | ≤ 2.0× LAMMPS |

Source: `.claude/specs/parallel-md-engine.acceptance.md` ac-022, ac-023, ac-024.

## Citations

- MSU ICME polymer-deformation reference:
  https://icme.hpc.msstate.edu/mediawiki/index.php/Deformation_of_Amorphous_Polyethylene
- LAMMPS: A. P. Thompson et al., *Comp. Phys. Comm.* 271 (2022) 108171.
- Dreiding force field: Mayo et al., *J. Phys. Chem.* 94 (1990) 8897.
