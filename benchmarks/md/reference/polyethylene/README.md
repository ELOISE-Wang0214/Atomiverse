# MSU amorphous-polyethylene reference data

This directory holds the LAMMPS reference data and reproduction recipe for
the parallel-md-engine Phase 5 acceptance gates (see
`.claude/specs/parallel-md-engine.md` Phase 5 + `.claude/specs/msu-pe-validation.md`).

## What's here

| File | Status | Description |
|---|---|---|
| `lammps.in` | committed | LAMMPS input mirroring `benchmarks/md/polyethylene_tension.cpp` 5-stage protocol. Self-contained — no external `msu_pe.data` needed; chain geometry comes from `pe_chain.mol` (molecule template) instantiated on a 4×4×4 lattice. |
| `pe_chain.mol` | committed | LAMMPS molecule template — 18 UA carbons, 17 bonds, 16 angles, 15 dihedrals. Geometry mirrors `make_pe_box()` in `polyethylene_tension.cpp` (zig-zag along +y, dy = 1.272798 Å, dz = 0.847247 Å). |
| `thermo_reference.csv` | committed | Populated by running `lammps.in`. Columns: `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`. 230 rows spanning the full 5-stage run. |
| `stress_strain_reference.csv` | committed | Populated by `lammps.in` stage 5. Columns: `strain,sigma_xx`. 50 rows from the tension stage covering ε ∈ [0, ~0.37]. |
| `profile_1.csv` / `profile_8.csv` | committed | Atomiverse per-component profile from `benchmarks/md/pe_profile_single.cpp` at 1 / 8 ranks on this host. |
| `profile_1.png` / `profile_8.png` | committed | `tools/plot_profile.py` stacked-bar PNGs. |
| `profile_32.csv` / `scaling_strong.png` / `scaling_weak.png` | **TBD** | Multi-node 32-rank profile + strong/weak scaling — needs cluster access. |
| `scaling_report.md` | **TBD (chain spec 04)** | Human-readable scaling/profile writeup; populated by `parallel-md-finish-04-scaling-report`. |

## Reproducing the LAMMPS reference run

### Required versions

The reference numbers in `thermo_reference.csv` and `stress_strain_reference.csv`
were produced with:

- **LAMMPS version**: *Large-scale Atomic/Molecular Massively Parallel Simulator — 22 Jul 2025 - Update 4* (`lammps-stable_22Jul2025_update4`), packaged in Homebrew Cellar at `/opt/homebrew/Cellar/lammps/20250722-update4`.
- **MPI**: Homebrew Open-MPI (`/opt/homebrew/bin/mpirun`).
- **Compiler / build flags**: as packaged by Homebrew (`PKG_MOLECULE` + `PKG_RIGID` + `PKG_KSPACE` enabled; AppleClang `-O3 -march=native`).

#### Hardware fingerprint

- CPU: Apple **M4** (`sysctl -n machdep.cpu.brand_string` output).
- Cores: 10 physical / 10 logical.
- OS: macOS 26.4.1.

### Step-count deviation from MSU paper

For acceptance-verification purposes, stage step counts are reduced from the
paper-scale schedule {10 000, 50 000, 50 000, 50 000, 200 000} to {1 000, 5 000,
5 000, 5 000, 5 000} (≈ 1 s wall-time on this host). The stage-5 strain rate is
bumped 150× from the paper rate (1 × 10⁻⁶/fs → 1.5 × 10⁻⁴/fs) so the shortened
run still covers ε ∈ [0, ~0.37]. Paper-scale runs are owned by
`.claude/specs/msu-pe-validation.md`, not this directory.

### Steps

```bash
cd benchmarks/md/reference/polyethylene
/opt/homebrew/bin/lmp_mpi -in lammps.in -log log.lammps

# Writes (into this directory):
#   thermo_reference.csv         (230 rows, schema in the table above)
#   stress_strain_reference.csv  (50 rows)
#   log.lammps                   (full LAMMPS log; not committed)
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
