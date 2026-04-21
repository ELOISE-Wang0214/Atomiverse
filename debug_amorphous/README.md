# debug_amorphous — Atomiverse vs LAMMPS amorphous-PE comparison

Reproduces the
[MSU amorphous polyethylene tension protocol][msu] (Stages 1 + 2 from
`amourphous.in`) on the Atomiverse CPU + MPI backend so the trajectories
can be put next to LAMMPS line-by-line.

[msu]: https://icme.hpc.msstate.edu/mediawiki/index.php/Deformation_of_Amorphous_Polyethylene.html

## Files

| File | Role |
|---|---|
| `PE_nc10_cl1000.dat` | LAMMPS data file: 10 chains × 1000 united atoms, 80.06 Å³ box |
| `amourphous.in` | LAMMPS input script (DREIDING-like force field, NVE/limit + Langevin Stage 1, NPT Stage 2) |
| `log.PE_nc10_cl1000.txt` | LAMMPS reference log (Stages 1 + 2 already run) |
| `lammps_loader.h` | C++ loader: LAMMPS data → `cpu::Frame` |
| `gather_frame.h` | Inter-stage all-gather so Stage 2 can resume from Stage 1's multi-rank state |
| `run_atomiverse.cpp` | Main MPI driver |
| `run.sh` | Build-and-launch helper |
| `compare.py` | Side-by-side T / P / Lx comparison + plot |
| `CMakeLists.txt` | Subdirectory build (gated on `-DATV_BUILD_DEBUG_AMORPHOUS=ON`) |

## Force-field mapping

Both engines use V = K(r − r₀)² (no ½ factor), V = K(θ − θ₀)², V = Σ Cₙ
cosⁿ(φ).  All `kcal/mol` ↔ `eV` and `degree` ↔ `rad` conversions live
inside `run_atomiverse.cpp` so the LAMMPS `.in` numbers can be pasted in
verbatim.

| | LAMMPS (real units) | Atomiverse (eV/Å/fs/amu/K) |
|---|---|---|
| Bond | K=350, r₀=1.53 | K=15.18 eV/Å², r₀=1.53 |
| Angle | K=60, θ₀=109.5° | K=2.602 eV/rad², θ₀=1.911 rad |
| Dihedral | A=1.73,−4.49,0.776,6.99,0 (kcal/mol) | C=0.0750,−0.1947,0.0337,0.3032,0 (eV) |
| LJ | ε=0.112, σ=4.01, cut=10.5 | ε=4.857e-3 eV, σ=4.01 Å, cut=10.5 Å |

## Stage protocol

| Stage | LAMMPS | Atomiverse mapping |
|---|---|---|
| 1 | `fix nve/limit 0.05` + `fix langevin T=500 τ=10 fs`, dt=1 fs, 10 000 steps | `cpu::md::NVELimit{1, 0.05}` + `cpu::md::FixLangevin{500, 10, 1, 904297}` |
| 2 | `fix npt T=500 P=0 τ_T=50 τ_P=1000`, dt=0.5 fs, 50 000 steps | `cpu::md::NoseHooverNPT{500, 100, 0, 1000, 0.5}` (single-chain MTK isotropic) |

## Build

```bash
cmake -B build -DATV_USE_MPI=ON -DATV_BUILD_DEBUG_AMORPHOUS=ON
cmake --build build --target run_atomiverse_amorphous -j
```

## Run

```bash
# Smoke test (1 000 + 500 steps, 4 ranks)
debug_amorphous/run.sh --np 4 --steps-stage1 1000 --steps-stage2 500

# Full LAMMPS-equivalent protocol on 8 ranks (~5 minutes)
debug_amorphous/run.sh --np 8

# Stage 1 only
debug_amorphous/run.sh --np 4 --stage1-only --steps-stage1 10000

# Single-rank reference
debug_amorphous/run.sh --np 1
```

The runner writes thermo every 100 steps to `log.atomiverse.txt` in a
LAMMPS-comparable format:

```
#       Step          Temp        Press[atm]      PotE[eV]      KinE[eV]        Lx        Ly        Lz
         0     8590.97  -1.7663e+08    14311222   11103.59    80.0586    80.0586    80.0586
       100      972.96  -1.7250e+11  1.3995e+10   1257.52    80.0586    80.0586    80.0586
       ...
```

## Compare

```bash
python3 debug_amorphous/compare.py
```

Output: side-by-side LAMMPS vs Atomiverse T(step), P(step), Lx(step) plots
in `compare.png` plus a text summary table on stdout.

## Known caveats

- **Stage 1 transients**: LAMMPS step 0 is reported BEFORE any
  integration; Atomiverse step 0 is reported AFTER the first step and
  reflects the strong relaxation kick from the over-compressed initial
  state. The two engines agree closely from step ≳ 200 onward.
- **Stage 2 transients**: Atomiverse's `NoseHooverNPT` is a
  single-chain MTK NHC; LAMMPS uses 3-chain by default. Single-chain
  oscillates more visibly during equilibration (T spikes of a few
  thousand K early on, plateau-and-overshoot pattern around the
  thermostat half-period of `tau_T`). Equilibrium ⟨T⟩ and ⟨P⟩ should
  match the LAMMPS canonical values within statistical error after
  ~5 × tau_T integration time.
- **Pressure sign**: Atomiverse stores stress as
  σ_phys×(−1) internally; the runner inverts on output so the printed
  number is in the LAMMPS sign convention (positive = compression).
- **Dynamics divergence**: LAMMPS Langevin and Atomiverse's custom
  Langevin use independent random streams, so individual trajectories
  differ; only ensemble averages are expected to match.
- **Inter-stage handoff**: each driver creates its own
  `DomainDecomposition`, so we explicitly `MPI_Allgatherv` the per-rank
  owners back to every rank between stages (`gather_frame.h`).
