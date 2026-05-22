# Upstream provenance of MSU paper-scale reference data

Two files in this directory come from the Mississippi State CAVS ICME wiki
and are the **canonical Tschopp 2010 / Hossain–Tschopp 2010 Polymer**
reference inputs for amorphous-polyethylene uniaxial tension:

- `msu_pe.data` — original filename `PE_nc10_cl1000.dat`. Saved from
  <https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/files/PE_nc10_cl1000.dat>
  on 2026-05-19. HTML wrapper stripped (lines bracketed by
  `<html><head>…</body></html>` removed); LAMMPS data payload is byte-
  identical to the upstream content.

- `lammps_msu_paper.in` — original filename `in.deform_dreiding.txt`.
  Inlined as a `<pre>` block on
  <https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>.
  HTML entities (`&#60;`, `&amp;`) decoded; tab/space whitespace
  preserved.

### System parameters (`msu_pe.data`)

| Item | Value |
|---|---|
| Number of chains | 10 |
| Monomers per chain | 1 000 |
| Total united atoms | 10 000 |
| Box (cubic, all axes) | 80.0586 Å |
| Atom type | 1 (CH₂ united atom, mass 14.02 amu) |
| Bonds / angles / dihedrals | 9 990 / 9 980 / 9 970 |

### Force field (`lammps_msu_paper.in`)

```
bond_style      harmonic
bond_coeff   1  350     1.53          # k [kcal/mol/Å²], r0 [Å]
angle_style     harmonic
angle_coeff  1   60   109.5            # k [kcal/mol/rad²], θ0 [°]
dihedral_style  multi/harmonic
dihedral_coeff 1  1.73 -4.49 0.776 6.99 0.0   # Dreiding-UA Tschopp 2010
pair_style      lj/cut 10.5
pair_coeff   1 1  0.112   4.01          # ε [kcal/mol], σ [Å]
```

### Protocol (5 stages, fixed by upstream script)

1. NVE/limit + Langevin @ 500 K, dt = 1.0 fs, 10 000 steps
2. NPT @ 500 K, dt = 0.5 fs, 50 000 steps
3. NVT cool 500 → 100 K, 50 000 steps
4. NPT @ 100 K, dt = 0.5 fs, 50 000 steps
5. Uniaxial tension at ε̇ = 1 × 10⁻⁵/timestep along x (≈ 2 × 10⁻⁵/fs at
   dt = 0.5 fs), 171 800 steps to ε ≈ 1.0

### Attribution + citation

Author: **Mark Tschopp, 2010**.

Cite when reporting results derived from this data:

> M.A. Tschopp, M.F. Horstemeyer, et al.; M.K. Hossain, M.A. Tschopp,
> et al., "Molecular dynamics simulations of deformation in amorphous
> polyethylene," *Polymer* (2010).

CAVS ICME source page (with method writeup and figures):
<https://www.cavs.msstate.edu/icme/code/lammps/atomistic_research/polymers/amorphous.php>

### Atomiverse FF discrepancy — resolved 2026-05-19

`benchmarks/md/polyethylene_tension.cpp` and 14 PE-related tests
previously used wrong (OPLS-AA-like) angle and dihedral parameters
that did not match this directory's MSU reference. Resolved in commit
`pe-ff-migrate-dreiding-ua`: the 15 source / test files now use the
Dreiding-UA constants verbatim from `lammps_msu_paper.in`.

#### Historical mismatch (pre-2026-05-19, kept for context)

| Term | Atomiverse code (before fix) | MSU paper (ground truth) |
|---|---|---|
| `bond r0` | 1.529 Å | 1.53 Å |
| `angle k` | 58.35 kcal/mol/rad² | **60 kcal/mol/rad²** |
| `angle θ0` | 112.7° (OPLS-AA-like) | **109.5°** (sp³ tetrahedral Dreiding-UA) |
| `dihedral K1..K4` | 2.0058, -2.5128, 0.2710, 6.2900 (OPLS) | **1.73, -4.49, 0.776, 6.99** (Dreiding-UA) |
| `LJ ε, σ` | 0.112, 4.01 | 0.112, 4.01 ✓ |

After the fix:
- Atomiverse FF and `lammps_msu_paper.in` exercise the same physics.
- `msu-pe-validation.md` Phase D (quantitative ±5% LAMMPS comparison)
  is unblocked.
- One test (`test_pe_single_chain.cpp::AllTransGeometryEnergyIsZero`)
  had its V≈0 assertion adjusted to the Dreiding-UA constant offset at
  trans (~4.4e-3 eV for 17 dihedrals; was an OPLS-only artifact).
- New pin test `tests/cpu/md/test_pe_ff_constants.cpp` guards against
  future drift.

**Do not edit the MSU files in this directory to "match" Atomiverse**.
Ground truth flows MSU → Atomiverse, never the other way.

---

## C20H42 bonded-energy cross-check reference (spec c20h42-energy-xcheck)

`c20h42_chain.data` + `in.c20h42_energy` are the LAMMPS reference for the
single-chain bonded-energy cross-check test
`tests/cpu/md/test_pe_energy_xcheck.cpp`.

- **Geometry** — a 20-atom all-trans united-atom C20H42 zigzag
  (r0 = 1.53 Å, θ0 = 109.5°) with a closed-form deterministic per-atom
  perturbation so every bond/angle/dihedral term is non-trivial
  (a perfect-equilibrium chain gives E_bond = E_angle = 0). Atom i
  (0-indexed) is displaced by
  `0.10 Å · (sin(2.4i+0.7), sin(3.1i+1.3), sin(1.7i+2.1))`.
  The same 20×3 coordinates are embedded as a `constexpr` array in the
  test — single source of truth, both engines evaluate identical coords.
- **Force field** — MSU Dreiding-UA constants verbatim from
  `lammps_msu_paper.in`: bond harmonic K=350 kcal/mol/Å² r0=1.53 Å;
  angle harmonic K=60 kcal/mol/rad² θ0=109.5°; dihedral multi/harmonic
  A1..A5 = 1.73 -4.49 0.776 6.99 0.0 kcal/mol.
- **Regenerate** the data file: the perturbation formula above is fully
  specified; a one-off script reproduces `c20h42_chain.data` exactly.
- **Regenerate** the reference energies:
  `cd benchmarks/md/reference/polyethylene && lmp_serial -in in.c20h42_energy`
  → reads the `Step E_bond E_angle E_dihed` thermo line at step 0.
- **Reference energies** (LAMMPS 22 Jul 2025, `run 0`, kcal/mol):
  `E_bond = 122.63225361801`, `E_angle = 30.774982565814`,
  `E_dihedral = 3.5122664779742`.
- **Unit conversion** — LAMMPS `units real` is kcal/mol; Atomiverse
  energy is eV. The test converts with
  `1 kcal/mol = 1 / 23.060547830619026 eV` (≈ 0.0433641043 eV).

---

## Phase D.1 — Atomiverse-vs-LAMMPS comparison inputs (spec `msu-pe-phase-d-01-reference-runs`)

`lammps.in` is the LAMMPS deck used for the Phase D quantitative
comparison. Its **physics — every fix, integrator, ensemble, timestep,
force-field coefficient and neighbor setting — is byte-for-byte the
canonical upstream deck `lammps_msu_paper.in`** (Tschopp 2010, kept as
immutable ground truth). `lammps.in` differs from it ONLY in non-physics
output hooks:

- a persistent `fix print` that streams a canonical-schema per-step CSV
  `thermo.lammps.csv`, header `step,T,P,E_pot,E_kin,Pxx,Pyy,Pzz`;
- a stage-5 `fix print` streaming `stress_strain.lammps.csv`, header
  `strain,sigma_xx`;
- an XYZ trajectory dump every 10000 steps and an accumulated RDF.

Nothing about the dynamics differs from the reference. Only the `run`
step counts may be shortened for a quick check; the stages and
simulation parameters must never change.

### Protocol (`lammps.in` ≡ `lammps_msu_paper.in`)

| Stage | Fix / ensemble | T | P | Steps | dt |
|---|---|---|---|---|---|
| 1 | `nve/limit 0.05` + `langevin` (overlap relaxation) | 500 K | — | 10000 | 1.0 fs |
| 2 | `npt ... iso` | 500 K | 0 bar | 50000 | 0.5 fs |
| 3 | `npt temp 500 100 ... iso` (cooling ramp) | 500 → 100 K | 0 bar | 50000 | 0.5 fs |
| 4 | `npt ... iso` | 100 K | 0 bar | 50000 | 0.5 fs |
| 5 | `npt ... y ... z` + `deform 1 x erate` | 100 K | 0 bar (y,z) | 171800 | 1.0 fs |

Stage 1 uses `nve/limit` because `msu_pe.data` is a built-but-unrelaxed
configuration with non-bonded overlaps (initial `E_pair ≈ 3.4e8
kcal/mol`); the displacement-capped integrator squeezes them out
gently. A plain `nvt`/Velocity-Verlet integrator blows the system up —
the stage-1 integrator choice is load-bearing, not cosmetic.

Force-field coefficients are identical to `lammps_msu_paper.in` and the
FF table above.

### Total step count

The five-stage protocol is exactly **331800 steps**:

```
10000 + 50000 + 50000 + 50000 + 171800 = 331800
```

(stage 1 + stage 2 + stage 3 + stage 4 + stage 5). An earlier "~360000"
figure here was an overestimate and has been corrected — the `run`
step-counts above are the authoritative source. The full run is still
~24h-class wall time, inside the PDC 24h job cap.

### Canonical reference artifacts (filename scheme)

`tools/run_msu_lammps.sh` and `tools/run_msu_atomiverse.sh` publish their
output as engine- and rank-suffixed CSVs directly under this directory:

| File | Producer |
|---|---|
| `thermo_lammps_<N>r.csv` | LAMMPS, N ranks |
| `thermo_atomiverse_<N>r.csv` | Atomiverse, N ranks |
| `stress_strain_lammps_<N>r.csv` | LAMMPS, N ranks |
| `stress_strain_atomiverse_<N>r.csv` | Atomiverse, N ranks |

Both run scripts route raw engine output through the shared, pytest-covered
converter `tools/canonicalize_msu_output.py` (no inline `awk`/`cp`), so
LAMMPS and Atomiverse use one tested canonicalization implementation:

```bash
# raw engine output  ->  canonical-schema CSV
python3 tools/canonicalize_msu_output.py <engine> <kind> <in> <out>
#   engine ∈ {lammps, atomiverse}   kind ∈ {thermo, stress-strain}
```

The Atomiverse benchmark currently logs an isotropic pressure only; the
`Pxx/Pyy/Pzz` columns of `thermo_atomiverse_*.csv` are filled with that
isotropic value as a **Phase G.2 stand-in** (umbrella spec
`msu-pe-validation.md` Phase G.2 wires per-axis pressure into
`benchmarks/md/polyethylene_tension.cpp`). `canonicalize_msu_output.py`
preserves this stand-in deliberately. Downstream consumers of per-axis
Atomiverse stress must treat those three columns as isotropic placeholders
until Phase G.2 lands.

### Staged-landing decision (per spec `msu-pe-phase-d-01-reference-runs`)

- **Now** — `stress_strain_lammps_8r.csv` is committed: it was produced by
  `canonicalize_msu_output.py` from the already-completed 8-rank paper-deck
  run in `runs/lammps_n8_paper/PE_nc10_cl1000.def1.txt`, giving
  `msu-pe-phase-d-04`'s load-bearing stress-strain gate its first real
  LAMMPS-side reference data without a re-run.
- **After the operator runs** — `thermo_lammps_<N>r.csv` and
  `stress_strain_lammps_<N>r.csv` (N ∈ {1,8}) from the four full PDC runs.
- **After Phase G.2** — `thermo_atomiverse_<N>r.csv` /
  `stress_strain_atomiverse_<N>r.csv` with real per-axis pressure. If an
  Atomiverse reference is run before G.2, its `Pxx/Pyy/Pzz` are isotropic
  stand-ins and must be flagged as such here.

### Deck-fidelity audit

`tools/audit_lammps_deck.sh` is the machine-verifiable gate behind the
"`lammps.in` physics ≡ `lammps_msu_paper.in`" claim above: it strips each
deck's non-physics output hooks (`fix print` / `dump` / `compute rdf` /
`log` / `write_restart` / `reset_timestep` / `thermo_style` ...) and
asserts the remaining physics lines — units, force field, every stage's
fixes / timesteps / run step-counts — are set-equal. It exits 0 on the
in-tree decks.

### Operator run procedure (NOT a CI step)

The full 331800-step protocol is ~24h-class wall time — an **operator
action, never a CI step**.

The LAMMPS launcher defaults to the user's single `lmp_cuda` binary run in
**CPU mode** (no GPU package args — plain MPI/OpenMP); `LMP_BIN` still
overrides it. On a local node:

```bash
./tools/run_msu_lammps.sh 1       # LAMMPS reference, 1 rank
./tools/run_msu_lammps.sh 8       # LAMMPS reference, 8 ranks
./tools/run_msu_atomiverse.sh 1   # Atomiverse reference, 1 rank
./tools/run_msu_atomiverse.sh 8   # Atomiverse reference, 8 ranks
```

On PDC Dardel, submit each as a SLURM job on the `main` partition via the
submission wrapper `tools/run_msu_pdc.sh` (account `naiss2025-1-25` for
the CPU arm; the `-A` flag is required, no default):

```bash
./tools/run_msu_pdc.sh lammps     1 -A naiss2025-1-25
./tools/run_msu_pdc.sh lammps     8 -A naiss2025-1-25
./tools/run_msu_pdc.sh atomiverse 1 -A naiss2025-1-25
./tools/run_msu_pdc.sh atomiverse 8 -A naiss2025-1-25
```

`run_msu_pdc.sh` generates a `main`-partition, 24h-cap SLURM script whose
body runs the matching launcher under `srun -n <ranks>`, then `sbatch`es
it. Each run also leaves a timestamped audit copy under `runs/`. Commit the
eight canonical CSVs once all four runs complete.
