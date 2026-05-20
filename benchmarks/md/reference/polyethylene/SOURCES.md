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
