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

### Known discrepancies vs. current Atomiverse code (2026-05-19)

`benchmarks/md/polyethylene_tension.cpp` and the smaller
`lammps.in` self-contained reference (committed by
`parallel-md-finish-02`) use **different angle and dihedral
parameters** than the MSU paper:

| Term | Atomiverse current | MSU paper (this dir) |
|---|---|---|
| `bond r0` | 1.529 Å | 1.53 Å (sub-Å rounding; effectively same) |
| `angle k` | 58.35 kcal/mol/rad² | 60 kcal/mol/rad² |
| `angle θ0` | **112.7°** | **109.5°** (sp³ tetrahedral) |
| `dihedral K1..K4` | 2.0058, -2.5128, 0.2710, 6.2900 (OPLS) | **1.73, -4.49, 0.776, 6.99** (Dreiding-UA) |
| `LJ ε, σ` | 0.112, 4.01 | 0.112, 4.01 (match) |

The 112.7° / OPLS dihedral set is what `make_pe_box()` and the
existing CPU PE tests in `tests/cpu/md/` use; it produces
self-consistent NVE / NVT trajectories but does **not** match the
Dreiding-UA potential the MSU paper validates against. Reproducing the
MSU stress-strain curve within the ±5% tolerance required by
`msu-pe-validation.md` Phase D will require **migrating
`polyethylene_tension.cpp` to the Dreiding-UA constants above** and
retuning the embedded tests' tolerances on the new reference numbers.
Tracked as a follow-up; not closed by simply staging these files.
