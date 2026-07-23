# Spec: `fd-phonon-supercell`

| Field | Value |
| ----- | ----- |
| Slug | `fd-phonon-supercell` |
| Branch | `wf/fd-phonon-supercell` |
| Domain | response (periodic FD phonons) |
| Method | Finite-displacement phonons on a supercell |
| Observable | Γ-point phonon wavenumbers \(\omega\) [cm⁻¹] |
| Geometry | Diamond Si (conventional or primitive + supercell) |
| Oracle | CP2K `VIBRATIONAL_ANALYSIS` when binary available offline; else literature/physical pins (acoustic sum rule + SW Si optical band) |
| Regression | `regressions/pw/fd-phonon-supercell/` |
| Queue row | `docs/developer/cp2k-method-queue.md` `## in_progress` `[~] fd-phonon-supercell` |
| Enables | `dfpt-phonon-si` |

## Goal

Ship a **public** finite-displacement **crystal phonon** path that is *not*
molecular `cpu::vib::Hessian` (which projects 3T+3R free-molecule modes and has
no supercell). This foothold:

1. Builds an \(n_a \times n_b \times n_c\) supercell from a unit-cell `Frame`.
2. Assembles Cartesian force constants by central differences of forces.
3. Diagonalizes the **mass-weighted dynamical matrix at Γ** (all \(3N_\mathrm{sc}\)
   modes; acoustic branches numerically ~0 at equilibrium).
4. Returns wavenumbers in **cm⁻¹** on a plain result type (no Frame pollution).

Force engines are **user-bound** (classical `Potential`, PW, GPW, …) on
`vib::Phase::force` — same composition rule as `cpu::vib::Hessian` / AIMD.

## Non-goals

- DFPT / analytic dynamical matrices (`dfpt-phonon-si`).
- Finite-\(q\) band structure beyond the supercell fold (phonopy-style path).
- Acoustic sum-rule projection as a required step (full \(3N\) spectrum is the
  product surface; acoustic near-zero is an acceptance check).
- Molecular thermochemistry / IR / free-molecule TR projection (stays on
  `cpu::vib::Hessian`).
- Public `*Config` bags, `make_*` / `build_*` factories, CUDA path.

## Public API (binding shape)

```cpp
// Expand unit cell → supercell Frame (rows of simbox.cell = lattice vectors).
cpu::Frame sc = cpu::vib::expand_supercell(unit, /*na=*/2, /*nb=*/2, /*nc=*/2);

// Driver: parameters = ctor args / public members (no PhononConfig bag).
cpu::vib::Phonon phonon(/*na=*/1, /*nb=*/1, /*nc=*/1);  // Γ of the given cell
// or: phonon.na = 2; phonon.nb = 2; phonon.nc = 2;

// Classical force twin (SW Si):
phonon.set_neighbor(cpu::md::NeighborList<>{cutoff, skin, 1, 0, true});
phonon.on(ATV::vib::Phase::force)
    .add(cpu::StillingerWeber{cpu::SWPotential{molrs::Element::Si, "stillinger-weber"}});

// Plane-wave (future e2e / CP2K twin):
// PW method; method.k_mesh = 1; ...
// phonon.on(vib::Phase::force).add(cpu::pw::PWEvaluator{std::move(method)}.with_gradient(true));

phonon.build();
cpu::vib::PhononResult r = phonon.run(unit_cell_frame);
// r.frequencies [cm^-1], length 3 * N_supercell, ascending
// r.n_atoms, r.na/nb/nc, r.hessian [eV/Å²]
```

Rules:

- **No** `PhononConfig` / nested `::Config`.
- Parameters: constructor `(na, nb, nc)` and public members `na`, `nb`, `nc`.
- Kind: preset over `Driver<vib::Phase>` (same phase enum as Hessian).
- Units: public \(\omega\) in **cm⁻¹**; force constants in **eV/Å²** on the result;
  forces on Frame remain eV/Å.

## Physics

### Supercell

Given unit-cell lattice vectors \(\mathbf{a},\mathbf{b},\mathbf{c}\) (rows of
`frame.simbox.cell()`):

\[
\mathbf{a}' = n_a\mathbf{a},\quad
\mathbf{b}' = n_b\mathbf{b},\quad
\mathbf{c}' = n_c\mathbf{c},
\]

positions of image \((i,j,k)\) are \(\mathbf{R} + i\mathbf{a} + j\mathbf{b} + k\mathbf{c}\)
for \(i\in[0,n_a)\), etc. Types/masses replicate; tags renumber uniquely.

### Force constants (Γ of the supercell)

Central difference with \(h = 10^{-3}\) Å (same as `kHessianStepAngstrom`):

\[
\Phi_{i\alpha,j\beta}
  = -\frac{F_{j\beta}(+\mathbf{e}_{i\alpha}h) - F_{j\beta}(-\mathbf{e}_{i\alpha}h)}{2h}
  \quad[\mathrm{eV}/\mathrm{\AA}^2].
\]

Symmetrize \(\Phi \leftarrow (\Phi+\Phi^T)/2\).

### Dynamical matrix at Γ

Mass-weight in Hartree/Bohr²/amu (reuse `MassWeighting` unit path), then

\[
D_{i\alpha,j\beta}
  = \frac{\Phi_{i\alpha,j\beta}}{\sqrt{m_i m_j}},
\]

diagonalize \(D\), convert eigenvalues to wavenumbers via `nist::au2wavenumber()`
(same CODATA-2014 path as molecular vib). Report **all \(3N\)** modes (imaginary
as negative). At a stress-free crystal minimum the three lowest \(|\omega|\) are
acoustic and must be near zero.

## Difference from `cpu::vib::Hessian`

| | Molecular Hessian | Phonon (this spec) |
| - | --- | --- |
| System | free molecule | crystal (PBC forces) |
| Cell | ignored for modes | supercell expansion |
| Diagonalizer | TR project \(3N-6/5\) | full \(3N\) crystal spectrum |
| Output | vib + thermochem | \(\omega\) [cm⁻¹] only |

## Regression twin

- **Geometry:** diamond Si, conventional cubic \(a=5.431\) Å, 8 atoms (1×1×1
  supercell of the conventional cell ≡ Γ of that cell).
- **Force engine (offline pin host, no CP2K binary this shift):** Stillinger–Weber
  Si (`cpu::StillingerWeber`, parameters from Stillinger & Weber, *Phys. Rev. B*
  **31**, 5262 (1985) — already in-tree). Documents FD machinery with a periodic
  many-body force; CP2K Quickstep twin input skeleton may live under
  `tools/cp2k_oracle/inputs/` for later wall timing.
- **Pins:**
  1. `n_modes == 24`.
  2. Acoustic sum rule: three lowest \(|\omega| < 30\) cm⁻¹.
  3. Highest optical \(\omega\) in the published SW-Si optical band
     **[500, 900]** cm⁻¹ (literature SW overestimates experiment ~520 cm⁻¹;
     see REFERENCE.md).

## Files (expected)

| Path | Role |
| ---- | ---- |
| `include/atomiverse/cpu/vibration/supercell.h` | `expand_supercell` |
| `include/atomiverse/cpu/vibration/phonon_result.h` | result type |
| `include/atomiverse/cpu/vibration/phonon_diagonalizer.h` | crystal diagonalizer |
| `include/atomiverse/cpu/vibration/phonon.h` | `cpu::vib::Phonon` preset |
| `src/cpu/vibration/phonon.cpp` | wiring + `run` |
| `src/cpu/vibration/phonon_diagonalizer.cpp` | \(3N\) eigensolve |
| `src/cpu/vibration/supercell.cpp` | expansion impl |
| `tests/cpu/vibration/test_phonon_supercell.cpp` | unit |
| `regressions/pw/fd-phonon-supercell/` | e2e + REFERENCE.md |
