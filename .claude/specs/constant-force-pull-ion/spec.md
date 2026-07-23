# Spec: `constant-force-pull-ion`

| Field | Value |
| ----- | ----- |
| Slug | `constant-force-pull-ion` |
| Branch | `wf/constant-force-pull-ion` |
| Domain | classical MD / shared electrochem-adjacent (EDL transport twin) |
| Method | Constant-force steered ion pull (noneq MD foothold) |
| Observable | per-atom \(\mathbf{F}^{\mathrm{add}}\), bias energy \(U = -\mathbf{F}\cdot\mathbf{r}^{\mathrm{u}}\); short-MD \(\langle z\rangle\) / work proxy \(W=F_z\Delta z^{\mathrm{u}}\) |
| Geometry | single (or type-selected) mobile ion in a cubic metal-unit box; optional inert background types; no pair force for the pin |
| Oracle | **LAMMPS** `fix addforce` (mandatory; feature exists in LAMMPS). `fix smd` is out of scope for this foothold. |
| Regression | `regressions/md/constant-force-pull-ion/` |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] constant-force-pull-ion` |

## Goal

Ship a public CPU component that adds a **constant external force** to a
selected atom type (the “pulled ion”), matching LAMMPS `fix addforce fx fy fz`
(metal units) on forces and bias energy, and verify short NVE displacement /
work proxy against an offline LAMMPS twin. Enables EDL transport / steered-ion
classics without a new Driver category.

## Non-goals

- LAMMPS `fix smd` (spring / constant-velocity steered MD with CV).
- Variable-style forces (`v_`), region keyword, RESPA-level selection.
- Per-atom variable energy (`energy v_…`) for minimization.
- Virial tally for the added force (can be added later if a pressure twin needs it).
- CUDA path; electrode const-pot; SCF external force.
- Public `*Config` bags or `make_*` / `build_*` factories.

## Public API (binding shape)

```cpp
// Constant force on all atoms of type type_id (matches kAtom_type).
// Force components in eV/Å (LAMMPS metal units).
cpu::md::FixAddForce pull(/*Fx=*/0.0, /*Fy=*/0.0, /*Fz=*/0.05, /*type_id=*/1);
// members may be reassigned: pull.Fz = ...; pull.type_id = ...

// Single-point force/energy:
cpu::Evaluate eval;
eval.on(evaluate::Phase::evaluate).add(std::move(pull));

// MD (default_bindings → post_force):
cpu::md::MD sim;
sim.add(cpu::md::FixAddForce{0.0, 0.0, 0.05, /*type_id=*/1});
```

- **No** `AddForceConfig` / nested `::Config`.
- Parameters are **constructor args** + public members.
- Kind: `cpu::Bias` (external energy/force; satisfies Evaluate / MD energy gate).
- Selection: single atom **type** (same pattern as `FixFreeze`) — maps to a
  LAMMPS type group for textbook ion-pull cases.
- Units (metal / public surface): \(\mathbf{F}\) in **eV/Å**, \(U\) in **eV**,
  positions in **Å**, time in **fs**, mass in **amu**.

## Physics (LAMMPS `fix_addforce.cpp` constant path)

For each owned atom with `type == type_id`:

\[
\mathbf{F}_i \mathrel{+}= \mathbf{F}^{\mathrm{add}},\qquad
U \mathrel{+}= -\mathbf{F}^{\mathrm{add}}\cdot\mathbf{r}_i^{\mathrm{u}}
\]

where \(\mathbf{r}^{\mathrm{u}} = \mathbf{r} + \mathrm{image}\cdot\mathbf{L}\)
(orthorhombic unwrap; triclinic throws in this foothold). Accumulates into
`kAtom_f*` and `meta.get(kEnergy)`.

Work proxy for a short NVE pull of one ion along \(z\):

\[
W = F_z\,\Delta z^{\mathrm{u}} = F_z\bigl(z^{\mathrm{u}}(t)-z^{\mathrm{u}}(0)\bigr)
\]

(equivalent to \(-\Delta U\) of the constant-force bias potential).

## Tasks

1. Write this spec + `acceptance.md` on `wf/constant-force-pull-ion` **before** production code.
2. Offline LAMMPS pin → `regressions/md/constant-force-pull-ion/REFERENCE.md`
   (`oracle: lammps`, `oracle_version: …`) + input under
   `tools/metalwalls_oracle/inputs/lammps/`.
3. Implement `cpu::md::FixAddForce` (`include/…/dynamics/fix_addforce.h`,
   `src/cpu/dynamics/fix_addforce.cpp`); wire CMake under `ATV_BUILD_MD`.
4. Unit test `tests/cpu/dynamics/test_fix_addforce.cpp` (F add, U formula, type filter, image unwrap).
5. E2E regression `reg_constant-force-pull-ion` vs offline pin (dual: force+energy and short-MD \(z\) / work).
6. MD-only build + run until PASS.
7. Perf gate vs same LAMMPS binary (or `perf_exempt` with reason).
8. Commit on `wf/constant-force-pull-ion` only.

## Acceptance criteria

See `acceptance.md` (binding).
