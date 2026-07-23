# Spec: `external-efield-md`

| Field | Value |
| ----- | ----- |
| Slug | `external-efield-md` |
| Branch | `wf/external-efield-md` |
| Domain | classical MD / electrochemistry (map **E3**) |
| Method | Applied uniform electric field on MD charges |
| Observable | \(F_i = q_i \mathbf{E}\), field energy \(U = -\sum_i q_i \mathbf{E}\cdot\mathbf{r}_i^{\mathrm{u}}\) |
| Geometry | 3 point charges (ions) in a cubic metal-unit box; no pair force for the pin |
| Oracle | **LAMMPS** `fix efield` (mandatory; feature exists in LAMMPS) |
| Regression | `regressions/md/external-efield-md/` |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] external-efield-md` |

## Goal

Ship a public CPU component that applies a **uniform external electric field** to
charged atoms in classical MD, matching LAMMPS `fix efield` (metal units) on
forces and field energy. Enables E3 / finite-field electrode setups later; this
slug is the **independent** E-field foothold (no const-pot / BO dependency).

## Non-goals

- Finite-field + single electrode (queue `finite-field-single-electrode`).
- Time-/space-varying or atom-style variable fields (LAMMPS `v_` forms).
- Dipole / torque (`atom_style` dipole) path.
- Region-restricted efield.
- CUDA path, SCF external field (already separate), const-pot electrodes (E1).
- Public `*Config` bags or `make_*` / `build_*` factories.

## Public API (binding shape)

```cpp
cpu::md::FixEfield efield(/*Ex=*/0.1, /*Ey=*/-0.05, /*Ez=*/0.2); // [V/Å]
// members may be reassigned: efield.Ex = ...
// Single-point:
cpu::Evaluate eval;
eval.on(evaluate::Phase::evaluate).add(std::move(efield));
// MD:
cpu::md::MD sim;
sim.add(cpu::md::FixEfield{Ex, Ey, Ez});  // default_bindings → post_force
```

- **No** `EfieldConfig` / nested `::Config`.
- Parameters are **constructor args** + public members.
- Kind: `cpu::Bias` (external energy/force; satisfies Evaluate / MD energy gate).
- Units (metal / public surface): \(\mathbf{E}\) in **V/Å**, \(q\) in \(e\),
  \(\mathbf{F}\) in **eV/Å**, \(U\) in **eV**. With LAMMPS metal `qe2f = 1`,
  \(\mathbf{F} = q\mathbf{E}\) and \(U = -\sum q\,\mathbf{E}\cdot\mathbf{r}^{\mathrm{u}}\).

## Physics (LAMMPS `fix_efield.cpp` constant field)

For each owned atom with charge \(q_i\):

\[
\mathbf{F}_i \mathrel{+}= q_i \mathbf{E},\qquad
U \mathrel{+}= -\mathbf{F}_i\cdot\mathbf{r}_i^{\mathrm{u}}
= -q_i\,\mathbf{E}\cdot\mathbf{r}_i^{\mathrm{u}}
\]

where \(\mathbf{r}^{\mathrm{u}}\) is the **unwrapped** coordinate
(\(\mathbf{r} + \mathrm{image}\cdot\mathbf{L}\) for orthorhombic). Writes
accumulate into `kAtom_f*` and `meta.get(kEnergy)`.

## Tasks

1. Write this spec + `acceptance.md` on `wf/external-efield-md` **before** production code.
2. Offline LAMMPS pin → `regressions/md/external-efield-md/REFERENCE.md`
   (`oracle: lammps`, `oracle_version: …`) + optional input under
   `tools/metalwalls_oracle/inputs/lammps/`.
3. Implement `cpu::md::FixEfield` (`include/…/dynamics/fix_efield.h`,
   `src/cpu/dynamics/fix_efield.cpp`); wire CMake under `ATV_BUILD_MD`.
4. Unit test `tests/cpu/dynamics/test_fix_efield.cpp` (F=qE, U formula).
5. E2E regression `reg_external-efield-md` vs offline pin (dual force+energy).
6. MD-only build + run until PASS.
7. Perf gate vs same LAMMPS binary (or `perf_exempt` with reason).
8. Commit on `wf/external-efield-md` only.

## Acceptance criteria

See `acceptance.md` (binding).
