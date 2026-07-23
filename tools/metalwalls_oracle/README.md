# Offline oracle for the MetalWalls capability loop

Mirrors `tools/pyscf_oracle/` and `tools/cp2k_oracle/`: produce **hard-coded
pins** for regressions under `metalwalls-method-loop`.

**ctest never launches LAMMPS or MetalWalls.**

## Binding oracle priority

| Condition | Primary offline oracle |
| --- | --- |
| Feature/algorithm **also exists in LAMMPS** | **LAMMPS** |
| Feature is **MetalWalls-only** (no LAMMPS twin) | **MetalWalls** |
| Neither binary available for the agent | **Literature** value with full citation; `perf_exempt` for speed |

Examples that **must** use LAMMPS when the binary is available offline:

- `pair_style lj/cut`, `lj/cut/coul/long`, `coul/long`
- `kspace_style pppm` / `ewald` / PME
- `bond_style harmonic`, angles, dihedrals, `special_bonds`
- `fix nve`, `fix nvt` (Nose–Hoover), `fix npt`
- `fix setforce 0` / freeze
- `fix efield`
- `pair_style lj/cut/coul/long` soft-core variants when matching ATV soft-core
- LAMMPS **ELECTRODE** / constant-potential packages **when installed** —
  prefer over MetalWalls for that twin

Examples that use **MetalWalls**:

- Gaussian electrode charge + BO / CG / MaZe as in MW
- Thomas–Fermi electrode metallicity
- PIM + electrode joint polarization
- Piston electrodes, MW 2D electrochemical demos without LAMMPS twin

## Contract

1. Human or agent runs **LAMMPS or MetalWalls offline** (login / batch) with a
   fixed input.
2. Record energy / forces / charges / temperature / capacitance into
   `REFERENCE.md` under the regression (cite code version, units, pair/kspace
   settings, timestep, ensemble).
3. Every `REFERENCE.md` must state:

   ```text
   oracle: lammps | metalwalls | literature
   oracle_version: <string>
   ```

4. Regression binary asserts Atomiverse vs those numbers within stated tol.
5. Optional fair wall: same geometry / comparable settings; write
   \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) into REFERENCE (not ctest).

## Suggested layout (grow as needed)

```text
tools/metalwalls_oracle/
  README.md           ← this file
  inputs/lammps/      ← minimal in.lammps twins (optional)
  inputs/metalwalls/  ← minimal MW runtime/data twins (optional)
  pin_notes.md        ← how pins were taken
  fair_bench.py       ← optional wall twin (later)
```

## MD-only build reminder

Workers in `metalwalls-method-loop` **must not** enable HF/PW/GPW/CUDA. See
`docs/developer/metalwalls-method-loop.md` Build section.

## Relation to capability map

Delivery slices **E1–E12**:
`.claude/notes/metalwalls-capability-map.md`.

Loop product pages: `docs/tutorial/methods/<slug>.md`.
