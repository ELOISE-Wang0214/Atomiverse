# Spec: `freeze-electrode-atoms`

| Field | Value |
| ----- | ----- |
| Slug | `freeze-electrode-atoms` |
| Branch | `wf/freeze-electrode-atoms` |
| Domain | classical MD (E1 electrode wall immobility) |
| Map slice | E1 — Freeze / immobile electrode atoms (capability map §1 P0) |
| Observable | frozen-set positions fixed; velocities zero after MD |
| Geometry | 4 wall atoms (type 2) + 4 fluid atoms (type 1), LJ, cubic 10 Å |
| Public API | `cpu::md::FixFreeze` Component on integrate path (no `*Config`, no factory) |
| Regression | `regressions/md/freeze-electrode-atoms/` |
| Oracle | **LAMMPS** `fix setforce 0.0 0.0 0.0` (mandatory — same algorithm) |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] freeze-electrode-atoms` |

## Goal

Land a **Component** that freezes a selected atom subset (electrode / wall
atoms) under classical MD so that positions stay at their build-time values and
velocities remain zero while mobile fluid atoms integrate normally. This is an
E1 dependency for constant-potential electrode MD (walls are almost always
immobile).

Architecture (binding):

- Component on existing `cpu::md::MD` — **not** a new Driver category.
- Parameters = constructor arguments / members only — **no** public `*Config`.
- Selection by atom `type` (electrode kind id), matching LAMMPS `group wall type 2`.
- LAMMPS twin: `fix setforce 0.0 0.0 0.0` on the wall group + initial
  `velocity wall set 0 0 0`.

## Non-goals

- Constant-potential charge solve, Gaussian electrodes, piston electrodes.
- CUDA / multi-GPU freeze.
- Generic group algebra beyond single-type selection.
- SCF / PW / GPW / Torch.
- Live LAMMPS at ctest time.

## Tasks

1. **Spec + acceptance** (this file + `acceptance.md`) before production code.
2. **Offline LAMMPS pin** under `tools/metalwalls_oracle/inputs/lammps/` and
   record numbers in `regressions/md/freeze-electrode-atoms/REFERENCE.md`
   (`oracle: lammps`, `oracle_version:`).
3. **Implement** `cpu::md::FixFreeze`:
   - Header `include/atomiverse/cpu/dynamics/fix_freeze.h`
   - Impl `src/cpu/dynamics/fix_freeze.cpp`
   - Wire into `CMakeLists.txt` (`ATV_BUILD_MD` sources)
   - Zero forces on frozen type at `post_force` (after physical + Langevin)
   - Zero velocities and restore positions after integrate halves
   - Report `3 * n_frozen` via `ConstraintProvider` for DOF accounting
4. **Unit test** `tests/cpu/dynamics/test_fix_freeze.cpp` — one method scope
   (launch zeros v/f; position restore).
5. **E2E regression** `regressions/md/freeze-electrode-atoms/` public API only;
   assert frozen \(x\) pin + \(v=0\) vs LAMMPS offline values; dual-gate notes.
6. **MD-only build** (HF/PW/GPW/CUDA/Torch OFF) + run unit + regression.
7. **Perf gate**: fair wall vs same LAMMPS input; record ratio in REFERENCE.
8. Commit on `wf/freeze-electrode-atoms` only.

## Acceptance criteria

See `acceptance.md` (binding).
