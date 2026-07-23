# Spec: `const-pot-md-nvt-short`

| Field | Value |
| ----- | ----- |
| Slug | `const-pot-md-nvt-short` |
| Branch | `wf/const-pot-md-nvt-short` |
| Domain | classical electrochemistry (E1 product) |
| Method | Constant-potential BO + freeze + NVT MD short trajectory |
| Chapter | E1 (capability map §1 P0 product) |
| Geometry | NaCl electrolyte + two parallel planar walls |
| Observable | finite \(E\), electrode \(Q\) sign vs \(\Psi\), no NaN |
| Tolerance | smoke stability + offline pin (see acceptance) |
| Public API | `cpu::ConstPotentialBO` + `cpu::md::FixFreeze` Components on `cpu::md::MD` |
| Regression | `regressions/md/const-pot-md-nvt-short/` |
| Oracle | **metalwalls** algorithm (Siepmann–Sprik / Reed / Scalfi BO); LAMMPS ELECTRODE not installed offline |
| Queue row | `docs/developer/metalwalls-method-queue.md` `## in_progress` `[~] const-pot-md-nvt-short` |
| Map | `.claude/notes/metalwalls-capability-map.md` E1 product row |
| Deps | BO dense inverse, freeze walls, `CoulLong` + `KSpacePME` |

## Goal

Land the **E1 integration product**: short NVT MD of a parallel-plate cell with
fluctuating electrode charges held at fixed potentials, frozen walls, and
classical real/reciprocal Coulomb for the electrolyte + electrode point-charge
force path. Proves the composition on `cpu::md::MD` (no new Driver category):

1. `cpu::ConstPotentialBO` solves \(A Q = \Psi - \Phi_{\mathrm{elec}}\) each
   step (dense BO inverse / Cholesky with electroneutrality), writes
   `kAtom_charge` on electrode atoms.
2. `cpu::md::FixFreeze` freezes wall atom types (positions fixed, \(v=0\)).
3. `CoulLong` + `KSpacePME` + soft LJ + `VelocityVerlet` + `NoseHoover` run a
   short NVT trajectory without NaN and with electrode \(Q\) polarity matching
   the applied \(\Delta\Psi\).

## Non-goals

- CG electrode solve, MaZe, charge predictor, Thomas–Fermi, piston.
- Gaussian charge **force** kernels (`gaussian-electrode-charge` separate;
  BO matrix uses Gaussian self/pair form; force path remains point-charge).
- 2D / slab Ewald (`slab-2d-ewald-foothold`).
- Dedicated charge observer / charging-work energy channel (sibling slugs).
- CUDA / multi-GPU; live MetalWalls or LAMMPS at ctest.
- Public `*Config` bags or `make_*` / `build_*` factories.
- Welding charge solve into `KSpace*`.

## Architecture (binding)

- **Components on `cpu::md::MD`** only — no electrode Driver category.
- Parameters = constructor arguments / members.
- Private workspace: electrode index map, dense \(A\), Cholesky factor, \(Q\).
- Bind `ConstPotentialBO` at `md::Phase::pre_force` (after Zero\*, before force).
- Bind `FixFreeze` at integrate + `post_force` (zero \(f\), restore \(x\), zero \(v\)).
- Frame is sole I/O; charges on `kAtom_charge`.

## Tasks

1. Write this spec + `acceptance.md` on `wf/const-pot-md-nvt-short`.
2. Offline pin in `REFERENCE.md` (`oracle: metalwalls` algorithm / literature
   BO; `oracle_version`; continuum capacitor check for sign + order of \(Q\);
   no live binary at ctest). Optional twin input under
   `tools/metalwalls_oracle/inputs/metalwalls/const-pot-md-nvt-short/`.
3. Implement `cpu::md::FixFreeze` and `cpu::ConstPotentialBO`; wire CMake MD
   sources; unit tests under `tests/cpu/dynamics/`.
4. E2E `reg_const-pot-md-nvt-short`: public headers only; NaCl + parallel walls;
   short NVT; assert finite \(E\), \(Q\) polarity, neutrality, freeze, no NaN.
5. MD-only configure/build; run unit + regression until PASS.
6. Perf gate: fair wall vs MetalWalls/LAMMPS ELECTRODE if available; else
   `perf_exempt` with reason.
7. Commit on `wf/const-pot-md-nvt-short` only.

## Acceptance criteria

See `acceptance.md` (binding).
