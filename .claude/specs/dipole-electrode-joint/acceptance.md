# Acceptance — dipole-electrode-joint

All criteria must pass for `acceptance_met=true`.

## Build

- [ ] MD-only configure: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`,
      `ATV_BUILD_TESTS=ON`, `ATV_BUILD_REGRESSIONS=ON`,
      **HF/PW/GPW/CUDA/Torch OFF**, `BLA_VENDOR=OpenBLAS`
- [ ] `atomiverse` + unit test + `reg_dipole-electrode-joint` link cleanly

## Architecture

- [ ] Public type is a **Component** on `cpu::md::MD` (default `pre_force`);
      **not** a new Driver category
- [ ] No public `*Config` bag; no `make_*` / `build_*` factory
- [ ] Charge / dipole joint solve is **not** inside `KSpace*`
- [ ] Frame I/O only (`kAtom_charge` write for electrode \(Q\); geometry from Frame)

## Physics / unit

- [ ] Unit tests under `tests/` only (single-function / method)
- [ ] \(\alpha=0\) reduction matches dense `ConstPotBO` \(Q\) on same geometry
      within \(10^{-10}\,e\)
- [ ] Joint residual \(\max|\Delta Q|\), \(\max|\Delta\mu|\) below `joint_tol`
      after convergence
- [ ] Electroneutrality residual \(|\sum Q_i| \le 10^{-10}\,e\) when enabled
- [ ] Coupling check: \(Q(\alpha>0) \neq Q(\alpha=0)\) (relative \(|\Delta Q|/|Q|
      \ge 10^{-6}\) on pin geometry)

## Regression (e2e)

- [ ] One case: `regressions/md/dipole-electrode-joint/`
- [ ] `REFERENCE.md` states `oracle: metalwalls` and `oracle_version: …`
- [ ] Offline pin (no live MW/LAMMPS at ctest)
- [ ] Dual gate (both must hold):
  1. Electrode charges: relative \(|Q_{\mathrm{top}}-Q_{\mathrm{pin}}|
     /|Q_{\mathrm{pin}}| \le 10^{-6}\) (and bot consistent with neutrality)
  2. Electrolyte dipoles: \(\max|\mu-\mu_{\mathrm{pin}}| \le 10^{-8}\) e·Å
     (or relative \(\le 10^{-6}\) if \(|\mu|\) large)
- [ ] Exit 0 on pass; non-zero on fail; ≤ 60 s

## Perf

- [ ] Fair wall vs MetalWalls if binary present; else
      `perf_exempt=true` with reason recorded in REFERENCE.md
- [ ] `status=green` requires `within_10x=true` **or** `perf_exempt=true`

## Git

- [ ] Commit on `wf/dipole-electrode-joint` only (never `dev`; never push)
