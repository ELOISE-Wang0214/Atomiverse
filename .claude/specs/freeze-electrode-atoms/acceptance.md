# Acceptance — `freeze-electrode-atoms`

Binding checklist. `acceptance_met=true` only when **every** box is verified.

## Oracle and pin

- [ ] `oracle: lammps` (feature exists as `fix setforce 0` / freeze path)
- [ ] `oracle_version:` recorded in `REFERENCE.md`
- [ ] Offline pin path: `tools/metalwalls_oracle/inputs/lammps/freeze-electrode-atoms/`
      (or documented twin input under that tree) + numbers in
      `regressions/md/freeze-electrode-atoms/REFERENCE.md`
- [ ] No live LAMMPS / MetalWalls launch at ctest time

## Observable and tolerance

- [ ] **Observable:** frozen-set Cartesian positions \(x,y,z\) [Å] and velocities
      \(v_x,v_y,v_z\) [Å/fs] after NVE MD with mobile fluid
- [ ] **Geometry:** 4 wall (type 2) + 4 fluid (type 1), box 10 Å cubic, LJ
      \(\varepsilon=0.010341\) eV, \(\sigma=3.4\) Å, cutoff 8 Å, \(dt=1\) fs,
      1000 steps (matches offline pin)
- [ ] Wall atom positions match initial lattice pin within **1e-10 Å** (component)
- [ ] Wall velocities magnitude max component ≤ **1e-12** Å/fs (or 0 exactly)
- [ ] Fluid atoms are **not** frozen (at least one fluid atom displaces >
      1e-3 Å) — dual check that freeze is selective

## Architecture

- [ ] Public type is a **Component** usable on `cpu::md::MD` integrate path
- [ ] No public `*Config` bag; no `make_*` / `build_*` factory
- [ ] Frame I/O only (no parallel `*System`)
- [ ] Parameters as ctor args / members (`type` id of frozen atoms)

## Build / tests

- [ ] MD-only configure: `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`,
      `ATV_BUILD_TESTS=ON`, `ATV_BUILD_REGRESSIONS=ON`,
      `ATV_BUILD_HF=OFF`, `ATV_BUILD_PW=OFF`, `ATV_BUILD_GPW=OFF`,
      `ATV_BUILD_CUDA=OFF`, `ATV_BUILD_TORCH=OFF`
- [ ] Unit test under `tests/` passes
- [ ] Regression `reg_freeze-electrode-atoms` passes vs offline pin

## Dual-gate (correctness + perf)

- [ ] Correctness: regression PASS
- [ ] Perf: fair wall ratio \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\) ≤ 10.0
      **or** `perf_exempt=true` with reason in REFERENCE when oracle binary
      unavailable for fair twin

## Delivery

- [ ] Branch `wf/freeze-electrode-atoms`; never commit on `dev`; never push
- [ ] Commit includes spec, component, unit test, regression, REFERENCE, CMake wire
