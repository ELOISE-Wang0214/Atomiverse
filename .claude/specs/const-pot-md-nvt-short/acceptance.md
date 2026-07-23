# Acceptance — `const-pot-md-nvt-short`

Binding checklist. `acceptance_met=true` only when **every** box is verified.

## Oracle and pin

- [ ] `oracle: metalwalls` (MW-only physics path; LAMMPS ELECTRODE not available
      offline on the worker binary — installed packages lack ELECTRODE)
- [ ] `oracle_version:` recorded in
      `regressions/md/const-pot-md-nvt-short/REFERENCE.md`
- [ ] Offline pin documents: geometry, \(\Psi\), Gaussian width \(w\), algorithm
      (Siepmann–Sprik / Reed BO + electroneutrality), continuum capacitor
      \(Q\) sign/order, smoke MD criteria
- [ ] No live MetalWalls / LAMMPS launch at ctest time

## Observable and tolerance

- [ ] **Geometry:** parallel-plate cell — left wall type 2, right wall type 3,
      Na+/Cl− electrolyte (types 0/1), orthorhombic box; walls frozen
- [ ] **Static BO (before / early MD):** with \(\Psi_L < \Psi_R\),
      \(Q_L < 0 < Q_R\) (higher-potential plate carries positive charge)
- [ ] **Electroneutrality:** \(|Q_L + Q_R| \le 1\times 10^{-8}\) e (electrolyte
      charge-neutral) after each solve in the regression
- [ ] **MD smoke (NVT, ≥ 200 steps, dt ≤ 1 fs):** all of
      `kEnergy`, `kKinetic`, every `kAtom_charge` finite (no NaN/Inf)
- [ ] **Freeze:** wall atom positions match build-time pin within \(10^{-10}\) Å;
      wall \(|v|\) components ≤ \(10^{-12}\) Å/fs
- [ ] **Q sign stable:** \(\mathrm{sign}(Q_L)<0\) and \(\mathrm{sign}(Q_R)>0\)
      at end of short trajectory (same \(\Psi\) polarity)

## Architecture

- [ ] `cpu::ConstPotentialBO` is a Component on `cpu::md::MD` (pre_force);
      private \(A,Q\); writes `kAtom_charge`; **not** inside `KSpace*`
- [ ] `cpu::md::FixFreeze` freezes by atom type; ConstraintProvider DOF
- [ ] No public `*Config` bag; no `make_*` / `build_*` factory
- [ ] Frame I/O only

## Build / tests

- [ ] MD-only configure:
      `ATV_BUILD_MD=ON`, `ATV_BUILD_KSPACE=ON`, `ATV_BUILD_TESTS=ON`,
      `ATV_BUILD_REGRESSIONS=ON`,
      `ATV_BUILD_HF=OFF`, `ATV_BUILD_PW=OFF`, `ATV_BUILD_GPW=OFF`,
      `ATV_BUILD_CUDA=OFF`, `ATV_BUILD_TORCH=OFF`
- [ ] Unit tests under `tests/` pass (FixFreeze and/or ConstPotentialBO)
- [ ] Regression `reg_const-pot-md-nvt-short` passes; ≤ 60 s single-rank

## Dual-gate (correctness + perf)

- [ ] Correctness: regression PASS vs offline pin / acceptance bounds
- [ ] Perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\) **or**
      `perf_exempt=true` with reason (no MetalWalls / ELECTRODE twin binary)

## Delivery

- [ ] Branch `wf/const-pot-md-nvt-short`; never commit on `dev`; never push
- [ ] Commit includes spec, components, unit tests, regression, REFERENCE, CMake
