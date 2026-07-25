# Acceptance — pw-uks-fe-smearing

## Correctness

- [ ] `PW::run` accepts `spin_polarized && smearing != none` (no throw).
- [ ] Restricted smear path (`pw-smearing-al` settings) and fixed-occ UKS
      (`pw-spin-polarized-o2` / `pw-uks-kmesh-fe` knobs) remain green if built
      (no intentional bit-break of existing paths).
- [ ] Regression `reg_pw-uks-fe-smearing`:
  - SCF `converged == true`
  - `|n_elec_occ - Z_val|` within REFERENCE tol (\(Z_{\mathrm{val}}=16\))
  - free energy within REFERENCE pin tol
  - magnetization finite and ferro-signed as specified in REFERENCE
  - finite `e_fermi`
- [ ] Unit tests for spin-orbital Fermi / occupations added or extended under
      `tests/cpu/planewave/test_pw_smearing.cpp` (or sibling).
- [ ] Offline pin documented in `regressions/pw/pw-uks-fe-smearing/REFERENCE.md`
      with settings, units, and CP2K algorithm twin note.
- [ ] No public `*Config` / `make_*` / `build_*` factories.
- [ ] Public energy path still writes eV to Frame when used via evaluator;
      internal Ha/Ry cited on `PWResult`.

## Perf

- [ ] Wall \(T_{\mathrm{atv}} \le 60\) s for the regression binary (OMP default).
- [ ] `within_10x` vs CP2K fair twin **or** `perf_exempt` with reason
      (expected: no fair NC-UPF PW twin wall — same class as
      `pw-uks-kmesh-fe` / `pw-smearing-al` energy twin notes).

## Process

- [ ] Spec under `.claude/specs/pw-uks-fe-smearing/`
- [ ] Branch `wf/pw-uks-fe-smearing`, commit with clear message
- [ ] Scoped build only: `atomiverse`, touched `test_*`, `reg_pw-uks-fe-smearing`
- [ ] Build trees wiped after dual gate
