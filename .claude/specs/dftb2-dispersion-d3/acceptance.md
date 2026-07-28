# Acceptance: dftb2-dispersion-d3

## Dual gate

### Unit (`tests/cpu/semiempirical/test_dftb2_dispersion_d3.cpp`)

- `d3=false`: `last_e_disp() == 0`, energy finite.
- `d3=true` on water dimer: `last_e_disp() < 0` (attractive), finite.
- Difference of total energies with/without `d3` equals `last_e_disp()` within 1e-9 eV.

### Regression (`regressions/semiempirical/dftb2-dispersion-d3/`)

- Hard-coded \(E_{\mathrm{disp}}\) [eV] matches offline pin in `REFERENCE.md`
  within absolute tol (≤ 1e-4 eV).
- \(E_{\mathrm{tot}}(d3)-E_{\mathrm{tot}}(off)=E_{\mathrm{disp}}\) within 1e-8 eV.
- Finite non-zero total energy with `d3=true`.
- Exit 0; public headers + link `atomiverse` only; ≤ 60 s.

### Build / process

- SE-only cmake (SEMIEMPIRICAL+EVALUATE+TESTS+REGRESSIONS ON; HF/MD/CUDA/… OFF).
- Scoped build: `--target atomiverse test_cpu_dftb2_dispersion_d3 reg_dftb2-dispersion-d3`.
- Scoped test: `ctest -R '^(test_cpu_dftb2_dispersion_d3)$'` + `run_reg_dftb2-dispersion-d3`.
- No public `*Config` / `make_*` / `build_*`.
- Spec present under `.claude/specs/dftb2-dispersion-d3/` until merge closes it.

## Done when

All of the above green on branch `wf/dftb2-dispersion-d3`; commit on that branch;
`acceptance_met=true`.
