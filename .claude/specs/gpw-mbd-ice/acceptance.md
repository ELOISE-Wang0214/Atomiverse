# Acceptance — `gpw-mbd-ice`

## Dual gate

1. **Unit** (`test_cpu_gpw_mbd`):
   - Isolated free-atom MBD@rsSCS water dimer E_disp matches offline libMBD
     Fortran `cmbd_mbd_scs_energy` / pure-Python `mbd_energy` to **1e-12 Ha**.
   - `GpwDft::mbd` folds `e_disp` into `e_total` (short SCF).

2. **Regression** (`reg_gpw-mbd-ice`):
   - Same geometry as `gpw-ts-ice` (periodic 12 Å cell for SCF; MBD pin is
     **isolated free-atom** MBD@rsSCS — geometry-only, basis-independent).
   - Hard gate: `|E_disp − E_ref| ≤ 1e-8 Ha`.
   - E_ref = **−1.4928804171407961×10⁻³ Ha** (libMBD Fortran, ATV Å→Bohr).

## Perf

`perf_exempt = true` — E_disp pin; full QS wall unfair (same as TS / D3BJ / D4).

## Done when

- Spec acceptance met; unit + `run_reg_gpw-mbd-ice` green.
- Commit on `wf/gpw-mbd-ice`.
