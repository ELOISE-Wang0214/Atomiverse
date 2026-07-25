# Acceptance — `gpw-mbd-hirshfeld`

## Dual gate

1. **Unit** (`test_cpu_gpw_hirshfeld`):
   - Free-atom **superposition** density \(\rho=\sum_A\rho_A^{\mathrm{free}}\)
     yields \(v_A\approx 1\) (to 1e-3) for the water-dimer geometry on a coarse grid.
   - Isolated MBD@rsSCS with a non-unit `vol_ratio` vector differs from free-atom
     and matches a dual recomputation to 1e-12 Ha.
   - `GpwDft::{mbd,hirshfeld}` reports finite attractive `e_disp` and
     `volume_ratios.size()==n_atoms` with every \(v_A\in(0.3,1.5)\).

2. **Regression** (`reg_gpw-mbd-hirshfeld`):
   - Same water-dimer geometry as `gpw-mbd-ice` / `gpw-ts-ice`.
   - Hard gate on **E_disp** (Hirshfeld-scaled MBD@rsSCS) to **1e-8 Ha** vs
     offline pin in `REFERENCE.md`.
   - Volume ratios reported and within physical bounds.

## Perf

`perf_exempt = true` — E_disp + volume-ratio pin; full QS wall unfair for a
dispersion correction (same policy as `gpw-mbd-ice` / `gpw-ts-ice`).

## Done when

- Spec acceptance met; unit + `run_reg_gpw-mbd-hirshfeld` green.
- Commit on `wf/gpw-mbd-hirshfeld`.
