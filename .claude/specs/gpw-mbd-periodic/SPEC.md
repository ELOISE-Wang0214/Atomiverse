# Spec — `gpw-mbd-periodic`

**Slug:** `gpw-mbd-periodic`  
**Branch:** `wf/gpw-mbd-periodic`  
**Method:** Periodic Ewald + k-mesh MBD@rsSCS  
**Loop:** CP2K capability (GPW/PW/periodic)  
**Oracle:** offline pure-Python pymbd (libMBD algorithm; CP2K stock has no MBD)

## Status: acceptance_met = true (2026-07-26)

Unit `test_cpu_gpw_mbd` PASS; regression `run_reg_gpw-mbd-periodic` PASS
(E_disp = −1.49657207049714×10⁻³ Ha vs pin −1.49657207049891×10⁻³, tol 1e-8).
**perf_exempt** (E_disp pin; same policy as free-atom MBD).

## Delivered

- Kernel `mbd_rsscs_pbe(..., periodic=true, k_mesh)` — Ewald + k-mesh ZPE.
- `GpwDft::mbd_k_mesh` / `GPWEvaluator::mbd_k_mesh`: `{0,0,0}` free-atom;
  non-zero → periodic (requires `!nonperiodic`).
- Regression `regressions/gpw/gpw-mbd-periodic/` + oracle note under
  `tools/cp2k_oracle/inputs/gpw-mbd-periodic/`.
