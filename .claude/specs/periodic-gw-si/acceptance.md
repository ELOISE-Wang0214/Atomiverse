# Acceptance: `periodic-gw-si`

## Dual gate

1. **Unit:** `test_cpu_gpw_g0w0`
   - dRPA amplitudes path: energy matches `direct_rpa_rhf_energy` to 1e-12 Ha
     and \(X^\top X - Y^\top Y \approx I\) (norm check).
   - `GpwDft::g0w0` guard throws (before build / unconverged / UKS /
     multi-k) without running SCF.
2. **Regression:** `reg_periodic-gw-si` (`regressions/gpw/periodic-gw-si/`)
   - SCF converges; `gw.converged == true`.
   - QP gap finite and **positive**; HOMO QP below LUMO QP.
   - \(|\varepsilon_p^{\mathrm{QP}}-\varepsilon_p^{\mathrm{KS}}|\) finite for
     HOMO and LUMO.
   - Frozen pin band on **QP gap [eV]** (offline CP2K when a fair spectral
     twin exists; else first-green setup-integrity band) — see REFERENCE.md.
   - Wall ≤ 60 s single process.

## Perf

`T_atv / T_oracle ≤ 10` when a fair CP2K wall exists; else `perf_exempt`
with reason recorded in REFERENCE.md (CP2K RI-GW frequency path ≠ spectral
plasmon G0W0; MT exchange vs periodic HF).

## Done when

- Spec closed (delete this directory on merge).
- `acceptance_met=true` in worker schema.
- Commit on `wf/periodic-gw-si`.
