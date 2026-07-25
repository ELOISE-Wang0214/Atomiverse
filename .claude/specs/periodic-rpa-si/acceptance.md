# Acceptance: `periodic-rpa-si`

## Dual gate

1. **Unit:** `test_cpu_gpw_rpa`
   - Plasmon formula on synthetic `(ia|jb)` matches an independent dense
     A/B eig path to 1e-12 Ha.
   - `GpwDft::direct_rpa` guard throws (before build / unconverged /
     UKS / multi-k) without running SCF.
2. **Regression:** `reg_periodic-rpa-si` (`regressions/gpw/periodic-rpa-si/`)
   - SCF converges; `E_c` finite and **strictly negative**.
   - `E_total ≡ E_SCF + E_c` to 1e-12.
   - RI vs Poisson foothold `|ΔE_c| ≤ 1e-10` Ha (same kernel).
   - Frozen pin band on `E_c` (offline CP2K when available; else
     first-green setup-integrity band ≤ 1e-3 Ha) — see REFERENCE.md.
   - Wall ≤ 60 s single process.

## Perf

`T_atv / T_oracle ≤ 10` when a fair CP2K wall exists; else `perf_exempt`
with reason recorded in REFERENCE.md.

## Done when

- Spec closed (delete this directory on merge).
- `acceptance_met=true` in worker schema.
- Commit on `wf/periodic-rpa-si`.
