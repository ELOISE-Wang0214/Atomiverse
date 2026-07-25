# Acceptance — `periodic-ump2-si`

## Must pass

1. **Unit** `test_cpu_gpw_mp2`:
   - UMP2 sum matches explicit spin-orbital formula on synthetic α≠β tensors
   - UMP2 ≡ RMP2 on identical α=β spatial integrals
   - cross-spin RI-GPW ≡ Poisson (`mo_eri_ovov_ab_ri` vs `mo_eri_ovov_ab`)
   - UKS without converged `run()` throws `logic_error` (not the old
     closed-shell-only `invalid_argument`)
2. **Regression** `reg_periodic-ump2-si` / `run_reg_periodic-ump2-si`:
   - Γ UKS-PBE SCF converges on sc-Si₂ (DZVP-MOLOPT-GTH, GTH-PBE-q4, 30 Ry)
   - `GpwDft::mp2` returns finite **negative** E_corr
   - OS+SS ≡ E_corr and E_total ≡ E_SCF+E_corr to 1e-12
   - virtual-space monotonicity with `mp2_n_virt = 8`
   - \|E_corr + 0.2764036524\| ≤ 2e-3 Ha (RMP2 dual-gate / offline band)
   - wall ≤ 60 s
3. **Architecture**: no public `*Config`, no `make_*` / `build_*` factories;
   Frame-adjacent I/O still atoms/GpwDft members only.
4. **Oracle**: offline pin notes + `tools/cp2k_oracle/inputs/periodic-ump2-si.inp`;
   never live CP2K in ctest.
5. **Perf**: `perf_exempt` if no CP2K wall twin; else
   \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\).

## Done when

All of the above green on `wf/periodic-ump2-si`; commit message cites slug;
spec ready for harness delete-on-merge.
