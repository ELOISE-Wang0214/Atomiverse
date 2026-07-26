# Acceptance — second-gen-cpmd-kpoints

## Correctness

- [x] `ASPCPredictorCplx` unit: record/predict/correct match real ASPC on
      real-valued matrices (complex imag = 0) and extrapolate a complex
      linear trend. (`test_cpu_scf_accelerators` ASPCCplx_*)
- [x] `test_cpu_gpw_second_gen_cpmd`: `build()` accepts multi-k
      `second_gen_cp`; rejects multi-k + `use_ot` + `second_gen_cp`.
- [x] `reg_second-gen-cpmd-kpoints` (dynamics):
  - multi-k RKS SGCP NVE on periodic H₂ (`k_mesh={1,1,2}`, unreduced)
  - after ASPC warm-up, production SCF iters ≤ `second_gen_cp_max_iter`
  - `n_kpts >= 2`, `P_k.size() == n_kpts`
  - relative NVE drift envelope ≤ 0.10 (soft; partial-SCF heat)
    (measured ~5.9e-4 on pin host)
  - soft H–H window; finite PE/KE

## Perf

- [x] Offline CP2K twin sketch under `tools/cp2k_oracle/inputs/second-gen-cpmd-kpoints/`
- [x] `perf_exempt = true` if no `cp2k` binary on pin host (same inventory as
      `xl-bomd-kpoints` / `second-gen-cpmd-water`); else
      \(T_{\mathrm{atv}}/T_{\mathrm{oracle}} \le 10\)
      → **perf_exempt** (no CP2K on host 2026-07-26)

## Architecture

- [x] No public `*Config` / `make_*` / `build_*`
- [x] Frame I/O only; energy on Frame in eV
- [x] Spec under `.claude/specs/second-gen-cpmd-kpoints/`; deleted on harness close
