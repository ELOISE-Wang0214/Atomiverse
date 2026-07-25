# Acceptance — `gpw-kpoint-uks-stress`

## Gates

| ID | Check | Tol |
|----|-------|-----|
| A1 | Multi-k UKS SCF converges; `C_occ_k` / `C_occ_k_beta` present | — |
| A2 | Analytic vs multi-SCF FD max\|\σ_an−σ_fd\| on diag | ≤ 1×10⁻⁴ Ha/Bohr³ |
| A3 | Off-diagonal Voigt product ~ 0 | ≤ 1×10⁻⁸ |
| A4 | Closed-shell multi-k UKS (nα=nβ=1) vs multi-k RKS max\|\dσ\| | ≤ 5×10⁻⁴ Ha/Bohr³ |
| A5 | max\|\σ_diag\| magnitude floor | > 1×10⁻⁸ Ha/Bohr³ |
| A6 | Finite energy; soft \|E − first-green pin\| when frozen | ≤ 1×10⁻³ Ha |

## Unit

Replace `GpwStress.UksKpointThrows` with dual-gate unit on cheap high-spin H₂
+ k_mesh (or finite/no-throw + identity); keep Γ-UKS and multi-k RKS green.

## Perf

`perf_exempt` if no fair CP2K NUMERICAL wall on pin host (same policy as
`gpw-kpoint-stress` / `gpw-uks-stress`). Target r ≤ 10 when twin available.
Product wall should stay reasonable for dual-gate on H₂ (≪ O₂ water twins).

## Done when

- `gpw_stress` no longer throws multi-k UKS
- `ctest -R test_cpu_gpw_stress` green
- `run_reg_gpw-kpoint-uks-stress` PASS
- Commit on `wf/gpw-kpoint-uks-stress`
