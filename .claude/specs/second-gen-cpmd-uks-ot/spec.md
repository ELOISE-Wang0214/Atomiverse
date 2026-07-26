# Spec — `second-gen-cpmd-uks-ot`

## Goal

Productize **second-generation Car–Parrinello MD (SGCP)** on **Γ-UKS** with the
**Orbital Transformation (OT)** SCF as the partial-SCF corrector:

- Dual Kolafa ASPC predictors on \(P^\alpha, P^\beta\)
- Warm-start OT manifolds from **spin-resolved natural orbitals** of the ASPC
  densities (not silent ignore of the density guess; not prior-step \(C\) alone)
- Partial OT budget `second_gen_cp_max_iter` once ASPC is warm
- Analytic forces at the post-corrector density (`gpw_forces_from_scf`)
- NVE energy-drift envelope on O₂ triplet (same family as `second-gen-cpmd-uks`)

## Residual closed

| Before | After |
| ------ | ----- |
| RKS SGCP+OT ✓ (`second-gen-cpmd-ot`) | unchanged |
| UKS SGCP+DIIS ✓ (`second-gen-cpmd-uks`) | unchanged |
| UKS OT+ASPC **throws** at `GPWEvaluator::build` / `GpwDft::run` | **product** |

## Non-goals

- Multi-k SGCP; metallic smear + SGCP OT; CPMD / XL-BOMD hybrids
- Live CP2K in ctest (offline pin / sketch only)
- Molecular post-HF / classical MW / CUDA

## Acceptance

1. Unit: `build()` accepts `second_gen_cp + use_ot + n_alpha/n_beta`
2. Unit: `GpwDft::run` accepts UKS OT + dual density guess (NO seed)
3. Regression `reg_second-gen-cpmd-uks-ot`: NVE drift ≤ 0.10; production OT
   iters ≤ corrector budget; ⟨S²⟩ triplet window; O–O soft window
4. Offline CP2K twin sketch under `tools/cp2k_oracle/inputs/second-gen-cpmd-uks-ot/`
5. Dual gate: correctness green; `perf_exempt` if no CP2K binary for AIMD wall
   (same class as sibling SGCP dynamics slugs)
