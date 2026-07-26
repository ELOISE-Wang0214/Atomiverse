---
slug: second-gen-cpmd-uks
status: in_progress
created: 2026-07-26
domain: dynamics / gpw
queue: cp2k-method-loop
oracle: CP2K
depends_on: [second-gen-cpmd-water, xl-bomd-uks]
---

# Spec: `second-gen-cpmd-uks` — Second-gen CP (ASPC + partial SCF) on Γ-UKS

**Slug:** `second-gen-cpmd-uks`  
**Domain:** GPW AIMD (second-generation Car–Parrinello)  
**Oracle:** offline CP2K (ASPC + limited MAX_SCF UKS MD)  
**Status:** in progress (wave 7 worker)

## Summary

Extend product **second-generation Car–Parrinello** (`GPWEvaluator::second_gen_cp`)
from closed-shell RKS Γ to **collinear Γ-UKS**: dual Kolafa ASPC predictors on
\(P^\alpha, P^\beta\) + partial SCF corrector budget + analytic forces at the
post-corrector density. RKS path stays bit-compatible.

## Domain basis

Kühne et al., *Phys. Rev. Lett.* **98**, 066401 (2007); Kolafa, *J. Comput.
Chem.* **25**, 335 (2004). CP2K Quickstep: `EXTRAPOLATION ASPC` + limited
`MAX_SCF` on MD steps. Open-shell: independent spin-density histories, same
corrector SCF budget as RKS.

## Scope

**In**

- `second_gen_cp` with `n_alpha`/`n_beta` > 0 on diagonalization UKS path
- Dual `ASPCPredictor` (α total channel uses `last_.P`; β uses `last_.P_beta`)
- Partial-SCF export on UKS budget exhaust (P/C/eps for force assembly)
- Forces via `gpw_forces_from_scf` at corrector density (BO not required)
- Unit guards + `reg_second-gen-cpmd-uks` O₂ triplet NVE envelope

**Out**

- OT corrector, k-points, classical CPMD UKS (`cpmd-uks`), XL-BOMD (already product)
- Live CP2K in ctest

## Acceptance

| ID | Gate |
| -- | ---- |
| A1 | `build()` accepts `second_gen_cp` + UKS; still rejects OT / CPMD / XL |
| A2 | Dual ASPC warm steps use `second_gen_cp_max_iter` corrector budget |
| A3 | Unit: RKS SGCP still accepted; UKS SGCP accepted |
| A4 | `reg_second-gen-cpmd-uks`: NVE drift ≤ 0.10; production iters ≤ budget; ⟨S²⟩ triplet window |
| A5 | perf_exempt if no CP2K binary; else T_atv/T_oracle ≤ 10 |

## Public API

```cpp
cpu::gpw::GPWEvaluator eval(...);
eval.n_alpha = 7;
eval.n_beta = 5;
eval.second_gen_cp = true;
eval.second_gen_cp_k = 1;
eval.second_gen_cp_max_iter = 5;
eval.compute_forces = true;
```

No new public types / factories / Config bags.
