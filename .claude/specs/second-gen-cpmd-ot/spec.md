---
slug: second-gen-cpmd-ot
status: in_progress
created: 2026-07-26
domain: dynamics / gpw
queue: cp2k-method-queue
depends_on: [second-gen-cpmd-water, ot-minimizer-water]
oracle: CP2K
---

# Spec: `second-gen-cpmd-ot` — Second-gen CP with OT partial-SCF corrector

**Slug:** `second-gen-cpmd-ot`  
**Domain:** AIMD / GPW (CPU)  
**Queue:** `docs/developer/cp2k-method-queue.md`  
**Oracle:** offline CP2K (`tools/cp2k_oracle/`)

## Summary

Productize **second-generation Car–Parrinello MD (SGCP)** with the
**Orbital Transformation (OT)** SCF minimizer as the partial-SCF
corrector. Today `GPWEvaluator::{second_gen_cp, use_ot}` is rejected at
`build()` and `GpwDft::run` refuses ASPC density guesses when
`use_ot=true`. Enable closed-shell Γ-point **ASPC predictor + OT
corrector** on the same H₂O NVE twin as `second-gen-cpmd-water`.

## Domain basis

Kühne et al., *Phys. Rev. Lett.* **98**, 066401 (2007): ASPC density
extrapolation + limited SCF corrector (no fictitious electronic mass).

Kolafa ASPC (`cpu::scf::ASPCPredictor`) predicts \(P\). The corrector is
either diagonalization+DIIS (product) or **OT** (this slug):

1. **Warm-start:** natural orbitals of the ASPC \(P\) seed the OT occupied
   manifold (density → \(C_0\)); previous \(C\) is insufficient alone —
   ASPC *is* the predictor.
2. **Corrector budget:** `max_iter = second_gen_cp_max_iter` (≪ full BOMD).
3. **Forces:** analytic forces at the post-corrector density (partial SCF
   may not hit `tol`; same semantics as DIIS SGCP).
4. **Record:** Kolafa correct + ASPC history of the post-OT density.

OT formulation in-tree: VandeVondele–Hutter / Kühne OT manifold
(`cpu/gpw/internal/ot.h`). Same gapped fixed point as DIIS.

## Public API (members only)

```cpp
cpu::gpw::GPWEvaluator eval({O, H}, /*cutoff*/200, /*max_iter*/60, /*tol*/1e-5,
                            /*mix*/0.3, XcKind::PBE, /*nonperiodic*/true,
                            /*compute_forces*/true);
eval.second_gen_cp = true;
eval.second_gen_cp_k = 1;
eval.second_gen_cp_max_iter = 5;
eval.use_ot = true;   // NEW: accepted with second_gen_cp
```

No public `*Config`, no `make_*` / `build_*` factories. Frame in / energy
and forces out (eV, eV/Å).

## Files (expected)

| Path | Change |
| ---- | ------ |
| `include/.../gpw_evaluator.h` | docs: `second_gen_cp` allows OT corrector |
| `src/cpu/gpw/gpw_evaluator.cpp` | drop OT reject; force path OT-safe |
| `src/cpu/gpw/gpw_solver.cpp` | accept `guess_P` with OT via natural-orbital warm-start; export `P`/`H_ks`/`S` on OT exit |
| `include/.../gpw_solver.h` | docs for density guess + OT |
| `tests/cpu/gpw/test_gpw_second_gen_cpmd.cpp` | accept OT; keep XL/CPMD rejects |
| `regressions/dynamics/second-gen-cpmd-ot/` | e2e NVE pin |
| `tools/cp2k_oracle/inputs/second-gen-cpmd-ot/` | offline twin sketch |

## Out of scope

- UKS OT + SGCP (dual ASPC + dual OT manifolds — follow-on if needed)
- Multi-k / smear / DFT+U with SGCP+OT
- XL-BOMD + OT, classical CPMD + OT
- Live CP2K in ctest

## Dual gate

1. **Correctness:** NVE drift envelope + production OT iters ≤ corrector
   budget + O–H integrity (literature / method-identity class of
   `second-gen-cpmd-water`).
2. **Perf:** `perf_exempt` if no CP2K wall on pin host (same as sibling
   SGCP dynamics slugs).
