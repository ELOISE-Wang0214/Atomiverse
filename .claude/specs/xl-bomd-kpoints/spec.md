# Spec — xl-bomd-kpoints

**Slug:** `xl-bomd-kpoints`  
**Branch:** `wf/xl-bomd-kpoints`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Oracle:** offline CP2K pins + Niklasson literature (Table I)  
**Domain:** dynamics (AIMD propagator on multi-k RKS GPW forces)

## Problem

`xl-bomd-water` / `xl-bomd-uks` productized Niklasson dissipative XL-BOMD
only on the **Γ real-density** path:

- `XlBomd` propagates `linalg::RealMatrix` only
- `GpwDft::run(..., guess_P)` is applied only after the multi-k branch
  returns early — `run_gpw_scf_kpoints` always cold-starts \(P_k\) from \(H_0\)
- `GpwResult` does not store complex \(P_k\) for the next MD step
- `GPWEvaluator` has no `k_mesh` member — multi-k MD cannot engage XL

Gap census (wave 10): *“Γ RKS/UKS XL ✓; multi-k XL product depth/pin open”*.

## Goal

Productize **Niklasson dissipative XL-BOMD on multi-k closed-shell RKS GPW**:

1. Propagate **per-k** auxiliary densities \(P_k\) (complex Hermitian AO
   matrices) with the same Table-I \((\kappa,\alpha,c_k)\) as Γ XL.
2. Seed the multi-k SCF with those \(P_k\) each MD step.
3. Keep SCF-converged XL (full `tol` each step); OT still refused.
4. Expose Monkhorst–Pack `k_mesh` on `GPWEvaluator` for the product MD path.
5. Regression: short **periodic H₂ NVE** on a small multi-k mesh with
   energy-drift envelope (same diagnostic class as `xl-bomd-water`).

## Non-goals

- Multi-k UKS XL (dual-spin \(P_k^\sigma\)) — residual after this slug
- OT + XL hybrid (still throws)
- Linearized / zero-SCF XL functional
- Metallic smearing + XL
- Live CP2K in ctest

## Design

### Propagator

- Keep `cpu::gpw::XlBomd` for real Γ densities.
- Add `cpu::gpw::XlBomdCplx` with the same Table-I coefficients and recursion
  on `linalg::CplxMatrix` (real \(\kappa,\alpha,c_k\) act element-wise; a
  Hermitian history stays Hermitian).
- Multi-k: one `XlBomdCplx` per sampled k (shared K).

### Solver API

```cpp
GpwResult& run(const std::vector<GpwAtom>& atoms,
               const linalg::RealMatrix* guess_C = nullptr,
               const linalg::RealMatrix* guess_P = nullptr,
               const linalg::RealMatrix* guess_P_beta = nullptr,
               const std::vector<linalg::CplxMatrix>* guess_P_k = nullptr);
```

| Path | Density guess |
| ---- | ------------- |
| Γ RKS | `guess_P` total real \(P\) |
| Γ UKS | `guess_P` / `guess_P_beta` |
| multi-k RKS | `guess_P_k` size == `n_kpts`, each `nb×nb` complex |
| multi-k + real `guess_P` | throw (no silent Γ guess on k path) |
| OT + any density guess | existing rules (RKS OT real only; multi-k OT+XL throw at evaluator) |

On multi-k convergence, store `GpwResult::P_k` (complex densities rebuilt
from the same occupations as forces).

### Evaluator

- Members: `int k_mesh[3] = {1,1,1}`, `bool k_time_reversal = true`
- `xl_bomd` + multi-k RKS allowed; multi-k UKS + XL throws at `build`
- State: `std::vector<XlBomdCplx> xl_k_` (empty until first multi-k SCF)
- `launch`: if `xl_k_` ready, pass `guess_P_k`; after converge, advance each
  channel from `last_.P_k`; reset all on SCF fail

### Regression

`regressions/dynamics/xl-bomd-kpoints` — periodic H₂ (DZV-GTH-PADE), mesh
`1×1×2` with `k_time_reversal=false` (two k), PADE, `xl_bomd=true`, K=5,
short NVE, relative energy-drift envelope + SCF every sampled step +
finite coords + soft bond window.

## Acceptance

See `acceptance.md`.

## Files (expected)

- `include/atomiverse/cpu/gpw/{xl_bomd,gpw_solver,gpw_evaluator}.h`
- `src/cpu/gpw/{xl_bomd,gpw_solver,gpw_evaluator}.cpp`
- `tests/cpu/gpw/test_gpw_xl_bomd.cpp`
- `regressions/dynamics/xl-bomd-kpoints/`
- `tools/cp2k_oracle/inputs/xl-bomd-kpoints/`
- `regressions/CMakeLists.txt`
