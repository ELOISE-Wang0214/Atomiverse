# Spec — xl-bomd-uks-kpoints

**Slug:** `xl-bomd-uks-kpoints`  
**Branch:** `wf/xl-bomd-uks-kpoints`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Oracle:** offline CP2K pins + Niklasson literature (Table I)  
**Domain:** dynamics (AIMD propagator on multi-k UKS GPW forces)

## Problem

Product XL-BOMD covers:

- Γ RKS (`xl-bomd-water`) and Γ UKS dual-spin (`xl-bomd-uks`) via real `XlBomd`
- multi-k **RKS** (`xl-bomd-kpoints`) via per-k `XlBomdCplx` + `guess_P_k`

**Residual:** multi-k **UKS** + XL still throws at `GPWEvaluator::build`
(“dual-spin \(P_k^\sigma\) is a residual”) and `GpwDft::run` refuses
`guess_P_k` on the k-UKS path. Gap census wave 11: *multi-k RKS XL +
Γ-UKS XL ✓; multi-k UKS XL throw*.

## Goal

Productize **Niklasson dissipative XL-BOMD on multi-k collinear UKS GPW**:

1. Store dual-spin complex densities \(P_k^\alpha, P_k^\beta\) on
   `GpwResult` after multi-k UKS SCF (rebuild from final MOs, occ=1).
2. Accept dual-spin multi-k density guesses
   (`guess_P_k` = α, `guess_P_k_beta` = β; both or neither).
3. Propagate **two** `XlBomdCplx` histories per sampled k (α and β) with
   Table-I \((\kappa,\alpha,c_k)\); feed them as the next SCF guess.
4. SCF-converged XL only (full `tol`); OT / CPMD / SGCP still mutually
   exclusive with XL.
5. Regression: short multi-k H₂ NVE energy-drift envelope on the UKS path
   with closed-shell fillings (\(n_\alpha=n_\beta=1\); dual non-empty
   \(P_k^\sigma\); high-spin H₂ is unbound).

## Non-goals

- OT + multi-k UKS XL hybrid
- Metallic smearing + XL
- Linearized / zero-SCF XL functional
- Live CP2K in ctest
- Changing Γ XL or multi-k RKS XL behaviour beyond shared docs

## Design

### Propagator

Reuse `cpu::gpw::XlBomdCplx` (no new algebra type). Multi-k UKS:

- `xl_k_[ik]` advances \(D_k^\alpha\)
- `xl_k_beta_[ik]` advances \(D_k^\beta\)

Same \(K\) and Table I as RKS multi-k XL.

### Solver API

```cpp
GpwResult& run(..., 
  const std::vector<linalg::CplxMatrix>* guess_P_k = nullptr,
  const std::vector<linalg::CplxMatrix>* guess_P_k_beta = nullptr);
```

| Path | Density guess |
| ---- | ------------- |
| multi-k RKS | `guess_P_k` only (`guess_P_k_beta` must be null) |
| multi-k UKS | both `guess_P_k` (α) and `guess_P_k_beta` (β), size `n_kpts` each |
| multi-k UKS partial | throw (no α-only silent cold-start of β) |
| multi-k + real `guess_P` | throw (unchanged) |

`GpwResult::P_k` — multi-k RKS total **or** multi-k UKS \(P_k^\alpha\).  
`GpwResult::P_k_beta` — multi-k UKS \(P_k^\beta\) (empty for RKS / Γ).

### Evaluator

- Remove multi-k UKS + XL `build()` throw.
- State: `xl_k_` (α / RKS) + `xl_k_beta_` (UKS β only).
- `launch`: when both ready, pass dual `guess_P_k` / `guess_P_k_beta`;
  after converge, advance both from `last_.P_k` / `last_.P_k_beta`.

### Regression

`regressions/dynamics/xl-bomd-uks-kpoints` — periodic H₂ (DZV-GTH-PADE,
\(n_\alpha=1\), \(n_\beta=1\) on UKS path), mesh `1×1×2` unreduced,
PADE, `xl_bomd=true`, K=5, short NVE, relative energy-drift envelope +
dual-spin k payload diagnostics + SCF every sampled step.

## Acceptance

See `acceptance.md`.

## Files (expected)

- `include/atomiverse/cpu/gpw/{gpw_solver,gpw_evaluator,xl_bomd}.h`
- `src/cpu/gpw/{gpw_solver,gpw_evaluator}.cpp`
- `tests/cpu/gpw/test_gpw_xl_bomd.cpp`
- `regressions/dynamics/xl-bomd-uks-kpoints/`
- `tools/cp2k_oracle/inputs/xl-bomd-uks-kpoints/`
- `regressions/CMakeLists.txt`
