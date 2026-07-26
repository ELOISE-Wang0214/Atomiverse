# Spec — xl-bomd-uks

**Slug:** `xl-bomd-uks`  
**Branch:** `wf/xl-bomd-uks`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Oracle:** offline CP2K pins + Niklasson literature (Table I)  
**Domain:** dynamics (AIMD propagator on Γ-UKS GPW forces)

## Problem

`xl-bomd-water` productized **closed-shell RKS Γ** Niklasson dissipative
XL-BOMD (`XlBomd` + `GPWEvaluator::xl_bomd`). Open-shell is still refused:

- `GPWEvaluator::build` throws when `xl_bomd && (n_alpha+n_beta > 0)`
- `GpwDft::run(..., guess_P)` throws for UKS
- SCF init applies `guess_P` only to the total RKS density `P`, never to
  spin densities \(P_\alpha / P_\beta\)

Gap census (wave 6): *“XL-BOMD closed-shell RKS Γ only; UKS/`use_ot` throw”*.
Foothold `xl-bomd-water` ≠ open-shell capability.

## Goal

Productize **Niklasson dissipative XL-BOMD on Γ-UKS GPW**:

1. Propagate **independent** auxiliary densities \(P^\alpha\) and \(P^\beta\)
   with the same Table-I \((\kappa,\alpha,c_k)\) as RKS.
2. Feed both as SCF initial spin densities each MD step.
3. Keep SCF-converged XL (full `tol` each step); OT still refused.
4. Regression: short **O₂ triplet NVE** energy-drift envelope (literature /
   BOMD conservation diagnostic — same envelope class as `xl-bomd-water`).

## Non-goals

- OT + XL hybrid (still throws)
- Multi-k XL-BOMD
- Linearized / zero-SCF XL functional
- CPMD / second-gen CP open-shell
- Live CP2K in ctest

## Design

### Propagator

Reuse `cpu::gpw::XlBomd` **twice** (α and β). Coefficients unchanged
(Niklasson JCP 130, 214109 (2009) Table I). No public `*Config`.

### Solver API

```cpp
GpwResult& run(const std::vector<GpwAtom>& atoms,
               const linalg::RealMatrix* guess_C = nullptr,
               const linalg::RealMatrix* guess_P = nullptr,
               const linalg::RealMatrix* guess_P_beta = nullptr);
```

| Path | `guess_P` | `guess_P_beta` |
| ---- | --------- | -------------- |
| RKS  | total \(P\) | must be null |
| UKS  | \(P_\alpha\) | \(P_\beta\) (both required if either set) |
| OT   | any non-null → throw | — |

UKS init: `Pa=*guess_P`, `Pb=*guess_P_beta`, `P=Pa+Pb` before first Fock.

### Evaluator

- `xl_bomd` + UKS occupations **allowed** (diagonalization path only)
- Members: `xl_` (RKS total or UKS α), `xl_beta_` (UKS β only)
- `launch`: pass both aux when ready; `advance` both from
  `last_.P` / `last_.P_beta` on converge; `reset` both on SCF fail

### Regression

`regressions/dynamics/xl-bomd-uks` — O₂ triplet (nα=7, nβ=5), MT box,
PADE or PBE, `xl_bomd=true`, K=5, short NVE, relative energy-drift envelope
+ SCF-converged-every-step + finite coords + soft O–O ceiling.

## Acceptance

See `acceptance.md`.

## Files (expected)

- `include/atomiverse/cpu/gpw/{gpw_solver,gpw_evaluator,xl_bomd}.h`
- `src/cpu/gpw/{gpw_solver,gpw_evaluator}.cpp`
- `tests/cpu/gpw/test_gpw_xl_bomd.cpp`
- `regressions/dynamics/xl-bomd-uks/`
- `tools/cp2k_oracle/inputs/xl-bomd-uks/`
- `regressions/CMakeLists.txt`
