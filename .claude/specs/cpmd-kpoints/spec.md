# Spec — cpmd-kpoints

**Slug:** `cpmd-kpoints`  
**Branch:** `wf/cpmd-kpoints`  
**Loop:** CP2K capability (GPW / periodic DFT)  
**Oracle:** offline CP2K pins + Car–Parrinello literature  
**Domain:** dynamics (AIMD propagator on multi-k RKS GPW forces)

## Problem

Γ-point Car–Parrinello MD productizes as `cpu::gpw::Cpmd` +
`GPWEvaluator::{cpmd, cpmd_mu, cpmd_dt}` on closed-shell RKS (`reg_cpmd-water`,
RATTLE, electronic thermostat) and Γ-UKS (`reg_cpmd-uks`). The electronic
state is a **real** occupied block \(C\) with real \(H,S\).

Multi-k RKS SCF already stores Bloch \(C_k\) / \(P_k\), and multi-k XL-BOMD
propagates complex \(P_k\). **Multi-k CPMD** (fictitious dynamics of Bloch
\(C_k\)) is unpinned:

- `Cpmd` is real-only (`linalg::RealMatrix`)
- Fixed-density KS eval (`max_iter == 0`) throws on the k-path
- `GpwResult` does not export multi-k \(H_k,S_k\) for the electronic EOM
- `GPWEvaluator` seeds CPMD only from real `last_.C`

Gap census (wave 12): *“Γ RKS/UKS CPMD ✓; seed needs real C (multi-k unpinned)”*.

## Goal

Productize **multi-k closed-shell RKS Car–Parrinello MD (Bloch \(C_k\))**:

1. Propagate **per-k** occupied Bloch coefficients \(C_k\) with fictitious
   mass \(\mu\), complex electronic RATTLE, and S-tangent force
   \(\mu\,\ddot C_k = -(I - S_k C_k C_k^\dagger)\,(2\cdot\mathrm{occ}\,H_k C_k)\).
2. Cold multi-k SCF seed; production steps use **fixed-density** multi-k
   KS eval (`iters == 0`) + analytic k-summed forces at that density.
3. Conserved diagnostic \(E_{\mathrm{cons}} = T_n + T_e + E_{\mathrm{KS}}\)
   with \(T_e = \sum_k w_k\,\tfrac12\mu\lVert\dot C_k\rVert_F^2\).
4. Regression: short **periodic H₂ NVE** on mesh `1×1×2`
   (`k_time_reversal=false`), CP conserved-energy drift envelope.

## Non-goals

- Multi-k UKS CPMD (dual-spin \(C_k^\sigma\)) — residual
- Multi-k CPMD + electronic Nosé–Hoover thermostat product pin
- Second-gen multi-k CP (`second-gen-cpmd-kpoints` — separate slug)
- OT + CPMD hybrid
- Live CP2K in ctest

## Design

### Propagator

- Keep `cpu::gpw::Cpmd` for real Γ coefficients.
- Add `cpu::gpw::CpmdCplx` on `linalg::CplxMatrix` with the same VV +
  electronic RATTLE (Hermitian constraint \(C^\dagger S C = I\); velocity
  kills \(\mathrm{Herm}(C^\dagger S\dot C)\)).
- Multi-k: one `CpmdCplx` per sampled k (shared \(\mu\)); \(T_e\) is
  k-weight-summed.

### Solver

- Multi-k fixed-density: `max_iter == 0` + `guess_P_k` builds \(E[\{P_k\}]\),
  exports `GpwResult::{H_ks_k, S_k, P_k, n_kpts, k_frac, k_weights}`
  (`iters = 0`, `converged = true`). No diagonalization / mixing.
- Γ fixed-density path unchanged.

### Evaluator

- `cpmd` + multi-k RKS allowed (UKS multi-k + CPMD throws at `build`).
- State: `std::vector<CpmdCplx> cpmd_k_` (sized after cold multi-k SCF).
- Cold seed: full multi-k SCF → seed each channel from occupied \(C_k\).
- Warm step: build \(P_k = 2 C_k C_k^\dagger\) → fixed-density multi-k eval →
  align to \(S_k\) → forces via `gpw_forces_from_scf` with
  \(C_{\mathrm{occ},k}\) + \(\varepsilon_k = \mathrm{diag}(C_k^\dagger H_k C_k)\)
  → advance each \(C_k\).

### Regression

`regressions/dynamics/cpmd-kpoints` — periodic H₂ (DZV-GTH-PADE), mesh
`1×1×2` unreduced, PADE, `cpmd=true`, short NVE, relative \(E_{\mathrm{cons}}\)
drift envelope + production `iters==0` + `n_kpts>=2` + soft bond window.

## Acceptance

See `acceptance.md`.

## Files (expected)

- `include/atomiverse/cpu/gpw/{cpmd,gpw_solver,gpw_evaluator}.h`
- `src/cpu/gpw/{cpmd,gpw_solver,gpw_evaluator}.cpp`
- `tests/cpu/gpw/test_gpw_cpmd.cpp`
- `regressions/dynamics/cpmd-kpoints/`
- `tools/cp2k_oracle/inputs/cpmd-kpoints/`
- `regressions/CMakeLists.txt`
