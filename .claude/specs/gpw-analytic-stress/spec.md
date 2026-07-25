# Spec — `gpw-analytic-stress`

**Status:** in progress  
**Branch:** `wf/gpw-analytic-stress`  
**Queue:** `docs/developer/cp2k-method-queue.md` row `gpw-analytic-stress`  
**Dep:** `gpw-stress-water` ✓ (strain-FD oracle + CP2K NUMERICAL pin)

## Problem

GPW cell stress today is orthorhombic **multi-SCF strain-FD** only
(`gpw_stress` / `GPWEvaluator::compute_stress`; `reg_gpw-stress-water` vs
CP2K `STRESS_TENSOR NUMERICAL`). Product cost is **6+1 SCF** per stress
call. Analytic / Hellmann–Feynman cell stress (one SCF + frozen-density
strain response) is the product path for forces (`gpw_forces_analytic`)
and for plane-wave stress (`pw_stress`); GPW must match that product shape.

## Goal

Ship **analytic GPW stress** as the default path for orthorhombic diagonal
Voigt components, with multi-SCF strain-FD retained as oracle / twin.

## Science

At an SCF stationary point, cell derivatives of the total energy reduce to
**Hellmann–Feynman** derivatives of the energy functional at fixed density
matrix \(P\) (AO basis). For orthorhombic diagonal strain \(\varepsilon_{aa}\)
at **fixed fractional coordinates**:

\[
\sigma_{aa} = -\Omega^{-1}\,\partial E/\partial\varepsilon_{aa}\big|_{P},
\quad
E[P;\varepsilon] = \mathrm{Tr}(P\,T) + \mathrm{Tr}(P\,V_{\mathrm{pp}})
  + E_{\mathrm{es}}[\rho_P,n_c] + E_{\mathrm{xc}}[\rho_P] + E_{\mathrm{disp}}.
\]

Product implementation (GPW analogue of PW Nielsen–Martin in this tree):

1. One full SCF at the reference cell → converged \(P\).
2. Central FD of the **non-self-consistent** total energy \(E[P;\varepsilon]\)
   under diagonal strain \(\pm h\) (rebuild geometry-dependent \(T\), \(V_{\mathrm{pp}}\),
   grid, Poisson, XC; **do not** re-converge SCF).
3. Off-diagonal shear remains **0** (orthorhombic GPW cell only).

This is a **real** fixed-\(P\) derivative, not multi-SCF strain-FD under a
renamed API. Multi-SCF FD stays as `gpw_stress_fd` for dual-gate / oracle.

Sign convention unchanged (CP2K `STRESS_TENSOR`):
\(\sigma_{aa} = -\Omega^{-1}\partial E/\partial\varepsilon_{aa}\)
(E Hartree, \(\Omega\) Bohr³). Public Frame `kStress` is eV/Å³.

### Included (this slug)

- Closed-shell Γ RKS orthorhombic diagonal stress.
- `local_pp_on_grid` true/false consistent with SCF energy functional.
- Optional Grimme D3/D3BJ energy contribution under the same fixed-\(P\)
  geometry strain (geometry-only; no density response).
- Axes mask (`stress_axes`) — masked axes report exactly 0 (no FD samples).
- Default product path on `GPWEvaluator::compute_stress`.

### Deferred (document honestly)

- Fully closed-form Pulay AO-derivative stress without fixed-\(P\) energy FD
  (term-by-term Nielsen–Martin decomposition tables).
- Off-diagonal / triclinic shear.
- UKS, k-point GPW, DFT+U, D4, hybrid XC stress.
- Matching CP2K `STRESS_TENSOR ANALYTICAL` on nonperiodic MT (known to
  disagree grossly with NUMERICAL on the water twin — product pin remains
  NUMERICAL / multi-SCF FD, not CP2K ANALYTICAL).

## API (no public Config bags, no factories)

| Symbol | Role |
| ------ | ---- |
| `gpw_stress` | **Product** analytic fixed-\(P\) Voigt [Ha/Bohr³] |
| `gpw_stress_fd` | Multi-SCF strain-FD oracle (former `gpw_stress` body) |
| `GPWEvaluator::compute_stress` | Enable stress write |
| `GPWEvaluator::stress_use_fd` | `false` default → analytic; `true` → FD oracle |
| `GPWEvaluator::stress_strain_h` | Dimensionless step for fixed-\(P\) FD **and** multi-SCF FD |
| `GPWEvaluator::stress_axes` | Per-axis mask (unchanged) |
| `last_stress_ha_bohr3()` | Voigt [Ha/Bohr³] |
| Frame `kStress` | Voigt [eV/Å³] |

Internal: `internal/gpw_stress.h` exposes `gpw_stress_analytic_impl` and
`gpw_stress_fd_impl` (or equivalent names).

## Tests

| Tier | Path | Gate |
| ---- | ---- | ---- |
| Unit | `tests/cpu/gpw/test_cpu_gpw_stress.cpp` | Analytic vs FD diagonals on cheap H₂; axes mask; evaluator writes `kStress`; off-diag 0 |
| E2E | `regressions/gpw/gpw-analytic-stress/` | Water MT PBE twin: analytic vs multi-SCF FD within **1e-5 Ha/Bohr³**; energy pin vs existing CP2K E; optional soft CP2K NUMERICAL pin |

Existing `reg_gpw-stress-water` remains the offline CP2K NUMERICAL pin (may
run product analytic path once default flips — residual must stay in gate).

## Docs

- `regressions/gpw/gpw-analytic-stress/REFERENCE.md`
- `docs/tutorial/methods/gpw-analytic-stress.md` + index/README rows
- Roadmap / queue checkmarks on close

## Non-goals

- Public `*Config` / `make_*`
- CUDA GPW stress
- Changing Frame units or Voigt order
