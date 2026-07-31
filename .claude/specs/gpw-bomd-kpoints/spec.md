# gpw-bomd-kpoints — Multi-k Born–Oppenheimer MD on GPW

Status: code-complete, acceptance met (wave 1, 2026-07-31) — all gates green:
unit 0.4 s (`test_cpu_gpw_bomd` ×3 incl. two new warm-start gates),
`reg_gpw-bomd-kpoints` PASS 2.14 s (E0 Δ=4.07e-4 Ha, F0 Δ=5.86e-3 Ha/Bohr =
Γ-model band, trace Δ=2.33e-4 Ha, drift 4.25e-3), perf ratio 0.064
(within_10x). Delete this spec on merge.
Branch: `wf/gpw-bomd-kpoints`
Chapter: dynamics
Oracle: CP2K 2026.1 Quickstep AIMD (`RUN_TYPE MD` + `&KPOINTS`), offline via
`tools/cp2k_oracle/` — never at ctest time.

## Problem

Γ-point BOMD (`gpw-bomd-water`) and multi-k **XL**-BOMD / CPMD / SGCP are
product, but **plain** multi-k BOMD (full SCF each step on a k-mesh, no
auxiliary propagator) is unpinned vs CP2K Quickstep AIMD with `&KPOINTS`,
and the plain multi-k MD path in `GPWEvaluator::launch` has **no warm start
at all**: after a multi-k SCF (`last_.n_kpts > 0`) neither the Γ `guess_C`
nor the Γ `guess_P` branch applies, and `guess_p_k` engages only through the
XL propagators. Every plain BOMD step is a cold multi-k SCF — correct, but a
gratuitous wall-time and convergence-robustness gap vs CP2K (whose MD always
seeds SCF from the previous step's density; `EXTRAPOLATION USE_PREV_P` is
exactly "previous converged density").

## Change

1. `src/cpu/gpw/gpw_evaluator.cpp`, plain (non-XL/CPMD/SGCP) path: when no
   XL k guess is active and the previous launch produced a **converged**
   multi-k result with a complete `P_k` payload, pass `&last_.P_k` (and
   `&last_.P_k_beta` for UKS — both or neither, no partial spin guess) to
   `GpwDft::run` as the initial k densities. Mirrors the existing Γ
   `guess_C`/`guess_P` warm-start block; solver-side `guess_P_k` validation
   (shape/count/partial-spin) already exists and is unit-tested.
2. No public-surface change: `GPWEvaluator` members are unchanged; the warm
   start is an internal property of the MD binding, exactly like the Γ one.

## Acceptance

- A1 (unit, `tests/cpu/gpw/test_gpw_bomd.cpp` → `test_cpu_gpw_bomd`):
  plain multi-k `GPWEvaluator` launched twice at the same geometry —
  second launch converges with **fewer SCF iterations** than the cold first
  and reproduces the cold energy to ≤ 1e-6 Ha (warm start reaches the same
  fixed point; regression against "cold every step").
- A2 (unit, same file): UKS multi-k second launch passes both spin guesses
  (converged, fewer iters than cold, same fixed point) — no partial-spin
  throw from the solver.
- A3 (e2e, `regressions/gpw/gpw-bomd-kpoints/`): step-0 potential energy of
  the periodic H₂ / DZV-GTH-PADE / PADE / 100 Ry / MP 1×1×2 twin matches the
  **CP2K 2026.1** MD step-0 potential (offline pin in REFERENCE.md) within
  the tolerance recorded there (same-method twin band, ~1e-3 Ha class).
- A4 (e2e): 4-step NVE trace — VV bootstrap seeded by a public single-point
  (forces at t=0, matching CP2K's step-0 force evaluation), per-step
  E_pot within the REFERENCE band of the CP2K `.ener` trace.
- A5 (e2e): NVE relative drift envelope ≤ 0.05 over 1 + 4 steps (family
  pattern of `gpw-bomd-water` / `xl-bomd-kpoints`), measured ≪ envelope.
- A6 (e2e): multi-k payload gates — `n_kpts == 2`, SCF converged at every
  sampled step, production-step `iters` ≤ cold step-0 `iters` (warm start
  never regresses the SCF; the strict same-geometry iteration win is A1's
  unit gate — at a *moved* geometry ATV's cold guess is already near-optimal
  for this 4-AO system, so strict-fewer does not hold per-step), finite
  coords, soft H–H window.
- A7 (process): regression is an out-of-library consumer, ≤ 60 s, no live
  oracle; registered `atv_add_regression(gpw gpw-bomd-kpoints THREADS 8)`
  behind `ATV_BUILD_GPW` + `ATV_BUILD_MD`.
- A8 (perf): fair wall vs CP2K on the same 4-step MD twin (same host, same
  OMP budget): `T_atv/T_cp2k <= 10`, or documented `perf_exempt` with
  reason.

## Twin definition (both engines)

| Knob | Value |
| ---- | ----- |
| System | H₂, bond 1.0 Bohr along x, cubic 12 Bohr cell, periodic XYZ |
| Basis / PP | DZV-GTH-PADE / GTH-PADE-q1 (CP2K `GTH_BASIS_SETS` / `GTH_POTENTIALS`) |
| XC | PADE (LDA) |
| CUTOFF | 100 Ry, NGRIDS 1 |
| k mesh | MONKHORST-PACK 1×1×2 (ATV `k_time_reversal=false` ⇔ CP2K `FULL_GRID T`) |
| SCF | tol 1e-7 (ATV) / `EPS_SCF 1e-7` (CP2K), DIIS diagonalization |
| MD | NVE, dt 0.5 fs, 4 steps, v₀ = ∓0.001 Å/fs antisymmetric stretch |
| CP2K guess | `EXTRAPOLATION USE_PREV_P` (twin of the previous-P_k warm start) |

## Notes

- `xl-bomd-linearized-uks` may touch the same file (linearized-force path);
  the plain-path warm-start block is disjoint from the linearized branch —
  coordinate at merge.
- Delete this spec on close (merge with acceptance met).
