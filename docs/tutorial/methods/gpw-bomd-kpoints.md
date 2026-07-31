# Plain multi-k Born–Oppenheimer MD on GPW (H₂ NVE vs CP2K)

*(Technical-report page; CP2K capability loop wave 1, 2026-07-31.)*

| Field | Value |
| ----- | ----- |
| Slug | `gpw-bomd-kpoints` |
| Regression | `regressions/gpw/gpw-bomd-kpoints/` |
| Public API | `cpu::md::MD` + `VelocityVerlet` + `cpu::gpw::GPWEvaluator{kinds, cutoff_ry, …, compute_forces=true}` with `k_mesh = {1,1,2}`, `k_time_reversal = false`, bound via `set_evaluator_ref`; step-0 forces seeded by `cpu::single_point(frame, gpw)` |
| Oracle | CP2K 2026.1 Quickstep AIMD (`RUN_TYPE MD` + `&KPOINTS`), offline via `tools/cp2k_oracle/inputs/gpw-bomd-kpoints/` |
| Roadmap / queue | `docs/developer/cp2k-method-queue.md` wave-13 row |

## Summary

Plain Born–Oppenheimer MD on a Monkhorst–Pack k-mesh: **full SCF to `tol`
each step** on the multi-k RKS path with analytic k-summed forces — no
auxiliary propagator (contrast `xl-bomd-kpoints`, `cpmd-kpoints`,
`second-gen-cpmd-kpoints`). The implementation delta of this slug is the
**previous-step \(P_k\) warm start** on the plain multi-k MD path of
`GPWEvaluator::launch`: after a converged multi-k step, `{P_k}` (UKS: both
\(P_k^\sigma\) — the solver rejects a partial spin guess) seeds the next
SCF, mirroring the Γ `guess_C`/`guess_P` block and CP2K's
`EXTRAPOLATION USE_PREV_P`. Previously every plain multi-k BOMD step was a
cold SCF.

## Twin and pins

Periodic H₂ (bond 1.0 Bohr, cubic 12 Bohr cell), DZV-GTH-PADE /
GTH-PADE-q1 / PADE, CUTOFF 300 Ry NGRIDS 1, MP 1×1×2 unreduced (2 k),
NVE dt 0.5 fs × 4 steps. The CP2K data-file basis `BASIS_SET → H
DZV-GTH-PADE` is numerically identical to ATV's in-tree table, so the twin
is same-basis / same-PP. VV bootstrap parity: a public single point seeds
forces at \(t = 0\), so ATV step \(n\) ↔ CP2K `.ener` step \(n\).

Pinned vs **CP2K 2026.1** (offline): step-0 energy (Δ = 4.1e-4 Ha measured),
step-0 force (Δ = 5.9e-3 Ha/Bohr — the Γ-point periodic-GPW model band; the
Γ control shows identical deltas, i.e. the k machinery adds < 1e-5), the
per-step E_pot trace (worst Δ = 2.3e-4 Ha), and the NVE drift envelope
(ATV 4.3e-3 vs CP2K 4.4e-3 relative at this deliberately aggressive dt).

## Gate

Unit (`test_cpu_gpw_bomd`): same-geometry relaunch converges in strictly
fewer iterations to the same fixed point (RKS and UKS multi-k); MD binding
gate on cheap periodic H₂. E2E (`reg_gpw-bomd-kpoints`): the four pins
above + `n_kpts == 2`, converged every step, warm-start-never-worse, soft
H–H window, ≤ 60 s.

Ratio: **0.064** (ATV 2.14 s vs CP2K 33.4 s, same host, OMP=8, single
rank) — within_10x, no exemption.

See `regressions/gpw/gpw-bomd-kpoints/REFERENCE.md`.
