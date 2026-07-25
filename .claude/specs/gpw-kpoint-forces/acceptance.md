# Acceptance — gpw-kpoint-forces

## Dual gate

### Unit (`tests/cpu/gpw/test_cpu_gpw_forces.cpp`)

| ID | Check | Tol |
| -- | ----- | --- |
| U1 | H₂ periodic PADE, `k_mesh={1,1,2}`: max\|F_analytic − F_fd\| | ≤ 5e-4 Ha/Bohr |
| U2 | Same system at explicit Γ through k machinery matches Γ assemble within | 1e-5 Ha/Bohr |
| U3 | `gpw_forces_analytic` with `k_mesh` does not throw / converges | — |

### Regression (`regressions/gpw/gpw-kpoint-forces/`)

| ID | Check | Tol |
| -- | ----- | --- |
| A1 | SCF converged, finite forces | — |
| A2 | Analytic vs central-FD energy gradient (ATV self) | max\|dF\| ≤ 1e-3 Ha/Bohr |
| A3 | Absolute forces / energy vs offline CP2K pin when present; else frozen first-green band | ≤ 2e-3 Ha/Bohr (F); 1e-3 Ha (E) |
| A4 | Translational residual \|∑F\| | ≤ 1e-3 Ha/Bohr |
| A5 | No external process at ctest; ≤ 60 s | — |

## Perf

- Fair wall vs CP2K ENERGY_FORCE when binary available: ratio ≤ 10 → `within_10x`.
- Else `perf_exempt=true` with reason (no module / no fair twin).

## Done when

- Spec acceptance green on `wf/gpw-kpoint-forces`.
- Commit on branch; `acceptance_met=true`.
