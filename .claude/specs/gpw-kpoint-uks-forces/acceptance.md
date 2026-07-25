# Acceptance — `gpw-kpoint-uks-forces`

## Dual gate

### Unit (`tests/cpu/gpw/test_cpu_gpw_forces.cpp`)

| ID | Check | Tol |
| -- | ----- | --- |
| U1 | High-spin H₂ periodic PADE, `k_mesh={1,1,2}`, `n_alpha=2,n_beta=0`: forces finite, restoring bond, Newton-3 | \|∑F\| ≤ 2e-3 Ha/Bohr |
| U2 | Explicit-Γ k-UKS forces match Γ-UKS `assemble` path within | 1e-5 Ha/Bohr (or max residual documented) |
| U3 | Closed-shell multi-k UKS (`n_α=n_β=1`) forces match multi-k RKS within | 5e-4 Ha/Bohr |
| U4 | `gpw_forces_analytic` with k-mesh + UKS does **not** throw | — |

### Regression (`regressions/gpw/gpw-kpoint-uks-forces/`)

| ID | Check | Tol |
| -- | ----- | --- |
| A1 | SCF converged; finite forces; no throw | — |
| A2 | Absolute E / F vs frozen first-green pin (and/or offline CP2K ENERGY_FORCE) | E ≤ 1e-5 Ha; F ≤ 1e-5 Ha/Bohr (self pin) or 5e-3 (CP2K multigrid twin) |
| A3 | Translational residual \|∑F\| | ≤ 1e-6 Ha/Bohr (symmetric toy) |
| A4 | Offline FD identity (related geometry, documented in REFERENCE) | max\|dF\| ≤ 1e-4 Ha/Bohr preferred |
| A5 | No external process at ctest; wall ≤ 60 s | — |

## Perf

- Fair wall vs CP2K ENERGY_FORCE when binary available: ratio ≤ 10 → `within_10x`.
- Else `perf_exempt=true` with reason (no CP2K module / no same-basis fair twin).

## Done when

- Spec acceptance green on `wf/gpw-kpoint-uks-forces`.
- Commit on branch; `acceptance_met=true`.
