# Acceptance — `periodic-wannier-gpw`

## Dual gate

| Gate | Criterion |
| ---- | --------- |
| Correctness | `reg_periodic-wannier-gpw` exits 0; unit `test_cpu_gpw_periodic_wannier` green |
| Speed | `perf_exempt` (no CP2K binary on pin host) **or** \(T_\mathrm{atv}/T_\mathrm{cp2k}\le 10\) |

## Hard checks (regression)

1. `scf_converged == true`
2. `n_wann == 8` and `centers.size() == 8`
3. max assigned min-image distance to BCT bond midpoints `< 0.90` Å
4. `total_spread < 0.5 * total_spread_initial` (and both finite, init > 0)
5. No external process; wall ≤ 60 s on Release

## Unit checks

1. `diamond_bct_bond_midpoints` / cubic SMV weight / `center_from_diagonal_phases`
2. `compute` throws on `k_mesh < 1`, `tol <= 0`, metallic smearing, UKS

## Done when

- Spec files present; implementation on `wf/periodic-wannier-gpw`
- Scoped build: `atomiverse`, `test_cpu_gpw_periodic_wannier`, `reg_periodic-wannier-gpw`
- Both gates green → `acceptance_met=true`
