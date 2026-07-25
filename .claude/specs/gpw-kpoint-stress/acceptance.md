# Acceptance — gpw-kpoint-stress

## Gates

| ID | Check | Tolerance |
| -- | ----- | --------- |
| U1 | Unit: k-mesh H₂ analytic vs multi-SCF FD max diag \(\|\sigma^{\mathrm{an}}-\sigma^{\mathrm{FD}}\|\) | ≤ **5×10⁻⁵** Ha/Bohr³ |
| U2 | Unit: off-diagonal product stress exactly 0 on orthorhombic cell | hard 0 / soft 1e-12 |
| U3 | Unit: multi-k no longer throws from `gpw_stress` | must return finite |
| A1 | Regression: SCF converges, finite energy and stress | hard |
| A2 | Regression: analytic vs multi-SCF FD max diag | ≤ **1×10⁻⁴** Ha/Bohr³ |
| A3 | Regression: off-diagonal ~ 0 | ≤ 1e-8 |
| A4 | Regression: absolute energy / stress pin (REFERENCE.md) | see REFERENCE |
| P1 | Perf: \(T_{\mathrm{atv}}/T_{\mathrm{oracle}}\le 10\) **or** `perf_exempt` | — |

## Dual gate

1. **Correctness:** offline pin + analytic↔FD identity (never live CP2K in ctest).
2. **Speed:** fair twin wall if CP2K present; else `perf_exempt` with reason
   (no CP2K binary on pin host — same as `gpw-kpoints-si`).

## Done when

- Spec acceptance green on `wf/gpw-kpoint-stress`.
- Scoped unit `test_cpu_gpw_stress` + `reg_gpw-kpoint-stress` pass.
- Commit on branch; worktree build wiped.
