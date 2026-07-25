# Acceptance — `cpmd-rattle-water`

## Unit (`tests/cpu/gpw/test_gpw_cpmd.cpp` + RATTLE cases)

| ID | Gate |
| -- | ---- |
| A1 | Existing foothold guards still green (`Cpmd` mass, seed, advance-before-seed, GPWEvaluator rejects) |
| A2 | After `advance` with nontrivial \(H\): \(\lVert C^\top S C - I\rVert_\infty \le 1\mathrm{e}{-9}\) |
| A3 | After `advance`: \(\lVert\mathrm{sym}(C^\top S \dot C)\rVert_\infty \le 1\mathrm{e}{-9}\) (velocity constraint) |
| A4 | `last_position_iters() ≥ 0` and finite; at-rest \(H=0\) keeps \(T_e\approx 0\) |
| A5 | `last_constraint_residual()` matches \(\lVert C^\top S C - I\rVert_\infty\) after advance |

## Regression (`regressions/dynamics/cpmd-rattle-water`)

| ID | Gate |
| -- | ---- |
| A10 | Public API: `MD` + `VelocityVerlet` + `GPWEvaluator{…, compute_forces, cpmd=true}` (same engagement as foothold) |
| A11 | After warm-up: finite PE \(<-100\) eV, finite nuclear KE ≥ 0 |
| A12 | Production: finite coords; O–H \(< 1.8\) Å |
| A13 | \(E_{\mathrm{cons}}=KE_{\mathrm{n}}+PE+T_{\mathrm{e}}\): relative drift \(\le 0.05\) (tighter than Löwdin foothold 0.10) |
| A14 | Production force steps: `last_result().iters==0` (not BOMD) |
| A15 | Wall ideally ≤60 s; else `perf_exempt` documented in REFERENCE |
| A16 | REFERENCE cites Andersen RATTLE + Car–Parrinello; offline CP2K twin input under `tools/cp2k_oracle/inputs/` |

## Docs / queue

| ID | Gate |
| -- | ---- |
| A17 | Tutorial `docs/tutorial/methods/cpmd-rattle-water.md` (honest RATTLE scope) |
| A18 | On merge: roadmap CPMD RATTLE ✓; queue row done; spec deleted |

## Green criteria (worker schema)

- Dual gate: unit tests + `run_reg_cpmd-rattle-water`
- `within_10x` or `perf_exempt` (literature envelope pin when no CP2K wall twin)
- `acceptance_met=true` only if A1–A16 hold
