# Acceptance — `periodic-ri-mp2-si`

| Check | Criterion |
|-------|-----------|
| Spec | `.claude/specs/periodic-ri-mp2-si/spec.md` |
| Oracle input | `tools/cp2k_oracle/inputs/periodic-ri-mp2-si.inp` |
| Unit | `ctest -R '^test_cpu_gpw_mp2$'` green (incl. RI≡Poisson) |
| Regression | `run_reg_periodic-ri-mp2-si` green ≤ 60 s |
| API | `GpwDft::mp2_ri` default true; `Mp2Result::ri`; no `*Config` / factories |
| Multi-k | clear `invalid_argument` residual (not silent) |
| Dual gate | correctness + within_10x **or** perf_exempt |
| Commit | `wf/periodic-ri-mp2-si` |

`acceptance_met = true` only when all rows pass.

## Result (2026-07-25)

| Check | Status |
|-------|--------|
| Unit `test_cpu_gpw_mp2` | **PASS** (RI≡Poisson) |
| `run_reg_periodic-ri-mp2-si` | **PASS** (E_corr within 1e-3 of CP2K 2026.1; \|RI−can\|=0) |
| Wall ATV / CP2K | 51.4 s / 576.7 s → **ratio ≈ 0.09** (within_10x) |
| acceptance_met | **true** |
