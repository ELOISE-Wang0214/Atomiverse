# Acceptance — `driver-api-02-aimd-wrappers`

## Criteria

- [x] #ac-001 type: api verified_by: agent-auto  
  Public `GPWEvaluator` exposes **no** `cpmd` / `xl_bomd` / `second_gen_cp` assignable members. Accessors `is_cpmd()` / `is_xl_bomd()` / `is_second_gen_cp()` only.

- [x] #ac-002 type: api verified_by: agent-auto  
  `cpu::electronic::CarParrinello` constructs with method + μ (optional thermostat); implements `TimestepConsumer`; registers as `GPW` kind.

- [x] #ac-003 type: api verified_by: agent-auto  
  `XlBomdForce` and `SecondGenCp` ctor-only; mutual exclusion by type.

- [x] #ac-004 type: behavioral verified_by: agent-auto  
  `test_cpu_gpw_cpmd` and `test_cpu_gpw_xl_bomd` green (worktree 2026-07-27).

- [x] #ac-005 type: api verified_by: agent-auto  
  Mean-field: `CarParrinelloMethod<M>` / `XlBomdMethod<M>`; AIMD flags protected on `MeanFieldMethod`.

- [ ] #ac-006 type: scrub  
  Full-tree `\.cpmd\s*=` residual on non-product docs/regressions may remain until sweep.

- [x] #ac-007 type: compile verified_by: agent-auto  
  Worktree CPU build of CP/XL unit tests succeeds.
