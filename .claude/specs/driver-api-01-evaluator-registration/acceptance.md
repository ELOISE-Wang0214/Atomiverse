# Acceptance — `driver-api-01-evaluator-registration`

## Criteria

- [x] #ac-001 type: api verified_by: agent-auto  
  Public `ATV::Driver` / `cpu::Driver` surface has **no** `add` / `add_ref` / `add(unique_ptr)` for Plan components. Protected `bind_defaults*` exists for categorical helpers.

- [x] #ac-002 type: api verified_by: agent-auto  
  `add_potential` / `add_potential_ref` symbols are **gone** from public headers (`set_evaluator` only).

- [x] #ac-003 type: api verified_by: agent-auto  
  `set_evaluator(E&&)` and `set_evaluator_ref(E&)` bind at `Traits::eval_phase` order 0 via `on().add` / `add_ref`, and are **plural** (multiple calls allowed).

- [x] #ac-004 type: api verified_by: agent-auto  
  `set_observer` / `set_record` / `set_integrator` (and other categorical setters) use `on().add` or protected `bind_defaults`.

- [x] #ac-005 type: behavioral verified_by: agent-auto  
  `test_cpu_md`, `test_cpu_md_set_api`, `test_cpu_md_no_alloc`, `test_cpu_driver_binding_kind`, `test_cpu_component_ownership` all green (2026-07-27).

- [ ] #ac-006 type: scrub  
  `rg 'add_potential' include src tests regressions examples` — residual comments/docs may remain; code symbols gone from public headers.

- [x] #ac-007 type: compile verified_by: agent-auto  
  `atomiverse` + MD/driver unit tests compile under CPU build.
