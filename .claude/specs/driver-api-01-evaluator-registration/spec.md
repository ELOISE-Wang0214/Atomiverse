---
slug: driver-api-01-evaluator-registration
status: approved
created: 2026-07-27
revised: 2026-07-27
domain: driver / public-api
queue: driver-api
deferred: false
depends_on: []
unblocks: [driver-api-02-aimd-wrappers]
grilled: true
---

# Spec: `driver-api-01-evaluator-registration` — one evaluator slot, no public `add`

**Slug:** `driver-api-01-evaluator-registration`
**Domain:** CPU (and matching CUDA surface) Driver registration API
**Queue:** driver-api
**Status:** approved (operator design 2026-07-27). **No backward compatibility.**

## Summary

Hard-cut the Driver registration surface so that:

1. **Primitive** Plan registration is only `on(phase).add` / `on(phase).add_ref`.
2. **Public `Driver::add` / `add_ref`** (implicit `default_bindings()` placement) are **removed** from the public API (protected `bind_defaults` remains for categorical helpers and criteria that declare stage bindings).
3. Energy/force sources register through one plural categorical entry: **`set_evaluator` / `set_evaluator_ref`** → `on(Traits::eval_phase).add(...)`.
4. **`add_potential` / `add_potential_ref` are deleted** (replaced by `set_evaluator`).
5. Categorical setters (`set_integrator`, `set_neighbor`, `set_observer`, `set_record`, `set_thermo`, …) **only** call `on().add` or protected `bind_defaults` — never a parallel public path.

Classical LJ, QMMM, GPW BOMD, nested HF methods, bias — all use `set_evaluator`. Identity is the **type**, not the setter name.

## Domain basis

Driver = staged loop + termination. Force/energy sources are components (or nested method-drivers) that write Frame channels. There is no scientific distinction that warrants separate `set_potential` vs `set_method` vs `set_force` registration APIs — only different **kinds** (`Potential`, `HF`, `GPW`, `Bias`, …) validated at `build()`.

Rejected: three parallel setters; public `add` that guesses phase; mode flags on methods (next spec).

## Design

### Public surface (presets: MD, Minimize, Evaluate, MC, …)

| API | Role | Plural? | Implementation |
| --- | ---- | ------- | -------------- |
| `on(phase).add(rvalue, order)` | primitive | yes | Plan pending |
| `on(phase).add_ref(lvalue, order)` | primitive borrow | yes | Plan pending |
| `set_evaluator(E&&)` | energy/force source @ `eval_phase` order 0 | **yes** | `on(eval_phase).add` |
| `set_evaluator_ref(E&)` | borrow | yes | `on(eval_phase).add_ref` |
| `set_neighbor` / `set_integrator` / … | categorical singular | no | `on(...).add` or `bind_defaults` |
| `set_observer` / `set_record` | categorical plural | yes | `bind_defaults` (requires `default_bindings`) |
| `set_thermo(N)` | convenience | singular | `on(end_phase).add(Thermo{N})` |
| `add_criterion` / `set_outcome` | policy slots | (as today) | not Plan energy |

### Deleted (no shim, no alias)

- Public `Driver::add` / `add_ref` / `add(unique_ptr)`
- `add_potential` / `add_potential_ref`
- Docs/examples/python that teach `sim.add(...)` for Plan components

### Protected

```cpp
// core/driver_template.h — was public add/add_ref
void bind_defaults(Comp&&);           // uses default_bindings()
void bind_defaults_ref(Comp&);
void bind_defaults(std::unique_ptr<Comp>);
```

`set_singular_default_owned_`, `add_criterion` (when criterion has phase bindings), `set_observer` call `bind_defaults*`.

### `set_evaluator` acceptance

Template forwards to `on(Traits::eval_phase).add`. Accepts:

- Any `Component` (pair, QMMM, GPWEvaluator, Bias, …)
- Nested pure `DriverNode` energy methods (HF tag methods that are not Components), via existing PhaseBuilder driver overloads

`build()` validation (existing kind walks) still requires ≥1 energy/force source for MD/Evaluate/Minimize/MC as today. No new kind enum.

### Preset ctors

Replace any `add(Migrate{...})` / `add(Zero*)` with explicit `on(phase).add(...)`.

## Files

| Path | Change |
| ---- | ------ |
| `include/atomiverse/core/driver_template.h` | public `add*` → protected `bind_defaults*` |
| `include/atomiverse/cpu/driver/driver.h` | `set_evaluator*`; delete `add_potential*`; observer/record/singular helpers |
| `include/atomiverse/cpu/dynamics/md.h` | ctor + docs; no public `add` usage |
| `include/atomiverse/cpu/minimize/minimize.h` | same |
| `include/atomiverse/cpu/evaluate/evaluate.h` | docs |
| `include/atomiverse/cpu/montecarlo/*` | `add_potential` → `set_evaluator` |
| CUDA mirror headers if present | same surface |
| tests / regressions / examples / python / docs / AGENTS.md / CLAUDE.md | mechanical migrate |

## Tasks

1. **API** — rename/protect `add` → `bind_defaults`; add `set_evaluator`; delete `add_potential`.
2. **Presets** — MD/Minimize/MC/Evaluate/Path ctors use only `on().add` / categorical.
3. **Migrate tree** — every `add_potential` → `set_evaluator`; every Plan `driver.add(comp)` → `set_*` or `on(phase).add`.
4. **Scrubs** — architecture greps fail on `add_potential` and public `\.add\(` misuse if any scrub exists; update teaching docs.
5. **Build** — CPU configure + compile tests that exercise MD/Evaluate registration.

## Testing

- Unit: existing MD/neighbor/thermo tests compile and pass with new API.
- Grep gate (optional test script): no `add_potential` in `include/` / `src/` / `tests/` / `regressions/` / `examples/`.
- No e2e under `tests/`.

## Out of scope

- Extracting CPMD/XL/SGCP flags into wrapper types (`driver-api-02-aimd-wrappers`).
- Changing kind hierarchy or zeroing contract.
- ForceField data-table redesign.
