# Frame is the only public token (no user Context)

**Status:** binding (2026-07-29). Router: `CLAUDE.md` iron law
`frame-only-token`.

## Rule

Public callers exchange **one** system object: **`Frame`**. They do **not**
author or pass a parallel Context bag. Parameters live on the method/component;
orchestration is composition of Drivers/Components; membership is Frame fields.

## Why

A user `*Context` reintroduces the deleted parallel world-model (`*System`,
ambient “session” structs): every new feature grows another
field on the bag, Frame stops being the contract, and examples stop looking
like the architecture. Engine-side Context is a **Runtime publication** for
`launch()` — that is not a license for a second user-facing token.

## Allowed engine types (not user API)

| Type | Role |
|------|------|
| `cpu::Context` / `cuda::Context` / `torch::Context` | Flat core + service slices; Runtime fills; components read |
| `LaunchContext` / `BuildContext` / `DestroyContext` | Component lifecycle view of the above |
| `SolverContext` (CUDA SCF isolation) | Narrow engine-only solver handle |

Users never `#include` these to *drive* a calculation. Tutorials and
`regressions/` construct Frame + methods + presets only.

## Forbidden shapes

- `struct FooContext { … }` / `class BarContext` on the public surface or in
  user-facing examples as the way to carry system + params
- `run(Context&)` / `run(Frame&, Context&)` / `set_context`
- Parking region masks, charges, link lists, or results on a context instead
  of Frame columns / meta
- “Just for this method” private bags that leak into public headers under a
  Context name

## Placement checklist

| Question | Answer |
|----------|--------|
| Is it system state or a result of the step? | **Frame** (`blocks` / `meta` / `box`) |
| Is it an algorithm knob? | ctor arg or **member** on the owning type |
| Is it another algorithm piece? | categorical setter / composition |
| Is it neighbor list / DD / pair list / loop status? | engine Context slice (not user) |

## molrs exchange (symmetric names)

molrs.Frame is the **exchange** layer (string block/column API). `cpu::Frame`
is the **execution** layer (`span` / `meta.get`). Boundary free functions
(once per `run`, never in the step loop):

| Direction | API |
|-----------|-----|
| molrs → `cpu::Frame` | `cpu::from_molrs(FrameRef)` |
| `cpu::Frame` → molrs | `cpu::to_molrs(Frame, FrameRef)` |

Same pair on CUDA/torch: `cuda::from_molrs` / `cuda::to_molrs`. Header:
`cpu/frame/molrs_exchange.h`. Do **not** reintroduce asymmetric
`materialize` / `write_back` as the public exchange names.

## Gate

`tests/architecture/no_user_context_scrub.sh` (`ctest -R test_no_user_context_scrub`):

- Bans **construction** of `cpu::Context` / `LaunchContext` / `cpu::Runtime`
  under `regressions/`, `examples/`, `app/`, `docs/tutorial/`
- Bans direct `#include` of `cpu/component/context.h` there
- Bans hand-driving: `.launch(`, brace `build({…})`, `setup_pass` /
  `step_once` / `teardown_pass`
- Allows Component **override signatures** (`void launch(const LaunchContext&)`)
  and Frame-only CV `build(const Frame&)`
- Does **not** gate `tests/` (unit lifecycle) or `docs/developer/`

## Related

- `.claude/notes/frame-region-fields.md` — which atoms → Frame fields
- `CLAUDE.md` iron law `no-public-config-bags` — no `*Config` bags either
- Architecture: Frame sole I/O; Context tiers are engine-internal
