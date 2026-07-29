# Prefer composition over inheritance

**Status:** binding (2026-07-29). Router: `CLAUDE.md` iron law
`prefer-composition`.

## Rule

Workflows are **outer Driver + registered Components/evaluators**, not deep
inheritance or god drivers.

## Patterns

| Need | Compose | Do not |
|------|---------|--------|
| Classical MD | `MD` + pair/bonded evaluators | subclass `MD` per force field |
| AIMD | `MD` + `set_evaluator(RHF\|CarParrinello\|…)` | `AimdDriver : MD` |
| QM/MM | `MD` + MM evaluators + `Qmmm*` connection evaluator + `set_engine(QM)` | `QmmmDriver`, dual `MD::run` |
| NEB / path | `PathDriver` + `NebCoupling` + nested optimizers | bake NEB into `Minimize` |

## When inheritance / a new preset is OK

- **New loop shape** (phases + termination topology), e.g. coupled-replica
  `PathDriver`.
- **Kind base** one level under `Component` (`Potential`, `HF`, …) — depth ≤ 2.
- **Concrete method** types (`RHF`, `PairLJCut`) — domain types, not workflow
  taxonomies.

## Anti-patterns

- Category driver bases (deleted: `DynamicsDriver`, …)
- One class + mode enum for mutually different Hamiltonians (god)
- Region/partition API on the driver instead of Frame fields (see
  `frame-region-fields`)
