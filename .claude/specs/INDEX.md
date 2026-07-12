# Specs index

Breadth-first architecture campaign (2026-07-12). Per user directive: implement
methods broadly to surface architectural defects fast; deepen later. Each probe
targets a **different pillar** and must land a passing out-of-library textbook
case under `regressions/`.

| slug | status | summary |
|---|---|---|
| [regressions-harness](regressions-harness.md) | approved | `regressions/` — out-of-library textbook e2e vs the literature; **gate for all four probes** |
| [manybody-potentials](manybody-potentials.md) | approved | Probe 1 (Component/Neighbor): EAM + SW. Wires up the dead `refresh_full`; first production user of `CommField` |
| [md-constraints](md-constraints.md) | approved | Probe 2 (Phase/Tagging): SETTLE + RATTLE. Closes the `Capability` write-set hole, adds `ComponentRole::Constraint`, gives Ndof one owner |
| [path-driver-neb](path-driver-neb.md) | approved | Probe 3 (Driver taxonomy): CI-NEB + the 6th category. Extracts a reusable `ReplicaSet`; de-welds L-BFGS from the Frame |
| [mc-gcmc-npt](mc-gcmc-npt.md) | approved | Probe 4 (Frame): GCMC + NPT-MC. Lifts variable-N out of `comm/` instead of reinventing it |
| [frame-organization](frame-organization.md) | draft | **决策 spec**：Frame 按后端分类 vs 后端×领域分类。刻意排在四个探针之后 —— 探针就是证据 |

## Ordering

```
regressions-harness            (first — every probe's acceptance depends on it)
      |
      +-- manybody-potentials  ─┐  both touch comm/halo + neighbor/
      +-- md-constraints       ─┘  → SERIALIZE these two, do not run in parallel
      |
      +-- mc-gcmc-npt           (independent — touches frame/atom + mc/ only)
      |
      +-- path-driver-neb       (last probe — largest blast radius: new PhaseKind,
                                 new category, ReplicaSet, LBFGS refactor)
      |
      +-- frame-organization    (DECISION, after all four — consumes their findings)
```

`mc-gcmc-npt` may run concurrently with the manybody/constraints pair.

## Standing user directives (2026-07-12) — binding on every spec above

1. **Frame is the only input/output currency.** No bare `old_xyz` / `new_xyz`
   parallel scalar arrays anywhere. `mc::Transaction` is therefore rewritten as a
   **Frame patch**, not given verbs on top of its raw coordinate vectors.
2. **Frame organization is an open decision** — see `frame-organization`. Do not
   speculatively restructure Frame inside a probe; record what the probe needed and
   feed it to the decision spec.
3. **N replicas come from an expansion strategy, not from the user.** A Path/Ensemble
   driver takes **one Frame** and a `ReplicaExpansion` that expands it to N (NEB =
   endpoint interpolation, PIMD = bead replication, REMD = replication + a temperature
   ladder, umbrella = replication + bias centers), then collapses back to one Frame.
   `Driver::run(frame, ...)` keeps its single-Frame signature.

## What the recon already proved (before a line of code)

- The breaks are **not** where the methods are — they are in `comm/` and `neighbor/`.
  A "one comm pass per step, and it only carries forces" assumption blocks both
  many-body (ghost ρ) and constraints (ghost velocities).
- `Capability` claims to be a write-set but has **no Position/Velocity bit**, so every
  component that actually moves atoms declares `Capability::None`.
- The Driver taxonomy has **one axis (termination)** but "N coupled replicas" is an
  orthogonal second axis — REMD/PIMD/umbrella have no home under it.
- The Frame's "fixed-N assumption" is a **myth**: `comm/` already implements variable N.
- Dead-but-written code to switch on, not rewrite: `CellList::refresh_full`,
  `HaloParticipant::CommField`, the `SamplingDriver` acceptance seam.
