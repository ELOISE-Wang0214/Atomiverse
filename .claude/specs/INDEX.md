# Specs index

Breadth-first campaign (2026-07-12 → 2026-07-14): **complete**. Four probes +
one decision, all closed (history in git log; findings in `.claude/NOTES.md`;
the Frame verdict and its research report live durably at
`docs/developer/frame-block-redesign.md`).

Next campaign: **frame-block** — execute the signed Frame verdict
(Frame = mapping[string, Block], molrs-aligned schema, independent per-backend
storage) via the staged variant-C route.

| slug | status | summary |
|---|---|---|

## Standing user directives (2026-07-12/14) — binding on every spec

1. **Frame is the only input/output currency.** No bare parallel scalar arrays.
2. **Frame = mapping[string, Block]** (signed 2026-07-14): molrs-aligned schema
   and keys, independent per-backend Column storage, fixed typed attrs
   (box/scalars/pop/tag-map), dtype follows ATV precision tiers. Resolution
   stops at the launch boundary — in-loop key lookup is forbidden (13.7x) and
   machine-checked (debug resolution-count assertion + LJ ⚡ gate).
3. **N replicas come from an expansion strategy, not from the user.**
