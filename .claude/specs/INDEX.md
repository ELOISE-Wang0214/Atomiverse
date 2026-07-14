# Specs index

Breadth-first campaign (2026-07-12 → 2026-07-14): **complete**. Four probes +
one decision, all closed (history in git log; findings in `.claude/NOTES.md`;
the Frame verdict and its research report live durably at
`docs/developer/frame-block-redesign.md`).

Campaign **frame-block** — executing the signed Frame verdict (Frame =
mapping[string, Block], molrs-aligned schema, independent per-backend storage)
via the staged variant-C route. **Step 01 (registry + ColKey + mapping API +
typed meta) is done** (git log). Next: 02+ blockify directory by directory
(torch -> cuda -> cpu/md, minimizer, path -> cpu/potentials), deleting each
member once its directory is clear.

| slug | status | summary |
|---|---|---|

## Standing user directives (2026-07-12/14) — binding on every spec

1. **Frame is the only input/output currency.** No bare parallel scalar arrays.
2. **Frame = mapping[string, Block]** (signed 2026-07-14): molrs-aligned schema
   and keys, independent per-backend Column storage, bookkeeping in ONE typed
   `meta` region (string→MetaValue, extending molrs's string→string; Box rides
   simbox; revised 2026-07-14), dtype follows ATV precision tiers. Resolution
   stops at the launch boundary — in-loop key lookup is forbidden (13.7x) and
   machine-checked (debug resolution-count assertion + LJ ⚡ gate).
3. **N replicas come from an expansion strategy, not from the user.**
