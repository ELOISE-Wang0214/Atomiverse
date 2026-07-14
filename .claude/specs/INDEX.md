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
| [frame-block-01-registry](frame-block-01-registry.md) | draft | 变体 C 第一步：字段注册表 + ColKey + mapping API，成员不删，热路径零改动；消灭 7+ 处手抄字段清单 |

## Standing user directives (2026-07-12/14) — binding on every spec

1. **Frame is the only input/output currency.** No bare parallel scalar arrays.
2. **Frame = mapping[string, Block]** (signed 2026-07-14): molrs-aligned schema
   and keys, independent per-backend Column storage, bookkeeping in ONE typed
   `meta` region (string→MetaValue, extending molrs's string→string; Box rides
   simbox; revised 2026-07-14), dtype follows ATV precision tiers. Resolution
   stops at the launch boundary — in-loop key lookup is forbidden (13.7x) and
   machine-checked (debug resolution-count assertion + LJ ⚡ gate).
3. **N replicas come from an expansion strategy, not from the user.**
