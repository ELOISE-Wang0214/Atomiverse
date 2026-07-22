---
title: Docs/Doxygen honesty sweep from the closed HF campaign
slug: hf-parity-09-docs-honesty
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P1. Run `/mol:spec
  hf-parity-09-docs-honesty` to promote. type: docs (owed to a human reviewer /
  --manual close, not /mol:bench).
---

# Docs/Doxygen honesty sweep

## Summary
Sweep any remaining docs/Doxygen stubs from the closed `hf-capability-gaps-*`
chain that still claim **Handler** language, the old **`SCFSolver`** user path,
or are missing **Results** tables. Parent: hf-deferred-backlog P1. Pure docs —
no code behavior change.

## Domain basis
CLAUDE.md is the source of truth: Termination is a Criterion (no Handler);
concrete methods `RHF`/`UHF`/… are the user path (SCFSolver is transitional
debt); every delivered feature ships a Results/regression provenance table.

## Design (direction)
- Grep the tutorial + Doxygen for `Handler`, `SCFSolver` as a user API, and
  feature sections lacking a Results table; rewrite to the current model.

## Files (best-effort)
- `docs/tutorial/ch*.md`, public-header Doxygen under `include/atomiverse/cpu/scf/`

## Tasks (provisional)
- [ ] Grep + inventory dishonest/stale doc claims.
- [ ] Rewrite Handler / SCFSolver-user-path language.
- [ ] Add missing Results tables.

## Testing
Docs-only; verified by review. No ctest change (no e2e in tests/).

## Out of scope
Retiring the `SCFSolver` code itself ([[hf-parity-21-retire-scf-debt]]).

## Depends on
Independent.
