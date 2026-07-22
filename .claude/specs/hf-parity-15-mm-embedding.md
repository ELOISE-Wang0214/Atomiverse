---
title: Polarizable / periodic MM embedding; Ewald-corrected external charges
slug: hf-parity-15-mm-embedding
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P3. Speculative — promote when a real
  QM/MM application demands it. Run `/mol:spec hf-parity-15-mm-embedding`.
  science.required.
---

# Polarizable / periodic MM embedding

## Summary
Extend point-charge embedding to **polarizable MM** (induced dipoles) and
**periodic MM** with **Ewald-corrected** external charges. Parent:
hf-deferred-backlog P3.

## Domain basis
Polarizable MM adds self-consistent induced dipoles coupled to the QM density;
periodic embedding needs Ewald summation of the external-charge potential (ties
to the engine's existing `longrange/` Ewald/PME). Reference: AMOEBA/QM-MMPol;
Ewald QM/MM.

## Design (direction)
- Add an induced-dipole self-consistency loop coupled to SCF; reuse `longrange/`
  Ewald for the periodic external potential.

## Files (best-effort)
- `src/cpu/scf/` embed code, `longrange/` (Ewald reuse)
- `regressions/scf/mm_polarizable/` (new)

## Tasks (provisional)
- [ ] Polarizable-MM induced-dipole coupling.
- [ ] Ewald-corrected periodic external charges.
- [ ] Regression vs pinned reference.

## Testing
Regression vs pinned reference, ≤ 60 s.

## Out of scope
Non-polarizable point-charge embed ([[hf-parity-08-embed-gate]]).

## Depends on
[[hf-parity-08-embed-gate]]; existing `longrange/` Ewald/PME.
