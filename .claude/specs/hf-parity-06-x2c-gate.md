---
title: x2c-1e hard gate (pVp vs PySCF, decoupling identities, Ar energy)
slug: hf-parity-06-x2c-gate
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P1. Run `/mol:spec
  hf-parity-06-x2c-gate` to research + promote. science.required.
---

# x2c-1e hard gate

## Summary
x2c-1e (`regressions/scf/x2c_argon`) shipped without a fully verified gate. Pin:
pVp matrix vs PySCF `int1e_pnucp`; the X/R decoupling identities; an Ar (or
similar) energy vs PySCF `x2c1e` (~1e-5 Ha); and an honest ch20 relativity
statement (no picture-change on properties yet). Parent: hf-deferred-backlog P1.

## Domain basis
x2c-1e builds the 4→2 component decoupling from the one-electron Dirac
Hamiltonian; the spin-free pVp integral (`IntegralPack::W`, `int1e_pnucp`
convention) is the key new integral. Reference: Ilias/Saue x2c; PySCF `x2c`.

## Design (direction)
- Unit-pin pVp vs PySCF; assert the decoupling identities; regression Ar energy
  vs `x2c1e`. Fix ch20 wording (picture-change absent on properties).

## Files (best-effort)
- `src/cpu/scf/` x2c impl (pVp/W path), `regressions/scf/x2c_argon/` (extend),
  `docs/tutorial/ch20*.md`

## Tasks (provisional)
- [ ] Unit: pVp/W vs pinned PySCF `int1e_pnucp`.
- [ ] Decoupling identity checks.
- [ ] Regression: Ar x2c1e energy vs pinned PySCF ~1e-5 Ha.
- [ ] Honest ch20 relativity section.

## Testing
Unit pVp pin (seconds). Ar energy regression ≤ 60 s.

## Out of scope
x2c+ECP (mutual exclusion, non-goal); picture-change properties / DKH / SOC
([[hf-parity-16-relativity-depth]]); x2c forces
([[hf-parity-03-environment-forces]]).

## Depends on
Existing x2c-1e code (in tree).
