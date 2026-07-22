---
title: CFMM-J hard gate (energy identity + crossover + honest ch19)
slug: hf-parity-05-cfmm-j-gate
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P1. Code is present; the science/docs
  ledger was not fully verified. Run `/mol:spec hf-parity-05-cfmm-j-gate` to
  research + promote. Re-open only if values fail or docs are dishonest.
---

# CFMM-J hard gate

## Summary
The CFMM-J path (FMM kernel + CFMM-J hook, `regressions/scf/c10h22_cfmm_j`)
shipped without a fully verified hard gate. Pin: elementwise CFMM-J vs Schwarz
direct J; converged-RHF energy identity (~1e-7 Ha); an alkane crossover table;
and an honest ch19 statement ("PySCF has no molecular FMM"). Parent:
hf-deferred-backlog P1.

## Domain basis
CFMM splits Coulomb into near-field (direct) + far-field (multipole). The
converged RHF energy must be independent of the J engine to the CFMM tolerance.
Reference: White/Head-Gordon CFMM; Strain/Scuseria.

## Design (direction)
- Add a hard-gate regression asserting the J-engine-invariance of the energy
  and an elementwise J comparison; capture the alkane crossover as a small table
  in REFERENCE.md. Fix the tutorial ch19 wording.

## Files (best-effort)
- `regressions/scf/c10h22_cfmm_j/` (extend), `docs/tutorial/ch19*.md`

## Tasks (provisional)
- [ ] Elementwise CFMM-J vs Schwarz direct J unit/regression check.
- [ ] RHF energy identity CFMM vs direct ≤ ~1e-7 Ha.
- [ ] Alkane crossover table in REFERENCE.md.
- [ ] Honest ch19 "no molecular FMM in PySCF" statement.

## Testing
Regression: energy identity + elementwise J on the c10h22 case, ≤ 60 s.

## Out of scope
CFMM-K (permanent non-goal — K stays direct/RI).

## Depends on
Existing CFMM-J code (in tree).
