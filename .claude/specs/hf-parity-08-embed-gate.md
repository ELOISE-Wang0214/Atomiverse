---
title: Embed / external field hard gate (finite-field dipole, MM bookkeeping)
slug: hf-parity-08-embed-gate
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P1. Run `/mol:spec
  hf-parity-08-embed-gate` to research + promote. science.required.
---

# Embed / external field hard gate

## Summary
Point-charge embedding / external field (`regressions/scf/h2o_embedding_field`)
shipped without a fully verified gate. Pin: finite-field dipole ≡ analytic
dipole; point-charge H₂O vs pinned PySCF `qmmm.mm_charge` with correct MM–MM
bookkeeping. Parent: hf-deferred-backlog P1.

## Domain basis
External field F couples via −μ·F; the finite-field derivative of E vs F must
equal the analytic dipole. QM/MM point charges add nuclear–MM and electron–MM
terms; MM–MM energy is bookkeeping-only (constant to the SCF). Reference: PySCF
`pyscf/qmmm`.

## Design (direction)
- Regression asserting finite-field-dipole ≡ analytic-dipole and the point-charge
  energy vs PySCF `mm_charge`, with explicit MM–MM accounting.

## Files (best-effort)
- `src/cpu/scf/` embed/field impl, `regressions/scf/h2o_embedding_field/` (extend)

## Tasks (provisional)
- [ ] Regression: finite-field dipole ≡ analytic dipole.
- [ ] Point-charge H₂O energy vs pinned PySCF `mm_charge` (MM–MM correct).

## Testing
Regression: both identities on the existing embedding case, ≤ 60 s.

## Out of scope
Polarizable / periodic MM ([[hf-parity-15-mm-embedding]]); embed forces
([[hf-parity-03-environment-forces]]).

## Depends on
Existing embed/field code (in tree).
