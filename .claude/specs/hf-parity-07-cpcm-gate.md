---
title: C-PCM hard gate (solvation free energy vs PySCF + honest docs)
slug: hf-parity-07-cpcm-gate
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P1. Run `/mol:spec
  hf-parity-07-cpcm-gate` to research + promote. science.required.
---

# C-PCM hard gate

## Summary
C-PCM (`regressions/scf/h2o_rhf_631gs_cpcm`) shipped without a fully verified
gate. Pin: solvation free energy vs PySCF (~0.1 kcal/mol class); document the
absence of non-electrostatic terms and the cavity discontinuity. Parent:
hf-deferred-backlog P1.

## Domain basis
C-PCM: apparent surface charges on a molecular cavity screen the solute; ΔG_solv
= E(solvent) − E(gas). Non-electrostatic (cavitation/dispersion/repulsion) terms
are absent by design. Reference: Cossi/Barone C-PCM; PySCF `pyscf/solvent/pcm`.

## Design (direction)
- Regression pinning ΔG_solv vs PySCF C-PCM; REFERENCE.md documents the
  electrostatic-only scope + cavity-discontinuity caveat.

## Files (best-effort)
- `src/cpu/scf/` C-PCM impl, `regressions/scf/h2o_rhf_631gs_cpcm/` (extend)

## Tasks (provisional)
- [ ] Regression: H2O ΔG_solv vs pinned PySCF C-PCM (~0.1 kcal/mol).
- [ ] Document non-electrostatic absence + cavity discontinuity.

## Testing
Regression: ΔG_solv vs pinned PySCF, ≤ 60 s.

## Out of scope
IEF-PCM / SS(V)PE / SMD / ddCOSMO ([[hf-parity-14-advanced-solvent]]);
C-PCM forces ([[hf-parity-03-environment-forces]]).

## Depends on
Existing C-PCM code (in tree).
