---
title: Forces through environment members (C-PCM, x2c, embed/field, RI-JK, CFMM-J)
slug: hf-parity-03-environment-forces
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft carved from hf-deferred-backlog P0. Likely SPLITS into
  per-member children at /mol:spec time (each environment term has its own
  gradient theory + oracle). Run `/mol:spec hf-parity-03-environment-forces`.
  science.required.
---

# Forces through environment members

## Summary
Several environment members currently throw when analytic forces are requested:
**C-PCM**, **x2c**, **point-charge embed / external field**, **RI-JK**, and
**CFMM-J**. Add their nuclear-gradient contributions so geometry optimization
works with each. Parent: hf-deferred-backlog P0. This is an umbrella — expect a
split into one child per member (their gradient math + oracles differ).

## Domain basis (per member)
- **RI-JK forces**: 3-index derivative integrals + 2-index metric derivative
  (DF gradient). Reference: PySCF `pyscf/df/grad`.
- **Embed / external field forces**: point-charge–nucleus and
  point-charge–electron force terms + MM–MM bookkeeping. Reference: PySCF
  `qmmm.mm_charge` gradient. (Ties to [[hf-parity-08-embed-gate]].)
- **C-PCM forces**: cavity-surface + solvent-reaction-field gradient, incl.
  cavity-discontinuity caveat. Reference: PySCF `pyscf/solvent/grad`.
- **x2c forces**: h^X2C nuclear derivative (picture-change on the gradient);
  today analytic forces with x2c throw. Reference: PySCF `x2c` grad.
- **CFMM-J forces**: derivative of the FMM/near-field split J.

## Design (direction — research before impl)
- One gradient contribution per member, gated on the member being active;
  remove each throw only after its term is validated. No shared polymorphic
  "environment gradient".

## Files (best-effort)
- `src/cpu/scf/` environment member impls (embed/pcm/x2c/ri-jk), `gradient.cpp`
- `regressions/scf/<member>_grad/` (one per member)

## Tasks (provisional)
- [ ] Split into per-member children at spec time.
- [ ] Per member: RED FD check → analytic term → drop throw → pinned regression.

## Testing
Per member: FD unit check + one pinned-PySCF regression ≤ 60 s.

## Out of scope
CFMM-K (permanent non-goal). Advanced PCM families
([[hf-parity-14-advanced-solvent]]).

## Depends on
[[hf-parity-08-embed-gate]] / [[hf-parity-07-cpcm-gate]] /
[[hf-parity-06-x2c-gate]] (energy hard gates should land first).
