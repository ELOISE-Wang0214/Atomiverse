---
title: Advanced solvent (IEF-PCM / SS(V)PE / SMD, ddCOSMO, surface switching)
slug: hf-parity-14-advanced-solvent
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P3. Speculative — promote only when
  geometry+PCM becomes common. Run `/mol:spec hf-parity-14-advanced-solvent`.
  science.required.
---

# Advanced solvent models

## Summary
Beyond C-PCM: **IEF-PCM / SS(V)PE / SMD-class** (with non-electrostatic terms)
and optionally a linear-scaling **ddCOSMO** solver; plus **surface switching**
(SWIG / York–Karplus) if geometry+PCM optimization becomes common. Parent:
hf-deferred-backlog P3.

## Domain basis
IEF-PCM uses the full dielectric IEF operator (vs C-PCM's conductor
approximation); SMD adds parameterized non-electrostatic (CDS) terms. ddCOSMO
is a domain-decomposition linear-scaling COSMO. Switching functions make the
cavity surface differentiable for gradients. Reference: Tomasi PCM review;
Marenich/Cramer/Truhlar SMD; Cancès/Stamm ddCOSMO.

## Design (direction — research heavy)
- Each model is a distinct surface-charge operator; likely a per-model child at
  spec time. Non-electrostatic terms are additive corrections.

## Files (best-effort)
- `src/cpu/scf/` solvent code; `regressions/scf/<model>/` (one per model)

## Tasks (provisional)
- [ ] Split per model at spec time; per model: energy vs pinned PySCF/reference.

## Testing
Per model: ΔG_solv vs pinned reference, ≤ 60 s.

## Out of scope
C-PCM itself ([[hf-parity-07-cpcm-gate]]); solvent forces
([[hf-parity-03-environment-forces]]).

## Depends on
[[hf-parity-07-cpcm-gate]] (C-PCM baseline first).
