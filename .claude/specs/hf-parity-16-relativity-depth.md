---
title: x2c picture-change properties; DKH; SO-ECP; SOC / complex GHF
slug: hf-parity-16-relativity-depth
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P3. Speculative — "not parity for
  parity"; promote only when a real heavy-element application demands it. Run
  `/mol:spec hf-parity-16-relativity-depth`. science.required.
---

# Relativity depth

## Summary
Depth beyond x2c-1e energy: **picture-change property operators** (correct
properties under x2c), optional **DKH**, **SO-ECP**, lanthanide **4f-in-core**
tables, and **SOC / complex GHF** — each only if a real application demands it.
Parent: hf-deferred-backlog P3.

## Domain basis
Picture-change: property operators must be transformed by the same x2c
decoupling U as the Hamiltonian, else properties are inconsistent. SO-ECP adds
the spin-orbit projector. SOC/complex GHF needs a complex 2-component density.
Reference: Ilias/Saue; Reiher/Wolf DKH.

## Design (direction — research heavy; likely splits)
- Picture-change transform of property integrals; SO-ECP projector; complex GHF
  storage. Each is a distinct child at spec time.

## Files (best-effort)
- `src/cpu/scf/` x2c + ECP + GHF code; `regressions/scf/<case>/` per feature

## Tasks (provisional)
- [ ] Split per feature at spec time; per feature: property/energy vs pinned reference.

## Testing
Per feature: pinned reference regression ≤ 60 s.

## Out of scope
x2c-1e energy gate ([[hf-parity-06-x2c-gate]]); ECP energy/forces
([[hf-parity-01-ecp-forces]]); x2c+ECP double counting (non-goal).

## Depends on
[[hf-parity-06-x2c-gate]]; [[hf-parity-02-openshell-forces]] (GHF).
