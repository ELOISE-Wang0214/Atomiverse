---
title: Population analysis depth (NPA/NBO, IAO/IBO, RESP)
slug: hf-parity-17-population-analysis
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P3. Run `/mol:spec
  hf-parity-17-population-analysis` to research + promote. science.required.
---

# Population analysis depth

## Summary
Mulliken / Löwdin / MK / CHELPG exist; add **NPA/NBO**, **IAO/IBO**, and
**RESP**. Parent: hf-deferred-backlog P3. These are property outputs written to
Frame (never a "charge driver").

## Domain basis
NPA: natural atomic orbitals from the density's per-atom blocks; NBO builds
localized bonds. IAO/IBO: intrinsic atomic/bond orbitals (Knizia). RESP:
restrained fit of charges to the ESP (two-stage). Reference: Reed/Weinhold NBO;
Knizia IAO; Bayly RESP.

## Design (direction)
- Each is a post-SCF analysis over the converged density; output charges/orbitals
  to Frame meta. IAO/IBO share a localization core.

## Files (best-effort)
- `src/cpu/scf/` or `observer/` population code
- `regressions/scf/{npa,iao,resp}_charges/` (new)

## Tasks (provisional)
- [ ] RED: per-method charges vs pinned reference on a tiny molecule (unit-ish).
- [ ] NPA/NBO; IAO/IBO; RESP.
- [ ] Regressions vs pinned PySCF/reference charges.

## Testing
Regression: charges vs pinned reference ≤ 60 s per method.

## Out of scope
Existing Mulliken/Löwdin/MK/CHELPG (in tree).

## Depends on
Independent; reuses the converged density.
