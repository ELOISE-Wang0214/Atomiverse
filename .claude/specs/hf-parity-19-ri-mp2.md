---
title: RI-MP2 product path / DF-ao2mo completeness
slug: hf-parity-19-ri-mp2
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P5. Run `/mol:spec hf-parity-19-ri-mp2`
  to research + promote. science.required.
---

# RI-MP2 / DF-ao2mo completeness

## Summary
An MP2 teaching path exists; make **RI-MP2** a product path and complete
**DF-ao2mo** (the density-fitted AO→MO integral transform). Parent:
hf-deferred-backlog P5.

## Domain basis
RI-MP2 approximates (ia|jb) via 3-index DF integrals (ia|P)(P|Q)^{-1}(Q|jb),
turning the ao2mo bottleneck into DF contractions. Reference: Feyereisen RI-MP2;
Weigend def2 JK/RI-C aux. The DF aux basis (`Def2JKFIT`) is already in tree.

## Design (direction)
- Build the DF 3-index (ia|P) transform; assemble the RI-MP2 energy; reuse the
  spherical aux path already present.

## Files (best-effort)
- `src/cpu/scf/` MP2 + density-fitting code (`density_fitting.cpp`)
- `regressions/scf/rimp2_energy/` (new)

## Tasks (provisional)
- [ ] RED: DF-ao2mo transform vs exact ao2mo on a tiny system (unit).
- [ ] RI-MP2 correlation energy assembly.
- [ ] Regression: RI-MP2 energy vs pinned PySCF `mp2.RIMP2`.

## Testing
Unit: DF-ao2mo vs exact (seconds). Regression: RI-MP2 energy vs pinned PySCF ≤ 60 s.

## Out of scope
CCSD/CCSD(T)/CASSCF (only under a dedicated science campaign — separate draft
when that campaign starts).

## Depends on
Existing MP2 path + DF aux basis.
