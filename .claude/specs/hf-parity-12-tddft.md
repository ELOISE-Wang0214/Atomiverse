---
title: UHF-TDHF / spin-flip TDA; TDDFT (RKS/UKS A/B)
slug: hf-parity-12-tddft
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P2. Run `/mol:spec
  hf-parity-12-tddft` to research + promote. science.required.
---

# UHF-TDHF / spin-flip TDA; TDDFT

## Summary
RHF TDHF/TDA is in tree (`regressions/scf/ch2o_tdhf`). Add **UHF-TDHF** and
**spin-flip TDA**, and **TDDFT** (RKS/UKS A/B matrices with the XC kernel).
Parent: hf-deferred-backlog P2.

## Domain basis
Casida A/B matrices; TDA = B→0. UHF/UKS use per-spin blocks; spin-flip couples
α→β. TDDFT adds f_xc to the coupling matrix. Reference: Casida; PySCF
`pyscf/tdscf/{uhf,rks,uks}`.

## Design (direction)
- Extend the existing TDHF Davidson eigensolver to per-spin blocks (UHF/UKS) and
  add the XC kernel term for TDDFT; spin-flip variant of TDA.

## Files (best-effort)
- `src/cpu/scf/` TDHF/TDA code + XC kernel; `regressions/scf/{uhf_tda,tddft}/` (new)

## Tasks (provisional)
- [ ] RED: smallest excitation energy vs pinned PySCF (unit-ish, capped iters).
- [ ] UHF-TDHF / spin-flip TDA blocks.
- [ ] TDDFT A/B with f_xc.
- [ ] Regressions vs pinned PySCF excitation energies.

## Testing
Regression: excitation energies vs pinned PySCF ≤ 60 s (small molecule, few roots).

## Out of scope
Excited-state gradients / NACs (only if a concrete consumer needs them —
separate draft when that appears).

## Depends on
Existing RHF TDHF/TDA; the KS XC kernel from [[hf-parity-11-uhf-ks-hessian]].
