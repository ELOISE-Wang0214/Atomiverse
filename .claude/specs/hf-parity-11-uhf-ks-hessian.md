---
title: UHF CPHF / analytic Hessian; RKS/UKS second derivatives (XC kernel)
slug: hf-parity-11-uhf-ks-hessian
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P2. Run `/mol:spec
  hf-parity-11-uhf-ks-hessian` to research + promote. science.required.
---

# UHF CPHF / analytic Hessian; RKS/UKS second derivatives

## Summary
Extend the analytic Hessian (RHF only today) to **UHF** (per-spin CPHF) and to
**RKS/UKS** (adds the XC second-derivative kernel f_xc on the grid). Parent:
hf-deferred-backlog P2.

## Domain basis
UHF Hessian: per-spin CPHF response. KS Hessian: the XC contribution needs the
second functional derivative f_xc = δ²E_xc/δρ² collocated on the grid. The
current `IntegralSecondDerivPack` throws for spherical (see the shipped
spherical-1e work) — spherical KS Hessian composes with that. Reference: PySCF
`pyscf/hessian/{uhf,rks,uks}`.

## Design (direction)
- Add per-spin CPHF to the Hessian; add the f_xc grid term for KS. Reuse the
  c2s derivative sandwich for spherical (extend to 2nd derivatives).

## Files (best-effort)
- `src/cpu/scf/` Hessian/CPHF + XC grid; `include/.../internal/integral_derivatives.h`
- `regressions/scf/{uhf_freq,uks_freq}/` (new)

## Tasks (provisional)
- [ ] RED: FD-vs-analytic UHF/UKS Hessian on tiny systems.
- [ ] UHF per-spin CPHF Hessian.
- [ ] KS f_xc second-derivative grid term.
- [ ] Spherical 2nd-derivative c2s sandwich (lift the IntegralSecondDerivPack throw).
- [ ] Regressions vs pinned PySCF frequencies.

## Testing
Unit FD Hessian checks (seconds). UHF + UKS frequency regressions ≤ 60 s.

## Out of scope
IR intensities ([[hf-parity-04-ir-intensities]]); excited-state Hessians.

## Depends on
Existing RHF Hessian; the shipped spherical-1e derivative path (for the 2nd-deriv
c2s sandwich).
