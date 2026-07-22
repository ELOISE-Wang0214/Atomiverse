---
title: Static polarizability (electric-field CPHF); optional hyperpolarizability
slug: hf-parity-10-polarizability
status: draft
chain: hf-parity
created: 2026-07-22
note: >
  Roadmap draft from hf-deferred-backlog P2. Run `/mol:spec
  hf-parity-10-polarizability` to research + promote. science.required.
---

# Static polarizability (electric-field CPHF)

## Summary
Add the static dipole polarizability α = −d²E/dF² via electric-field CPHF;
optionally the first hyperpolarizability β. Parent: hf-deferred-backlog P2.
Reuses the CPHF solver already built for the RHF Hessian.

## Domain basis
α_ij = −∂²E/∂F_i∂F_j; solve the CPHF equations with the dipole operator as the
perturbation. β from the 2n+1 rule / third derivative. Reference:
Helgaker/Jørgensen/Olsen response theory; PySCF `pyscf/prop/polarizability`.

## Design (direction)
- Reuse the CPHF machinery with dipole-integral RHS; contract to α; write to
  Frame meta. β optional behind the same solver.

## Files (best-effort)
- `src/cpu/scf/` CPHF code; `regressions/scf/h2o_polarizability/` (new)

## Tasks (provisional)
- [ ] RED: α_zz vs FD of the dipole under a finite field (unit).
- [ ] Electric-field CPHF → α tensor; write to Frame.
- [ ] Regression: H2O α vs pinned PySCF.
- [ ] (Optional) β via 2n+1.

## Testing
Unit FD-vs-CPHF α_zz (seconds). Regression α vs pinned PySCF ≤ 60 s.

## Out of scope
Frequency-dependent (dynamic) polarizability ([[hf-parity-12-tddft]]).

## Depends on
Existing RHF CPHF solver.
