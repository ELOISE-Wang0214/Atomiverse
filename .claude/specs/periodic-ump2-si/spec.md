---
slug: periodic-ump2-si
status: in_progress
created: 2026-07-26
domain: gpw / post-hf
queue: cp2k-method-loop
depends_on: [periodic-mp2-si, periodic-ri-mp2-si, gpw-uks-pbe]
unblocks: [periodic-ump2-si]
---

# Spec: `periodic-ump2-si` — Unrestricted GPW MP2 (UKS reference)

## Summary

Ship **Γ-point unrestricted GPW MP2 (UMP2)** on a converged UKS reference via
the existing public `GpwDft::mp2(atoms)` entry (no new factory / Config bag).
Today `mp2` throws on UKS; this slug removes that guard for Γ, implements the
molecular/PySCF UMP2 doubles sum (same-spin αα/ββ + opposite-spin αβ) on
RI-GPW (product) or Poisson (foothold) \((ia|jb)\) tensors, dual-gates against
the closed-shell RMP2 Si₂ first-green band, and leaves multi-k UMP2 as a
follow-on.

## Science

- Pair densities and Coulomb path: same Quickstep RI-GPW / Poisson machinery
  as RKS MP2 (`mo_eri_ovov` / `mo_eri_ovov_ri`), plus cross-spin
  `mo_eri_ovov_ab` / `_ri`.
- Doubles: PySCF / molecular `ump2_energy` convention
  \(E_{\mathrm{SS}}^\sigma=\frac14\sum[(ia|jb)-(ib|ja)]^2/\Delta\),
  \(E_{\mathrm{OS}}=\sum(i_\alpha a_\alpha|j_\beta b_\beta)^2/\Delta\).
- Identity: when \(C^\alpha=C^\beta\) and \(n_\alpha=n_\beta\), UMP2 = RMP2.

## Public API

```cpp
GpwDft m;
m.n_alpha = 4; m.n_beta = 4;  // UKS
m.build(atoms); m.run(atoms);
const Mp2Result& r = m.mp2(atoms);  // UMP2@UKS at Γ
```

No new public types. Parameters remain members (`mp2_ri`, `mp2_n_virt`,
`mp2_os_scale`, `mp2_ss_scale`).

## Out of scope

- Multi-k UMP2 / UGW / URPA (hold behind this foothold)
- Analytic UMP2 forces / stress
- OT + MP2 (still requires diagonalization virtuals)

## Acceptance

See `acceptance.md`. Green ⇒ `acceptance_met=true`.
