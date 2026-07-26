---
slug: ot-dftu
status: done
created: 2026-07-26
domain: gpw
queue: cp2k-method-loop
depends_on: [ot-uks-o2, gpw-dftu-nio]
unblocks: []
---

# Spec: `ot-dftu` — OT + Dudarev DFT+U (Γ)

## Gap

Γ-UKS **OT** (`use_ot`) and **Dudarev DFT+U** (`plus_u`) are each product paths
(`reg_ot-uks-o2`, `reg_gpw-dftu-nio`). Combining them on the UKS OT path **throws**:

```text
GPW UKS OT: DFT+U is not wired on the OT path (use_ot=false or drop plus_u)
```

The KS builder used by UKS OT (`build_uks_H`) already applies Mulliken Dudarev
PLUS_U energy + potential; only the fail-fast gate blocks the combo. CP2K twins
for DFT+U already use `&OT` + `PLUS_U` (see `tools/cp2k_oracle/inputs/gpw-dftu-nio/`).

## Goal

Ship **Γ-UKS OT + Dudarev DFT+U** as a product path:

- `cpu::gpw::GPWEvaluator` / `GpwDft` with `use_ot=true` and non-empty `plus_u`
  at Γ (no k-mesh) converge to the same fixed point class as DIIS + PLUS_U.
- Offline CP2K pin (OT + PLUS_U) in `regressions/gpw/ot-dftu/`.
- No public `*Config`, no `make_*` factories; Frame I/O; public energy still eV.

## Design

1. **Remove** the UKS OT `do_dftu` throw; keep multi-k + OT rejection (out of
   scope — slug `ot-kpoints`).
2. Rely on existing `internal::apply_dft_plus_u` inside the UKS H/energy rebuild
   (already shared by DIIS UKS and the OT energy / Armijo line search).
3. OT gradient / FULL_SINGLE_INVERSE preconditioner see the DFT+U-augmented
   spin Fock matrices (Hellmann–Feynman consistent with \(E_U[n]\)).
4. Unit: small H₂ UKS + PLUS_U — OT converges, \(E_U > 0\), |E_OT − E_DIIS| tight.
5. Regression: NiO molecule (same twin as `gpw-dftu-nio`) with `use_ot=true`
   vs offline CP2K OT+PLUS_U pin; optional OT/DIIS parity.

## Out of scope

- OT + multi-k (`ot-kpoints`)
- Closed-shell RKS Γ DFT+U (Gamma RKS `build_H_and_energy` still omits PLUS_U;
  multi-k RKS +U and UKS Γ +U remain the product surfaces)
- Analytic forces/stress under OT+DFT+U (DIIS +U stress/forces remain separate)
- Live CP2K in ctest

## Acceptance

See `acceptance.md`. Green when dual gate passes and `acceptance_met=true`.
