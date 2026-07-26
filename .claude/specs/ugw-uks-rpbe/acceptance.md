---
spec: ugw-uks-rpbe
created: 2026-07-26
criteria:
  - id: ac-001
    summary: eval_xc_grid_spin(RPBE) product, no throw
    type: code
    pass_when: |
      src/cpu/gpw/g0w0.cpp maps XcKind::RPBE to xc_grid_rpbe_spin; the
      "UKS RPBE spin XC for UGW is not yet product" throw string is gone.
    status: met
  - id: ac-002
    summary: Γ UGW@UKS-RPBE regression dual-gates RKS G0W0@RPBE
    type: runtime
    pass_when: |
      regressions/gpw/ugw-uks-rpbe builds as reg_ugw-uks-rpbe, runs ≤ 60 s,
      never launches CP2K; |gap_UGW − gap_RKS| ≤ 0.05 eV on the Si₂ deck
      with nα=nβ=4; internal SCF+Newton gates pass.
    status: met
  - id: ac-003
    summary: RPBE v_xc is active (differs from PBE UGW pin)
    type: runtime
    pass_when: |
      Regression asserts |gap_UGW@RPBE − 5.719545 eV| ≥ 1e-3 eV (PBE UGW
      first-green of periodic-ugw-si) so the spin XC path is not silently
      still PBE.
    status: met
  - id: ac-004
    summary: Offline CP2K twin documented
    type: docs
    pass_when: |
      REFERENCE.md + tools/cp2k_oracle/inputs twin cite CP2K version,
      functional form (GGA_X_RPBE + GGA_C_PBE), cell, basis/PP, cutoff;
      pin strategy matches periodic-ugw-si (dual-gate + internal; no
      invented RI-GW absolute QP digits).
    status: met
  - id: ac-005
    summary: API iron laws
    type: code
    pass_when: |
      No public *Config, no make_*/build_* factories; GpwDft members only.
    status: met
  - id: ac-006
    summary: perf_exempt or within_10x
    type: runtime
    pass_when: |
      REFERENCE.md records perf_exempt with fair reason (spectral ≠
      CP2K RI-G0W0) or T_atv/T_oracle ≤ 10.
    status: met
out_of_scope:
  - Absolute CP2K RI-G0W0 QP digit match
  - Multi-k-only residual (Γ product closes XcKind throw for both paths)
  - CUDA
  - Live oracle in ctest
---

# Acceptance — `ugw-uks-rpbe`

**Done** means spectral UGW on UKS-RPBE no longer throws, is dual-gated
against RKS G0W0@RPBE on the same deck, documents the offline CP2K twin,
and keeps public API iron laws.

## AC-001 — product wire

The residual is a one-line XcKind hole, not a new algorithm.

## AC-002 — dual-gate

UGW ≡ RGW when α=β on closed-shell occupations.

## AC-003 — functional integrity

Guards silent PBE fallback.

## AC-004 — oracle contract

Offline pins only; spectral algorithm honesty.

## AC-005 — OOP / no Config

Constructor and members only.

## AC-006 — perf

Same exempt class as other periodic GW product rows.
