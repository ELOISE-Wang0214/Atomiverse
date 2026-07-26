---
spec: esmm-cpmd-03-observability
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Electronic energy meta keys registered and written each CPMD launch
    type: runtime
    pass_when: |
      kEnergyElectronKinetic and kEnergyElectronThermostat are registered
      serializable in core/frame_schema.h; after any GPW CPMD launch
      (Γ RKS, Γ UKS, multi-k) both exist on frame.meta, are finite, and
      T_e matches last_electronic_ke_ev() bit-for-bit (test in
      test_gpw_cpmd.cpp).
    status: pending
  - id: ac-002
    summary: CpmdConserved assembles E_cons on kConserved
    type: runtime
    pass_when: |
      tests/cpu/dynamics/test_cpmd_conserved.cpp shows kConserved ==
      kEnergy + kKinetic + kEnergyElectronKinetic + kEnergyElectronThermostat
      after launch; default binding is md::Phase::end_of_step @ -90.
    status: pending
  - id: ac-003
    summary: Missing electronic keys fail fast
    type: runtime
    pass_when: |
      CpmdConserved::launch throws a typed exception when the electronic
      keys are absent (no CPMD evaluator ran); EXPECT_THROW test present.
    status: pending
  - id: ac-004
    summary: cpmd-water drift pin flows through Frame meta
    type: runtime
    pass_when: |
      regressions/dynamics/cpmd-water wires CpmdConserved and pins the
      relative E_cons drift envelope read from frame.meta[kConserved];
      tolerance unchanged from the shipping pin.
    status: pending
  - id: ac-005
    summary: Multi-writer hazard documented
    type: docs
    pass_when: |
      CpmdConserved Doxygen and the cpmd-water tutorial page state the
      do-not-combine rule with kConserved-owning integrators (NPT, CSVR)
      and include the user-space StopWhenDrift criterion snippet.
    status: pending
out_of_scope:
  - SCF-side CPMD key writes (esmm-cpmd-04)
  - kConserved multi-writer runtime detection
  - New Thermo columns
---

# Acceptance — `esmm-cpmd-03-observability`

**Done** means a Thermo line shows Car–Parrinello E_cons drift and the
fictitious T_e is a Frame-visible adiabaticity monitor, with the cpmd-water
regression pinning drift through the public meta channel.

## AC-001 — the adiabaticity channel

T_e on Frame is what tells a user whether μ/dt are safe; readbacks alone
cannot feed observers.

## AC-002 — one assembler, correct phase

After KineticEnergy (@ -100), before Reduce/Publish/Thermo — the only spot
where every term is current.

## AC-003 — no fallback

An E_cons of "whatever was lying around" is worse than a throw.

## AC-004 — regression witness

The shipping drift pin now proves the public observability path.

## AC-005 — honest limits

Explicit wiring + documented multi-writer hazard, per the decided design.
