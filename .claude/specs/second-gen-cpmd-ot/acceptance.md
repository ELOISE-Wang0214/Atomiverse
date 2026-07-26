---
spec: second-gen-cpmd-ot
created: 2026-07-26
criteria:
  - id: ac-001
    summary: build() accepts second_gen_cp + use_ot (RKS Γ)
    type: runtime
    pass_when: |
      test_gpw_second_gen_cpmd BuildAcceptsWithOt passes; BuildRejectsWithOt removed
      or inverted.
    status: met
  - id: ac-002
    summary: GpwDft::run accepts guess_P with use_ot via natural orbitals
    type: code
    pass_when: |
      No throw on use_ot + guess_P; ASPC density is purified to C0 for OT
      (not silently ignored).
    status: met
  - id: ac-003
    summary: Production OT corrector respects second_gen_cp_max_iter
    type: runtime
    pass_when: |
      Regression asserts last_result().iters ≤ corrector budget after ASPC is warm.
    status: met
  - id: ac-004
    summary: NVE energy-drift envelope on H2O
    type: runtime
    pass_when: |
      reg_second-gen-cpmd-ot: relative |E−E0|/scale ≤ 0.10 over short production,
      PE < −100 eV after seed, O–H < 1.8 Å.
    status: met
  - id: ac-005
    summary: Offline CP2K twin input + REFERENCE.md
    type: docs
    pass_when: |
      tools/cp2k_oracle/inputs/second-gen-cpmd-ot/ has ASPC+OT+MAX_SCF twin;
      REFERENCE.md cites Kühne 2007 + OT path + pin status.
    status: met
  - id: ac-006
    summary: No public Config / factory
    type: code
    pass_when: |
      Members-only API on GPWEvaluator; no make_*/build_* for this feature.
    status: met
out_of_scope:
  - UKS OT SGCP
  - multi-k / smear / PLUS_U with SGCP+OT
  - live CP2K in ctest
---

# Acceptance — `second-gen-cpmd-ot`

**Done** means closed-shell Γ GPW second-gen CP can use OT as the partial
SCF corrector, with a short H₂O NVE regression and offline CP2K twin sketch.
