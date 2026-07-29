---
spec: dftb2-kpoints
created: 2026-07-29
criteria:
  - id: ac-001
    summary: Public k_mesh / k_time_reversal members (no Config / factories)
    type: code
    pass_when: |
      DFTB2 exposes std::array<int,3> k_mesh and bool k_time_reversal as public
      members; no *Config bag; no make_*/build_* for k-sampling.
    status: met
  - id: ac-002
    summary: Default k_mesh={1,1,1} reproduces real Γ energy
    type: runtime
    pass_when: |
      On the C-chain cell, k_mesh={1,1,1} energy matches the existing Γ path
      within 1e-8 eV (unit test).
    status: met
  - id: ac-003
    summary: Multi-k energy differs from Γ on the 1-D chain
    type: runtime
    pass_when: |
      k_mesh={4,1,1} |E_k − E_Γ| > 0.01 eV on the C-chain geometry (unit test
      and regression).
    status: met
  - id: ac-004
    summary: Hard energy pin in regression (ATV Hamiltonian twin)
    type: runtime
    pass_when: |
      reg_dftb2-kpoints |E − E_ref| ≤ 0.05 eV with REFERENCE.md pin;
      oracle: dftb+ documented; never shells to DFTB+.
    status: met
  - id: ac-005
    summary: Offline DFTB+ twin inputs land under tools/tb_oracle
    type: docs
    pass_when: |
      tools/tb_oracle/inputs/dftbplus/dftb2-kpoints/ has geo + HSD notes for
      KPointsAndWeights SupercellFolding / MP mesh.
    status: met
  - id: ac-006
    summary: Multi-k guards (forces/open-shell/ThreeOb throw)
    type: runtime
    pass_when: |
      with_forces or open shell or ThreeOb + multi-k throws a clear error;
      unit test covers at least with_forces multi-k.
    status: met
out_of_scope:
  - Analytic multi-k forces / stress
  - Metallic smearing
  - Live DFTB+ in ctest
  - CUDA
---

# Acceptance — `dftb2-kpoints`

**Done** means closed-shell MioSketch SCC-DFTB2 can sample a Monkhorst–Pack
k-mesh beyond Γ on the public Frame path, with a hard energy pin for one
1-D chain case, without breaking Γ / molecular pins.
