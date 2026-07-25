---
spec: periodic-mp2-kpoints
created: 2026-07-25
criteria:
  - id: ac-001
    summary: Multi-k SCF stores complex Bloch C_k
    type: code
    pass_when: |
      GpwResult carries per-k complex MO coefficients after a multi-k
      converged run; store_k_payload no longer discards C.
    status: met
  - id: ac-002
    summary: GpwDft::mp2 accepts multi-k references
    type: runtime
    pass_when: |
      mp2() no longer throws invalid_argument solely because n_kpts > 1;
      still throws for UKS, OT, unconverged, and vanishing gap.
    status: met
  - id: ac-003
    summary: Kernel Γ recovery (unit)
    type: runtime
    pass_when: |
      test_cpu_gpw_mp2: mo_on_grid_k at κ=0 ≡ real mo_on_grid; complex
      OS/SS ≡ real mp2_rhf_energy. (Full SCF k-Γ vs real-Γ blocked by
      k-path level_shift ban on artificial small-gap cells.)
    status: met
  - id: ac-004
    summary: Multi-k regression gates on bct-Si4 2×1×1
    type: runtime
    pass_when: |
      regressions/gpw/periodic-mp2-kpoints builds under ATV_BUILD_GPW,
      never launches CP2K, asserts e_corr < 0, finite, OS+SS and
      e_total identities, and first-green band (or CP2K pin when present).
    status: met
  - id: ac-005
    summary: Public API hygiene
    type: code
    pass_when: |
      No public *Config, no make_*/build_* factory for multi-k MP2;
      parameters remain GpwDft members.
    status: met
  - id: ac-006
    summary: Offline oracle input shipped
    type: docs
    pass_when: |
      tools/cp2k_oracle/inputs/periodic-mp2-kpoints.inp + REFERENCE.md
      document twin settings for re-pin.
    status: met
out_of_scope:
  - UMP2 k-points
  - RI_AUX Gaussian tables
  - MP2 forces/stress
  - Live CP2K in ctest
  - RPA/GW
---

# Acceptance — `periodic-mp2-kpoints`

**Done** means multi-k closed-shell GPW MP2 correlation energy is available on
`GpwDft::mp2` after a multi-k `run()`, with Γ recovery identity and a
regression deck under `regressions/gpw/periodic-mp2-kpoints/`.
