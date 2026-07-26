---
spec: periodic-ump2-kpoints
created: 2026-07-26
criteria:
  - id: ac-001
    summary: Multi-k UKS SCF stores full Bloch C_k and C_k_beta
    type: code
    pass_when: |
      GpwResult carries per-k full complex MO coefficients for α (C_k)
      and β (C_k_beta) after a multi-k UKS converged run, not only
      occupied C_occ_k / C_occ_k_beta.
    status: met
  - id: ac-002
    summary: GpwDft::mp2 accepts multi-k UKS references
    type: runtime
    pass_when: |
      mp2() no longer throws "multi-k UMP2 is a follow-on" for UKS with
      n_kpts > 0; still requires mp2_ri, periodic cell, converged run,
      virtuals, and finite HOMO–LUMO gap per spin.
    status: met
  - id: ac-003
    summary: Kernel units (no full SCF)
    type: runtime
    pass_when: |
      test_cpu_gpw_mp2: multi-k UMP2 α=β recovers multi-k RMP2 OS/SS on
      synthetic complex tensors; nk=1 real-phase UMP2 recovers real
      mp2_uhf_energy.
    status: met
  - id: ac-004
    summary: Multi-k UMP2 regression on bct-Si4 2×1×1
    type: runtime
    pass_when: |
      regressions/gpw/periodic-ump2-kpoints builds under ATV_BUILD_GPW,
      never launches CP2K, asserts e_corr < 0, finite, OS+SS and
      e_total identities, n_kpts/ri, and dual-gate band vs RKS multi-k
      first-green (or CP2K pin when present). Wall ≤ 60 s.
    status: met
  - id: ac-005
    summary: Public API hygiene
    type: code
    pass_when: |
      No public *Config, no make_*/build_* factory for multi-k UMP2;
      parameters remain GpwDft members.
    status: met
  - id: ac-006
    summary: Offline oracle input shipped
    type: docs
    pass_when: |
      tools/cp2k_oracle/inputs/periodic-ump2-kpoints.inp + REFERENCE.md
      document twin settings for re-pin.
    status: met
out_of_scope:
  - UMP2 forces/stress
  - RI_AUX Gaussian tables
  - UGW/URPA
  - Live CP2K in ctest
---

# Acceptance — `periodic-ump2-kpoints`

**Done** means multi-k unrestricted GPW MP2 correlation energy is available on
`GpwDft::mp2` after a multi-k UKS `run()`, with kernel dual-gate identities and
a regression deck under `regressions/gpw/periodic-ump2-kpoints/`.
