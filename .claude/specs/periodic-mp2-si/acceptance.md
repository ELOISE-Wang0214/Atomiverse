---
spec: periodic-mp2-si
created: 2026-07-25
criteria:
  - id: ac-001
    summary: Closed-shell doubles sum matches an explicit spin-orbital MP2 sum
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_mp2.cpp builds a small synthetic (ia|jb)
      tensor + eigenvalues and shows |E_closed_shell − E_spin_orbital| ≤ 1e-12
      Ha, and e_corr_os + e_corr_ss == e_corr (unscaled) to 1e-12.
    status: pass
  - id: ac-002
    summary: GPW MO-ERI machinery is consistent with the Poisson Hartree energy
    type: runtime
    pass_when: |
      On a toy periodic cell, the pair-integral route gives
      |(pp|pp) − 2·E_H[ψ_p²]| ≤ 1e-10 Ha where E_H comes from solve_poisson
      on the same grid, and (ia|jb) == (jb|ia) to 1e-12.
    status: pass
  - id: ac-003
    summary: Public surface is members + method, no bag, no factory
    type: code
    pass_when: |
      Public API is exactly GpwDft::mp2(atoms) + members
      (mp2_n_virt, mp2_os_scale, mp2_ss_scale) + gpw::Mp2Result value type in
      cpu/gpw/mp2.h. No public *Config, no make_*/build_* free function; ERI
      kernels live under cpu/gpw/internal/.
    status: pass
  - id: ac-004
    summary: Guards throw (unconverged SCF, UKS, vanishing gap)
    type: runtime
    pass_when: |
      Unit tests show GpwDft::mp2 throws std::logic_error before a converged
      run() and std::invalid_argument when configured for UKS.
    status: pass
  - id: ac-005
    summary: periodic-mp2-si regression green — Si2 MOLOPT+GTH Γ MP2
    type: runtime
    pass_when: |
      regressions/gpw/periodic-mp2-si builds via atv_add_regression under
      ATV_BUILD_GPW, uses public headers only, never launches CP2K, runs
      ≤ 60 s, and passes ALL of: (a) SCF converged, E_corr finite and < 0;
      (b) virtual-space monotonicity E_corr(all virt) ≤ E_corr(n_virt
      truncated) ≤ 0 (each added virtual can only lower the MP2 sum);
      (c) |e_corr_os + e_corr_ss − e_corr| ≤ 1e-12 Ha;
      (d) e_total == e_scf + e_corr to 1e-12;
      (e) E_corr within the frozen band pinned in REFERENCE.md (first-green
      ATV band until an offline CP2K twin re-pin; tol as stated there).
    status: pass
  - id: ac-006
    summary: REFERENCE.md documents the offline oracle contract
    type: docs
    pass_when: |
      regressions/gpw/periodic-mp2-si/REFERENCE.md records method, system,
      grid/SCF settings, the pin provenance (no CP2K binary on host — frozen
      first-green band + internal gates), and a CP2K input twin exists at
      tools/cp2k_oracle/inputs/periodic-mp2-si.inp for a reproducible re-pin.
    status: pass
---

# Acceptance — periodic-mp2-si

Green requires every criterion above pass. Perf: no CP2K binary on the
worker host (`cp2k`, `cp2k.psmp`, `cp2k.ssmp` absent; `module avail cp2k`
empty — same status as every prior wave) → `perf_exempt` with that reason,
wall time of the regression recorded in REFERENCE.md.
