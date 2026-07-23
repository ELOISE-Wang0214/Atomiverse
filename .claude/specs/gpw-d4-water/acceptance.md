---
spec: gpw-d4-water
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Offline pin documented (dftd4 / CP2K note)
    type: code
    pass_when: |
      regressions/gpw/gpw-d4-water/REFERENCE.md cites dftd4 version 3.7.0,
      PBE-D4 two-body (s9=0), geometry, cell, hard-coded E_disp pin, and
      notes CP2K 2026.1 Arrhenius build lacks DFTD4 lib. No live oracle in
      ctest/CMake.
    status: pass
  - id: ac-002
    summary: Public GpwDft::d4 / GPWEvaluator::d4 (no Config bag)
    type: code
    pass_when: |
      bool d4 on GpwDft and GPWEvaluator; mutually exclusive with d3/d3bj;
      e_disp on GpwResult; no public *Config; no make_*/build_* for D4.
    status: pass
  - id: ac-003
    summary: Unit test E_disp vs pin + optional FD if forces later
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_dftd4.cpp asserts periodic water-dimer
      E_disp within 1e-8 Ha of pin; nonperiodic path still matches
      molecular two-body within 1e-8 Ha on same geometry.
    status: pass
  - id: ac-004
    summary: E2E regression vs offline pin
    type: runtime
    pass_when: |
      regressions/gpw/gpw-d4-water exits 0: |E_disp - pin| ≤ 1e-8 Ha;
      e_total finite and includes e_disp (e_disp < 0). Public GPW headers
      + link atomiverse. Runtime ≤ 60 s.
    status: pass
  - id: ac-005
    summary: Build with ATV_BUILD_GPW=ON
    type: runtime
    pass_when: |
      Configure HF+PW+GPW+tests+regressions CUDA off; build reg_gpw-d4-water
      and test_cpu_gpw_dftd4.
    status: pass
  - id: ac-006
    summary: Perf gate within 10x or exempt
    type: runtime
    pass_when: |
      Fair wall recorded in REFERENCE.md, OR perf_exempt=true with reason
      (dispersion-only gate; full QS wall unfair — same as gpw-pbe-d3bj-ice).
    status: pass
out_of_scope:
  - ATM three-body
  - D4 forces/stress product gate
  - CP2K full e_tot
  - CUDA
---

# Acceptance — gpw-d4-water

Done means periodic PBE-D4 two-body E_disp on the water-dimer GPW path matches
the offline dftd4 3.7.0 pin to 1e-8 Ha, exposed as `GpwDft::d4` without public
Config bags, dual-gated in REFERENCE.md.
