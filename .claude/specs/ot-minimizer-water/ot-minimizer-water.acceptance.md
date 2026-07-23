---
spec: ot-minimizer-water
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Public use_ot member on GPWEvaluator (no Config bag)
    type: code
    pass_when: |
      GPWEvaluator has public bool use_ot (default false). build() sets
      GpwDft::use_ot from the member (env ATV_GPW_OT may still force on).
      No new public *Config bag; no make_*/build_* factory for OT.
    status: pending
  - id: ac-002
    summary: Unit tests for OT machinery and OT/DIIS energy parity
    type: runtime
    pass_when: |
      tests/cpu/gpw/test_cpu_gpw_ot (or equivalent) passes under ctest with
      ATV_BUILD_GPW=ON: (1) fixed-H OT reaches occupied band energy consistent
      with generalized diagonalization within 1e-6 Ha; (2) mini-H2 SCF with
      use_ot=true and use_ot=false agree on e_total within 1e-5 Ha when both
      converge; (3) GPWEvaluator with use_ot=true produces finite converged
      energy without requiring ATV_GPW_OT env.
    status: pending
  - id: ac-003
    summary: Offline CP2K OT pin documented
    type: code
    pass_when: |
      regressions/gpw/ot-minimizer-water/REFERENCE.md cites CP2K version,
      geometry, basis/PP, cutoff, OT settings, and hard-coded e_tot Ha pin.
      No live CP2K in ctest/CMake.
    status: pending
  - id: ac-004
    summary: E2E regression vs CP2K OT pin + OT/DIIS dual gate
    type: runtime
    pass_when: |
      regressions/gpw/ot-minimizer-water exits 0: |E_OT - E_CP2K| ≤ pin_tol
      (1.5e-3 Ha at 250 Ry per REFERENCE), and |E_OT - E_DIIS| ≤ 1e-4 Ha
      on the same geometry/settings (public API only). Runtime ≤ 60 s (OMP=1).
    status: pending
  - id: ac-005
    summary: Perf gate within 10x or exempt
    type: runtime
    pass_when: |
      Fair wall T_atv/T_cp2k ≤ 10 for OT twin recorded in REFERENCE.md, OR
      perf_exempt=true with explicit reason (e.g. no fair single-rank twin).
    status: pending
out_of_scope:
  - Sparse OT / linear scaling
  - Spin-polarized OT
  - k-points
  - MPI scaling
---

# Acceptance — ot-minimizer-water

Done means OT is a product-facing GPW switch with offline CP2K energy parity on
isolated water (MOLOPT), dual-checked against the in-tree DIIS path.
