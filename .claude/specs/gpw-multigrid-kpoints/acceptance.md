# Acceptance — gpw-multigrid-kpoints

All criteria verified 2026-07-31 on `wf/gpw-multigrid-kpoints` (see the
regression output and REFERENCE.md for the measured numbers).

- [x] `GpwDft::build` accepts `ngrids > 1 && rel_cutoff > 0` on the
      closed-shell multi-k path; UKS multi-k multigrid and `ngrids > 1`
      without a positive `rel_cutoff` throw with a how-to
      (`tests/cpu/gpw/test_cpu_gpw_solver.cpp` gate tests).
- [x] Unit identities (`tests/cpu/gpw/test_cpu_gpw_multigrid.cpp`,
      `GpwMultigridCells.*`): image-cell index bookkeeping; uniform-stack
      reduction to the Γ product path (bit-identical); per-cell integrate =
      exact transpose of per-cell collocate; D-stack density + phase-summed
      V_t match the independent `KpointBasisGrid` machinery at a complex k
      (Hermitian per-k grid KS); time-reversal D-stack folding exact.
- [x] Full-SCF bring-up identities, measured once and recorded in
      REFERENCE.md (not re-run per invocation — SCF cost vs the ≤ 60 s
      budget): Γ-limit 2.0e-14 Ha; TR full-vs-folded 2.8e-14 (95 Ry),
      2.0e-14 (60 Ry), 2.3e-14 (40 Ry).
- [x] Regression `regressions/gpw/gpw-multigrid-kpoints` (type: scientific,
      ~55 s at OMP 8):
  - M2 TR-folded multigrid SCF converges on the 40 Ry / REL 20 ladder
  - M3 multigrid vs single grid ≤ 1e-3 Ha (measured 1.1e-5)
  - M4 converged, finite, k payload shapes (150 Ry main)
  - M5 |E − CP2K 2026.1 pin| ≤ 1e-3 Ha (measured 4.7e-7); pin is an offline
    CP2K number (binary cp2k.psmp 2026.1 + committed twin input + logged
    output line), never seeded from ATV
- [x] Perf: T_atv/T_cp2k ≤ 10 on the NGRIDS 4 twin (measured < 1; REFERENCE).
- [x] Found-rot handled, not silenced:
  - periodic product screens lost integrated mass (0.09 e⁻ / 4.1e-2 Ha on Si
    bct) — fixed via `kProductMassMargin` in `internal/multigrid.cpp`;
    shipped `gpw-uks-multigrid` re-run green (|dE| vs CP2K now 1.26e-5 Ha)
  - k-path nonlocal vdW silently dropped E_c^nl — now rejected at build
  - `gpw-kpoints-si` self-golden A7 pin replaced by the genuine CP2K 2026.1
    number from its own committed twin input
