---
spec: buckingham-pair-lammps
created: 2026-07-23
criteria:
  - id: ac-001
    summary: Offline LAMMPS pin documented
    type: code
    pass_when: |
      regressions/forcefield/buckingham-pair-lammps/REFERENCE.md states
      oracle: lammps and oracle_version: <string>, geometry, A/rho/C/rc,
      pair_modify shift yes, metal units, hard-coded E and F pins from
      offline lmp. No live LAMMPS in ctest/CMake.
    status: pass
  - id: ac-002
    summary: Public PairBuck API (no Config bag / factory)
    type: code
    pass_when: |
      cpu::PairBuck(A, rho, C, cutoff) exists under
      include/atomiverse/cpu/potentials/pair_buck.h. Parameters are ctor
      args only. No public *Config; no make_*/build_* factory. Inherits
      Pair/PairPotential; name() returns a stable string (e.g. "buck").
    status: pass
  - id: ac-003
    summary: Unit test analytic E and F
    type: runtime
    pass_when: |
      tests/cpu/potentials/test_pair_buck.cpp (ctest target
      test_cpu_pair_buck) asserts energy and force match analytic shifted
      Buckingham within 1e-12 on a fixed dimer (same formula as LAMMPS).
    status: pass
  - id: ac-004
    summary: E2E regression vs offline pin
    type: runtime
    pass_when: |
      regressions/forcefield/buckingham-pair-lammps exits 0:
      |E - E_pin| ≤ 1e-10 eV and |Fx - Fx_pin| ≤ 1e-10 eV/Å.
      Public headers + link atomiverse only. Runtime ≤ 60 s.
    status: pass
  - id: ac-005
    summary: MD-only build
    type: runtime
    pass_when: |
      Configure with ATV_BUILD_MD=ON ATV_BUILD_KSPACE=ON ATV_BUILD_TESTS=ON
      ATV_BUILD_REGRESSIONS=ON and HF/PW/GPW/CUDA/Torch all OFF. Build
      atomiverse + reg_buckingham-pair-lammps and unit test.
    status: pass
  - id: ac-006
    summary: Perf gate within 10x or exempt
    type: runtime
    pass_when: |
      Fair wall T_atv/T_lammps ≤ 10 recorded in REFERENCE.md, OR
      perf_exempt=true with explicit reason (e.g. dimer-only O(N²) foothold
      not a fair bulk twin). Green requires within_10x or exempt.
    status: pass
out_of_scope:
  - Fumi–Tosi multi-type
  - buck/coul/long
  - CUDA
  - Multi-element coefficient tables
---

# Acceptance — buckingham-pair-lammps

Done means a public shifted Buckingham pair Potential matches offline LAMMPS
`pair_style buck` energy and force on the ion-dimer pin under an MD-only
build, with dual-gate (correctness + perf) recorded in REFERENCE.md.
