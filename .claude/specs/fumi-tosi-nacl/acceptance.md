---
spec: fumi-tosi-nacl
created: 2026-07-24
criteria:
  - id: ac-001
    summary: Offline LAMMPS pin documented
    type: code
    pass_when: |
      regressions/forcefield/fumi-tosi-nacl/REFERENCE.md states
      oracle: lammps and oracle_version: <string>, geometry (4-ion NaCl
      square), multi-type Born A/rho/sigma/C/D table, coul/cut, hybrid/overlay
      twin, metal units, hard-coded E and F pins from offline lmp. No live
      LAMMPS in ctest/CMake.
    status: pass
  - id: ac-002
    summary: Public multi-type PairBorn API (no Config bag / factory)
    type: code
    pass_when: |
      cpu::PairBorn exists under include/atomiverse/cpu/potentials/pair_born.h
      with (A, rho, sigma, C, D, cutoff) single-type ctor and
      (vector<BornTypeCoeff>, cutoff) multi-type ctor. Parameters are ctor
      args / table only. No public *Config; no make_*/build_* factory.
      Inherits Pair/PairPotential; name() returns "born".
    status: pass
  - id: ac-003
    summary: Unit test analytic E and F
    type: runtime
    pass_when: |
      tests/cpu/potentials/test_pair_born.cpp (ctest target
      test_cpu_pair_born) asserts energy and force match analytic shifted
      Born–Mayer–Huggins within 1e-12 on a fixed dimer, and multi-type path
      selects the correct type-pair coefficients.
    status: pass
  - id: ac-004
    summary: E2E regression Born+Coul vs offline pin
    type: runtime
    pass_when: |
      regressions/forcefield/fumi-tosi-nacl exits 0:
      |E - E_pin| ≤ 1e-9 eV and |F - F_pin| ≤ 1e-8 eV/Å on the 4-ion NaCl
      square using PairBorn multi-type + PairCoulCut. Public headers + link
      atomiverse only. Runtime ≤ 60 s.
    status: pass
  - id: ac-005
    summary: MD-only build
    type: runtime
    pass_when: |
      Configure with ATV_BUILD_MD=ON ATV_BUILD_KSPACE=ON ATV_BUILD_TESTS=ON
      ATV_BUILD_REGRESSIONS=ON and HF/PW/GPW/CUDA/Torch all OFF. Build
      atomiverse + reg_fumi-tosi-nacl and unit test.
    status: pass
  - id: ac-006
    summary: Perf gate within 10x or exempt
    type: runtime
    pass_when: |
      Fair wall T_atv/T_lammps ≤ 10 recorded in REFERENCE.md, OR
      perf_exempt=true with explicit reason (e.g. 4-ion single-point O(N²)
      foothold not a fair bulk twin). Green requires within_10x or exempt.
    status: pass
out_of_scope:
  - born/coul/long + PME product path
  - CUDA
  - Steele wall
  - Public Config bags
---

# Acceptance — fumi-tosi-nacl

Done means a public multi-type shifted Born pair Potential, composed with
`PairCoulCut`, matches offline LAMMPS `hybrid/overlay born + coul/cut` energy
and force on the multi-type NaCl square under an MD-only build, with dual-gate
(correctness + perf) recorded in REFERENCE.md.
