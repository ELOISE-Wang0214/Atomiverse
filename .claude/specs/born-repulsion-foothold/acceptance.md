---
spec: born-repulsion-foothold
created: 2026-07-24
criteria:
  - id: ac-001
    summary: Offline LAMMPS pin documented
    type: code
    pass_when: |
      regressions/forcefield/born-repulsion-foothold/REFERENCE.md states
      oracle: lammps and oracle_version: <string>, geometry, A/rho/sigma/C/D/rc,
      pair_modify shift yes, metal units, hard-coded E and F pins from
      offline lmp. No live LAMMPS in ctest/CMake.
    status: pass
  - id: ac-002
    summary: Public PairBorn API (no Config bag / factory)
    type: code
    pass_when: |
      cpu::PairBorn(A, rho, sigma, C, D, cutoff) exists under
      include/atomiverse/cpu/potentials/pair_born.h. Parameters are ctor
      args only. No public *Config; no make_*/build_* factory. Inherits
      Pair; name() returns a stable string (e.g. "born").
    status: pass
  - id: ac-003
    summary: Unit test analytic E and F
    type: runtime
    pass_when: |
      tests/cpu/potentials/test_pair_born.cpp (ctest target
      test_cpu_pair_born) asserts energy and force match analytic shifted
      Born–Mayer–Huggins within 1e-12 on a fixed dimer (same formula as LAMMPS).
    status: pass
  - id: ac-004
    summary: E2E regression vs offline pin
    type: runtime
    pass_when: |
      regressions/forcefield/born-repulsion-foothold exits 0:
      |E - E_pin| ≤ 1e-10 eV and |Fx - Fx_pin| ≤ 1e-10 eV/Å.
      Public headers + link atomiverse only. Runtime ≤ 60 s.
    status: pass
  - id: ac-005
    summary: MD-only build
    type: runtime
    pass_when: |
      Configure with ATV_BUILD_MD=ON ATV_BUILD_KSPACE=ON ATV_BUILD_TESTS=ON
      ATV_BUILD_REGRESSIONS=ON and HF/PW/GPW/CUDA/Torch all OFF. Build
      atomiverse + reg_born-repulsion-foothold and unit test.
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
  - born/coul/long
  - Steele wall
  - CUDA
  - Multi-element coefficient tables
---

# Acceptance — born-repulsion-foothold

Done means a public shifted Born–Mayer–Huggins pair Potential matches offline
LAMMPS `pair_style born` energy and force on the ion-dimer pin under an MD-only
build, with dual-gate (correctness + perf) recorded in REFERENCE.md.
